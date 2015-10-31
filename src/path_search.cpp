/*
 *
 * Copyright(c) 2012. Zheng Chang, 
 * The University of Georgia, USA.
 *                                                                                      
 * This file is part of program Bridger, which is used to de novo 
 * assemble transcriptome from RNASeq data.                                         
 *                                                                                      
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * Version : 0.1
 *
 * Problems, suggestions, comments, etc ? Contact :
 *    Zheng Chang  <changzmaths@gmail.com>
 *
 * Created on June 25, 2012
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <time.h>


#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

#include "common.h"
#include "loadreads.h"
#include "kmerhash.h"
#include "splicing_graph.h"
#include "compatibility_graph.h"
#include "transitive_closure.h"
#include "reachability_bp_graph.h"
#include "matching_merge.h"
#include "path_search.h"


int parse_options(int argc, char* argv[]) {
  int option_index = 0;
  int next_option;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
      case -1:     /* Done with options. */
        break;
      case 'k':
        g_kmer_length = atoi(optarg);
        break;
      case 'i':
        rg_file = optarg;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'h':
        g_help = true;
        break;
      case OPT_MIN_TRANSCRIPT_LENGTH:
        g_min_transcript_length = atoi(optarg);
        break;
      case OPT_PAIR_GAP_LENGTH:
        g_pair_gap_length = atoi(optarg);
        //g_max_pair_gap_length = 2* g_pair_gap_length;
        break;
      case OPT_MAX_PAIR_GAP_LENGTH:
	g_max_pair_gap_length = atoi(optarg);
        break;
      case OPT_DOUBLE_STRANDED_MODE:
        g_double_stranded_mode = true;
        break;
      case OPT_IS_PAIR_END:
        g_is_paired_end = true;
        break;
      case OPT_DEBUG:
	g_debug = true;
	break;
      default:
        std::cout << usage();
        exit(1);
    }
  } while(next_option != -1);
  if (rg_file == "") {
    std::cerr << "[Error] option '--input/-i' requires a argument!" << std::endl;
    std::cout << usage() ;
    exit(1);
  }
  if (g_help) {
    std::cout << usage() ;
    exit(0);
  }
  return 0;
}



std::string usage () {
  std::stringstream usage_info;
  usage_info
      << std::endl
      << "===============================================================================" << std::endl
      << " Usage: PathSearch --input/-i <raw_graph>  [opts] " << std::endl
      << "===============================================================================" << std::endl
      << " **Required :" << std::endl
      << " --input/-i <str>              " << ": basename of raw graph file in RawGraphs subdirectory" << std::endl;
  usage_info
      << std::endl
      << " ** Optional :" << std::endl
      << " --kmer_length/-k <int>        " << ": length of kmer, default: 25" << std::endl
      << " --double_stranded_mode        " << ": use it if double stranded mode" << std::endl
      << " --pair_end                    " << ": use it if paired reads" << std::endl
      << " --min_transcript_length       " << ": min length of transcript" << std::endl
      << " --pair_gap_length             " << ": gap length of paired reads, default : 200" << std::endl
      << " --pair_max_gap_length         " << ": max gap length of paired reads, default : 500" << std::endl
      << " --output/-o <str>             " << ": name of output file, default print on screen." << std::endl
      << " --help/-h                     " << ": display the help information."<< std::endl
      << std::endl;
  usage_info
      << " ** Note:" <<std::endl
      << " Plese ensure the a directory named RawGraphs is in the current path; otherwise, you will " << std::endl
      << " see error information like 'terminate called after throwing an instance of boost::archive::archive_exception"<< std::endl
      << "================================================================================" << std::endl
      << std::endl;
  return (usage_info.str());
}


int main(int argc, char* argv[]) {

  // process command line arguments
  int parse_ret = parse_options(argc, argv);
  if (parse_ret) 
    return parse_ret;
  
  time_t begin = time(NULL);

  SplicingGraph<info_with_read_set_t> splicing_graph;
  std::string infile = "RawGraphs/" + rg_file + ".rg";
  
  // check existence of infile
  namespace fs = boost::filesystem;
  if (!fs::exists("RawGraphs/")) {
    std::cerr << "File " << infile << " does Not exists! " << std::endl 
       <<"Make sure the directory 'RawGraphs',which contains all graphs built by Assemble, is in your current path. "<< std::endl;
    std::cerr << usage() << std::endl;
     return -1;
  }
  splicing_graph.load(infile);

  if (g_debug)
    std::cerr << splicing_graph.describe_graph();
  
  splicing_graph.compact_graph();

  // get a kmer_map based on reads of graph
  KmerMap<info_with_read_set_t> kmer_map(g_kmer_length, g_double_stranded_mode);

  kmer_map.get_hash(splicing_graph.reads_, false);
  
  //compute coverage of each edge and each node
  //std::map<pair_t, double> edge_cov_map; 
  //splicing_graph.get_coverage_of_edges(edge_cov_map);
  //std::map<int, double> node_cov_map;
  //splicing_graph.get_coverage_of_nodes(node_cov_map);


  // compute coverage of each edge and each node
  //std::map<pair_t, double> edge_cov_map;
  // trim graph by deleting small branches, which are derived from sequencing error
  splicing_graph.trim_graph(kmer_map, splicing_graph.reads_);  


  splicing_graph.refine_tips(kmer_map, splicing_graph.reads_);

  if (splicing_graph.get_size() == 1) {  // just one contig, trivial case
    splicing_graph.refine_tips(kmer_map, splicing_graph.reads_);
    std::cout <<  ">" << rg_file <<"_seq0 len="
      << splicing_graph.node_set_[0].sequence.length() << " path=[0:"
      << splicing_graph.node_set_[0].sequence.length() << "]" << std::endl
      << add_fasta_seq_line_breaks(splicing_graph.node_set_[0].sequence, 60) << std::endl;
    return 0;
  }


  //if (g_debug)
  //  std::cerr << splicing_graph.describe_graph();

  std::map<pair_t, double> edge_cov_map; 
  splicing_graph.get_coverage_of_edges(edge_cov_map);

  // get read set of each node by mapping reads to graph
  std::vector<std::set<size_t> > reads_mapped_to_node;
  splicing_graph.map_reads_to_graph(kmer_map, reads_mapped_to_node);

  // splicing graph -> compatibility graph
  directed_acyclic_graph_t compatibility_graph;
  std::vector<pair_t> splicing_graph_edges;
  create_compatibility_graph(compatibility_graph, splicing_graph, splicing_graph_edges, edge_cov_map, kmer_map, reads_mapped_to_node); 

  // add terminal nodes
  pair<directed_acyclic_graph_node_t, directed_acyclic_graph_node_t>
  terminal = add_terminal_nodes(compatibility_graph);
  directed_acyclic_graph_node_t source = terminal.first;
  directed_acyclic_graph_node_t sink = terminal.second;
  PercentSplicedInMap
    psi_map = boost::get(vertex_percent_spliced_in_t(), compatibility_graph);
  psi_map[source] = 1.0f;
  psi_map[sink] = 1.0f;

  // transitive closure
  boost::adjacency_list<> transitive_closure;
  get_transitive_closure(compatibility_graph, transitive_closure);

  // compatibility graph -> reachability bipartite graph
  reachability_bp_graph_t reachability_bp_graph;
  create_reachability_bp_graph(compatibility_graph, reachability_bp_graph, transitive_closure);

  if (g_debug)
    std::cerr << "closure edges :" << num_edges(transitive_closure) << std::endl;  //debug

  // add weights to reachability_bp_graph
  reachability_bp_graph_t::UEdgeMap<long long> weights(reachability_bp_graph);
  add_weights_to_reachability_bp_graph(reachability_bp_graph, psi_map, weights);

  // matching
  typedef lemon::MinCostMaxBipartiteMatching
    <reachability_bp_graph_t,reachability_bp_graph_t::UEdgeMap<long long> > Matcher;
  Matcher matcher(reachability_bp_graph, weights);
  matcher.run();

  // get paths from matching 
  std::vector<std::vector<directed_acyclic_graph_node_t> > chains;
  make_chains_from_matching<Matcher>(reachability_bp_graph, matcher, chains);

  std::vector<std::vector<directed_acyclic_graph_node_t> > paths;
  extend_chains_to_paths(compatibility_graph, chains, transitive_closure, source, sink, paths, splicing_graph_edges, splicing_graph, edge_cov_map);

  std::vector<std::vector<directed_acyclic_graph_node_t> > paths_checked;
  check_paths(splicing_graph, kmer_map, splicing_graph_edges, edge_cov_map, reads_mapped_to_node, paths, paths_checked);
 
  // get isoforms  
  get_isoforms_from_paths(splicing_graph, paths_checked, splicing_graph_edges);

  time_t end = time(NULL);
  std::cerr << "Paths search done! (" << (end-begin) <<" seconds)"<< endl;  

  return 0;
}





