/*
 * get_reads_from_graph.cpp
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <time.h>

#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include "loadreads.h"
#include "common.h"
#include "kmerhash.h"
#include "splicing_graph.h"
#include "get_reads_from_graph.h"



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
      case 'p':
	CPU =  atoi(optarg);
	break;
      case OPT_MIN_TRANSCRIPT_LENGTH:
        g_min_transcript_length = atoi(optarg);
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
    std::cerr << "Error : --input must have a argument!! " << std::endl;
    std::cout << usage() ;
    exit(1);
  }
  if (g_help) {
    std::cout << usage() ;
    exit(1);
  }
  return 0;
}



std::string usage () {
  std::stringstream usage_info;
  usage_info
      << std::endl
      << "===============================================================================" << std::endl
      << " Usage: get_read_from_Graph --input/-i <graph_file>  [opts] " << std::endl
      << "===============================================================================" << std::endl
      << " **Required :" << std::endl
      << " --input/-i <str>              " << ": file of raw graph" << std::endl;
  usage_info
      << std::endl
      << " ** Optional :" << std::endl
      << " --kmer_length/-k <int>        " << ": the length of kmer, default: 25" << std::endl
      << " --double_stranded_mode        " << ": please use it if double stranded mode" << std::endl
      << " --pair_end                    " << ": please use it if paired reads" << std::endl
      << " --output/-o <str>             " << ": the output file, default : paths.fa " << std::endl
      << " --help/-h                     " << ": display the help information"<< std::endl
      << std::endl;
  usage_info
      << "================================================================================" << std::endl
      << std::endl;
  return (usage_info.str());
}


int main(int argc, char* argv[]) {

  // process command line arguments
  int parse_ret = parse_options(argc, argv);
  if (parse_ret) {
    return parse_ret;
  }

  SplicingGraph<info_with_read_set_t> splicing_graph;
  //string infile = "RawGraphs/" + rg_file + ".rg";
  splicing_graph.load(rg_file);

  // output the reads that belong to current graph in fasta format
  splicing_graph.get_reads_from_graph();

  return 0;
}
