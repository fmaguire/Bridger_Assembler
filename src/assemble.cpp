/*
 *
 * Copyright(c) 2012,2013. Zheng Chang.
 *
 * This file is part of program Bridger, which is used to de novo
 * assemble transcriptome from RNASeq data.
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
 *
 * Problems, suggestions, comments, etc ? Contact :
 *    Zheng Chang  <changzmaths@gmail.com>
 *
 * Created on June 25, 2012.
 */

#include "loadreads.h"
#include "common.h"
#include "splicing_graph.h"
#include "kmerhash.h"
#include "assemble.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <errno.h>
#include <ctime>

#include <boost/filesystem.hpp>


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
      reads_file = optarg;
      break;
    case 'o':
      out_dir = optarg;
      break;
    case 'h':
      g_help = true;
      break;
    case 'p':
      CPU = atoi(optarg);
      break;
    case OPT_KMERS:
      kmers_file = optarg;
      break;
    case OPT_DOUBLE_STRANDED_MODE:
      g_double_stranded_mode = true;
      break;
    case OPT_IS_PAIR_END:
      g_is_paired_end = true;
      break;
    case OPT_MIN_KMER_COVERAGE:
      g_min_kmer_coverage = atoi(optarg);
      break;
    case OPT_MIN_KMER_ENTROPY:
      g_min_kmer_entropy = atof(optarg);
      break;
    case OPT_MIN_SEED_COVERAGE:
      g_min_seed_coverage = atoi(optarg);
      break;
    case OPT_MIN_SEED_ENTROY:
      g_min_seed_entropy = atof(optarg);
      break;
    case OPT_MIN_JUNCTION_COVERAGE:
      g_min_junction_coverage = atoi(optarg);
      break;
    case OPT_MIN_RATIO_NON_ERROR:
      g_min_ratio_non_error = atof(optarg);
      break;
    case OPT_MIN_EXON_LENGTH:
      g_min_exon_length= atoi(optarg);
      break;
    case OPT_MIN_KMERS_PER_GRAPH:
      g_min_kmers_per_graph = atoi(optarg);
      break;
    case OPT_FR_STRAND:
      g_fr_strand = atoi(optarg);
      break;
    case OPT_PAIR_GAP_LENGTH:
      g_pair_gap_length = atoi(optarg);
      break;
    case OPT_DEBUG:
      g_debug = true;
      break;
    default:
      std::cout << usage();
      exit(1);
    }
  } while(next_option != -1);

  if (g_help) {
    std::cout << usage() ;
    exit(1);
  }

  if (reads_file == "") {
    std::cerr << "Error : --input option needs an argument!! " << std::endl;
    std::cout << usage() ;
    exit(1);
  }

  if (g_kmer_length > 32) {
    errAbort(const_cast<char *>("Length of kmer can not be excess 32!\n"));
  }

  if (g_min_anchor_length >= g_kmer_length) {
    g_min_anchor_length = g_kmer_length - 2;
  } 

  if (g_kmer_length != 25) {
    g_min_anchor_length = g_kmer_length - 3;
  } 
 
  if (g_fr_strand != 1 && g_fr_strand != 2) {
    errAbort(const_cast<char *>("--g_fr_strand can only be 1 or 2!\n"));
  }

  return 0;
}


std::string usage () {

  std::stringstream usage_info;
  usage_info
    << std::endl
    << "===============================================================================" << std::endl
    << " Usage: Assemble [--reads/--kmers] <filename>  [opts] " << std::endl
    << "===============================================================================" << std::endl
    << " **Required :" << std::endl
    << " --reads/-i <string>           " << ": the name of the file containing reads" << std::endl;
  //usage_info
  //  << " or " << std::endl
  //  << " --kmers <string>	       " << ": the name of the file containing kmers" << std::endl;
  usage_info
    << std::endl
    << " ** Optional :" << std::endl
    << " --kmer_length/-k <int>        " << ": length of kmer, default: 25." << std::endl
    << " --double_stranded_mode	       " << ": set it true if double stranded mode." << std::endl
    << " --fr_strand<int>	       " << ": strand specific protocol, default: 1 " << std::endl
    << "                                       ( 1 : fr-firststrand, e.g. dUTP, NSR, NNSR " << std::endl
    << "                                         2 : fr-secondstrand, e.g. Strandard SOLID ) " << std::endl
    << " --paired_end		       " << ": set it true if paired reads." << std::endl
    << " --min_seed_coverage <int>     " << ": minimum coverage of seed kmer, default: 2." << std::endl
    << " --min_seed_entropy <float>    " << ": minimum entropy of seed kmer, default: 1.5" << std::endl
    << " --min_kmer_coverage <int>     " << ": minimum coverage of kmer used to extend, default: 1." << std::endl
    << " --min_kmer_entropy <float>    " << ": minimum entroy of kmer used to extend, default: 0.0" << std::endl
    << " --min_junction_coverage <int> " << ": minimum of the coverage of a junction, default: 2." << std::endl
    << " --min_ratio_non_error <float> " << ": min ratio for low/high alternative extension that is "<< std::endl 
    << "                                        not an error, default: 0.05." << std::endl
    << " --pair_gap_length	       " << ": gap length of paired reads, default: 200." << std::endl
    << " --out_dir/-o <string>         " << ": name of directory for output, default : ./RawGraphs " << std::endl
    << " --help/-h                     " << ": display the help information."<< std::endl
    << std::endl;
  usage_info
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

  if (out_dir != "") {
    //int retcode = mkpath(out_dir.c_str(), 0777);
    bool retcode = boost::filesystem::create_directories(out_dir);
    //if (retcode == -1) {
    if (retcode == false) {
      errAbort(const_cast<char *>("Cannot create directory %s !\n"), out_dir.c_str());
    }
  }
  
  // load data
  std::vector<std::string> data;
  if (reads_file != "")
    load_data(data, reads_file);
  
  KmerMap<info_with_read_set_t> kmer_map(g_kmer_length, g_double_stranded_mode);
  
  assembler(kmer_map, data);

  time_t end = time(NULL);
  std::cerr << "\nDone. " << rg_index << " splcicing grphs are built. Elapsed time: " << (end-begin) <<" s"<< std::endl;  // calculate used time

  data.clear();
  std::vector<std::string> s;
  swap(data,s);
  time_t end1 = time(NULL);
  std::cerr << "Clear memory. Elapsed time: " << (end1-end) <<" s"<< std::endl;
 
  return 0;
}

