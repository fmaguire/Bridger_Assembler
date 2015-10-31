#ifndef ASSEMBLE_H
#define ASSEMBLE_H
/*
 * assemble.h
 *
 */
#include <boost/serialization/set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "common.h"
#include "kmerhash.h"
#include "splicing_graph.h"

#ifdef WIN32
#define __STDC__ 1
#include <Windows.h>
#include "windows/getopt.h"
#else
#include <getopt.h>
#endif


static const char *short_options = "k:i:o:p:h";

static struct option long_options[] = {
  // general options
  {"kmer_length",                 required_argument,      0,      'k'},
  {"help",                        no_argument,            0,      'h'},
  {"CPU",			  required_argument,      0,      'p'},
  // assemble option
  {"min_kmer_coverage",           required_argument,      0,      OPT_MIN_KMER_COVERAGE},
  {"min_kmer_entropy",            required_argument,      0,      OPT_MIN_KMER_ENTROPY},
  {"min_seed_coverage",           required_argument,      0,      OPT_MIN_SEED_COVERAGE},
  {"min_seed_entropy",            required_argument,      0,      OPT_MIN_SEED_ENTROY },
  {"min_junction_coverage",       required_argument,      0,      OPT_MIN_JUNCTION_COVERAGE},
  {"min_ratio_non_error",         required_argument,      0,      OPT_MIN_RATIO_NON_ERROR},
  {"min_exon_length",             required_argument,      0,      OPT_MIN_EXON_LENGTH},
  {"min_kmers_per_graph",         required_argument,      0,      OPT_MIN_KMERS_PER_GRAPH},
  {"reads",                       required_argument,      0,      'i'},
  {"kmers",			  required_argument,	  0,	  OPT_KMERS},
  {"double_stranded_mode",	  no_argument,		  0,	  OPT_DOUBLE_STRANDED_MODE},
  {"pair_gap_length",             required_argument,      0,      OPT_PAIR_GAP_LENGTH},
  {"fr_strand",			  required_argument,      0,      OPT_FR_STRAND},
  {"out_dir",                     required_argument,      0,      'o'},
  {"pair_end",			  no_argument,            0,	  OPT_IS_PAIR_END},
  {"debug",			  no_argument,            0,      OPT_DEBUG},
//...
  {0,0,0,0} // terminator

};

int parse_options(int argc, char* argv[]);
std::string usage();

static size_t rg_index = 0;

std::string base_name() {

    std::stringstream idx ;
    idx << "comp" << rg_index ;
    rg_index++;

    return idx.str();
}


template <typename T>
void assembler(KmerMap<T>& kmer_map, std::vector<std::string>& data) {

  // build the hash 
  if (!data.empty()) {
    kmer_map.get_hash(data);
  } else if (kmers_file != "") {
    //kmer_map.get_hash_from_kmers(kmers_file, g_kmer_length);
  } else { 
    errAbort(const_cast<char *>("Building kmer hash failed."));
  }
 
  // prune errors
  kmer_map.remove_erroneous_kmers(g_min_ratio_non_error);
  kmer_map.prune_hash(g_min_kmer_coverage, g_min_kmer_entropy);

  // get seed kmers
  std::vector<kmer_int_type_t> seeds_sorted;
  kmer_map.get_seed_sort_descending_counts(seeds_sorted);
  //seeds_sorted = kmer_map.get_seed_sort_descending_counts();

  // assemble graphs one by one
  if (seeds_sorted.empty()) 
    errAbort(const_cast<char *>("No seeds available!\n"));

  const std::string & raw_graph_list = out_dir + "/" + rg_list;
  std::fstream out_rg;
  out_rg.open(raw_graph_list.c_str(), std::fstream::out);
  if (!out_rg.is_open()) {
    errAbort(const_cast<char *>("File %s can't be opened.\n"), raw_graph_list.c_str());
  }

  for (unsigned int i = 0; i < seeds_sorted.size(); ++i) {

    if (kmer_map.exists(seeds_sorted[i])) {

      SplicingGraph<T> splicing_graph;
      if (splicing_graph.build(kmer_map, seeds_sorted[i], data)) {

        const std::string & basename = base_name();
        // save each graph as one individual file
        const std::string name = out_dir + "/" + basename + ".rg";
	splicing_graph.save(name);

        // record the basename of graph
        out_rg << basename << std::endl;

        // output the graph in a readable format for debug
        if (g_debug)  
          splicing_graph.save_debug(name + ".debug"); 

        if (rg_index % 2000 == 0)
          std::cerr << rg_index << " graphs have been built." << std::endl;

        // remove kmers used in current graph from hash table
        splicing_graph.clear(kmer_map);
      }
    }
  } // for

  out_rg.close();
}


#endif
