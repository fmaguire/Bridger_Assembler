#ifndef GET_READS_H
#define GET_READS_H

// get_reads_from_graph.h

#include <stdlib.h>
#include <string>
#include <sstream>
#include "common.h"

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
  {"min_transcript_length",       required_argument,      0,      OPT_MIN_TRANSCRIPT_LENGTH},
  {"double_stranded_mode",        no_argument,            0,      OPT_DOUBLE_STRANDED_MODE},
  {"pair_end",                    no_argument,            0,      OPT_IS_PAIR_END},
  {"input",                       required_argument,      0,      'i'},
  {"output",                      required_argument,      0,      'o'},
  {"debug",			  no_argument,            0,      OPT_DEBUG},
//...
  {0,0,0,0} // terminator

};
                                      
int parse_options(int argc, char* argv[]);
std::string usage();


#endif
