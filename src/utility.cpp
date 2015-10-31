/*
   utinity.cpp
*/

// This file was modified from a file called sequenceUtil.cpp 
// in Inchworm modules of Trintiy.
// The original copyright info is listed below
// Copyright (c) 2010, The Broad Institute, Inc. 
// Distributed under the  Distributed under Trinity Software LICENSE
// (See accompanying file LICENSE)


#include <sstream>
#include "utility.h"
#include <stdlib.h>
#include <math.h>
#include <libgen.h>
//#include <sstream>

char _int_to_base [4] = {'G', 'A', 'T', 'C'};
unsigned char _base_to_int [256] = {
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0-19
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20-39
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40-59
  255, 255, 255, 255, 255,   1, 255,   3, 255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, //  60-79
  255, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   1, 255,   3, //  80-99
  255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   2, 255, 255, 255, // 100-119
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120-139
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140-159
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160-179
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180-209
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200-219
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220-239
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240-255
};

bool contains_non_gatc (const std::string& kmer) {
  for (unsigned int i = 0; i < kmer.size(); ++i) {
    unsigned char c = kmer[i];
    /*
    if (! (c == 'g' || c == 'G'
           || c == 'a' || c == 'A'
           || c == 't' || c == 'T'
           || c == 'c' || c == 'C')
       ) {
      return(true);
    }
   */
   if (_base_to_int[c] > 3)
     return(true);
  }
  return(false);
}

bool contains_non_gatc (const std::string& kmer, unsigned int kmer_length) {
  for (unsigned int i = 0; i < kmer_length; ++i) {
    unsigned char c = kmer[i];
    /*
    if (! (c == 'g' || c == 'G'
           || c == 'a' || c == 'A'
           || c == 't' || c == 'T'
           || c == 'c' || c == 'C')
       ) {
      return(true);
    }
    */
    if (_base_to_int[c] > 3)
     return(true);
  }
  return(false);
}

std::string revcomp (const std::string& kmer) {
  std::string revstring;
  for (int i = kmer.size() -1; i >= 0; --i) {
    char c = kmer[i];
    char revchar;
    switch (c) {
      case 'g':
        revchar = 'c';
        break;
      case 'G':
        revchar = 'C';
        break;
      case 'a':
        revchar = 't';
        break;
      case 'A':
        revchar = 'T';
        break;
      case 't':
        revchar = 'a';
        break;
      case 'T':
        revchar = 'A';
        break;
      case 'c':
        revchar = 'g';
        break;
      case 'C':
        revchar = 'G';
        break;
      default:
        revchar = 'N';
    }
    revstring += revchar;
  }
  return (revstring);
}

kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length) {
  kmer_int_type_t rev_kmer = 0;
  kmer = ~kmer;
  for (unsigned int i = 0; i < kmer_length; i++) {
    int base = kmer & 3;
    rev_kmer = rev_kmer << 2;
    rev_kmer += base;
    kmer = kmer >> 2;
  }
  return rev_kmer;
}

char int_to_base (int baseval) {
  if (baseval < 0 || baseval > 3) {
    std::cerr << "Error, baseval out of range 0-3" << std::endl;
    exit(1);
  }
  return(_int_to_base[baseval]);
}

int base_to_int (char nucleotide) {
  switch (nucleotide) {
    case 'G':
    case 'g':
      return(0);
    case 'A':
    case 'a':
      return(1);
    case 'T':
    case 't':
      return(2);
    case 'C':
    case 'c':
      return(3);
    default:
      return(-1);
  }
}


kmer_int_type_t kmer_to_intval (const std::string& kmer) {
  kmer_int_type_t kmer_val = 0;
  for (unsigned int i = 0; i < kmer.length(); ++i) {
    char c = kmer[i];
    int val = base_to_int(c);
    kmer_val = kmer_val << 2;
    kmer_val |= val;
  }
  return (kmer_val);
}

kmer_int_type_t kmer_to_intval(const std::string& kmer, unsigned int kmer_length) {
  kmer_int_type_t kmer_val = 0;
  for (unsigned int i = 0; i < kmer_length; ++i) {
    char c = kmer[i];
    int val = base_to_int(c);
    kmer_val = kmer_val << 2;
    kmer_val |= val;
  }
  return (kmer_val);
}

std::string intval_to_kmer (kmer_int_type_t intval, unsigned int kmer_length) {
  // std::string kmer = "";
  std::string kmer(kmer_length, ' ');
  for (unsigned int i = 1; i <= kmer_length; ++i) {
    int base_num = intval & 3ll;
    // char base = int_to_base(base_num);
    kmer[kmer_length-i] = _int_to_base[base_num];
    intval = intval >> 2;
    // kmer = base + kmer;
  }
  return (kmer);
}

float compute_entropy(const std::string& kmer) {
  std::map<char,int> char_map;
  for (unsigned int i = 0; i < kmer.length(); ++i) {
    char c = kmer[i];
    char_map[c]++;
  }
  float entropy = 0;
  char nucs[] = { 'G', 'A', 'T', 'C' };
  for (unsigned int i = 0; i < 4; ++i) {
    char nuc = nucs[i];
    int count = char_map[nuc];
    float prob = (float)count / kmer.length();
    if (prob > 0) {
      float val = prob * log(1.f/prob)/log(2.f);
      entropy += val;
    }
  }
  return(entropy);
}

float compute_entropy(kmer_int_type_t kmer, unsigned int kmer_length) {
  char counts[] = { 0, 0, 0, 0 };
  for (unsigned int i = 0; i < kmer_length; ++i) {
    int c = kmer & 3;
    kmer = kmer >> 2;
    counts[c]++;
  }
  float entropy = 0;
  for (unsigned int i = 0; i < 4; i++) {
    float prob = (float)counts[i] / kmer_length;
    if (prob > 0) {
      float val = prob * log(1.f/prob)/log(2.f);
      entropy += val;
    }
  }
  return(entropy);
}

// store a kmer as itself or its reverse complement, whichever
// corresponds to a larger 64-bit int.
kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length) {
  kmer_int_type_t rev_kmer = revcomp_val(kmer_val, kmer_length);
  if (rev_kmer > kmer_val)
    kmer_val = rev_kmer;
  return(kmer_val);
}

/*
std::string replace_nonGATC_chars_with_A (std::string& str) {

  std::stringstream newStr;

  for (unsigned int i = 0; i < str.size(); ++i) {
    char c = str[i];
    if (_base_to_int[c] > 3)
       c = 'A';
    newStr << c;
  }

  return(newStr.str());
}
*/

std::vector<std::string> potential_variations(const std::string& kmer) {
  std::vector<std::string> variations;
  for (unsigned int i = 0 ; i < kmer.length(); ++i) {
    std::string variation = kmer;
    char nucs[] = { 'G', 'A', 'T', 'C' };
    for (int j = 0; j < 4; ++j) {
      if (kmer[i] != nucs[j]) {
        variation[i] = nucs[j];
 	variations.push_back(variation);
      }
    }
  }
  return variations;
}


// Print progress message
void progress(char *format, ...) {
  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
}

// Print error message but do not exit
void err(char *format, ...) {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "[Warning] ");
  vfprintf(stderr, format, args);
  va_end(args);
}

// Print error message and exit
void errAbort(char *format, ...) {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "[Error] ");
  vfprintf(stderr, format, args);
  va_end(args);
  exit(1);
}

