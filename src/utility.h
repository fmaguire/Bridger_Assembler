/*
       utility.h
*/

// This file was modified from a file called sequenceUtil.hpp 
// in Inchworm modules of Trintiy.
// The original copyright info is listed below
// Copyright (c) 2010, The Broad Institute, Inc. 
// Distributed under the  Distributed under Trinity Software LICENSE
// (See accompanying file LICENSE)


#ifndef UTILITY_H
#define UTILITY_H


#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <cstdarg>

// Compatibility of __attribute__ with non-GNU
#ifndef __GNUC__
#define __attribute__(x) 
#endif

 
// code each kmer as a 64-bit unsigned integer
typedef unsigned long long kmer_int_type_t;

// check a kmer contains characters that not gatc
bool contains_non_gatc(const std::string& kmer);
bool contains_non_gatc(const std::string& kmer,unsigned int kmer_length);

// get the reverse complementary std::string of a kmer
std::string revcomp(const std::string& kmer);
kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length);

char int_to_base(int baseval); // 0 1 2 3 => G A T C
int base_to_int(char nucleotide); // (GATC) = {0 1 2 3}, others = -1

// convert a kmer into a 64-bit interger
kmer_int_type_t kmer_to_intval(const std::string& kmer); // must be less than 32 bases 
kmer_int_type_t kmer_to_intval(const std::string& kmer,unsigned int kmer_length);

// convert 64-bit integer to kmer
std::string intval_to_kmer (kmer_int_type_t intval, unsigned int kmer_length);

// store a kmer as itself or its reverse complement
kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length);

// compute entropy of a kmer
float compute_entropy(const std::string& kmer);
float compute_entropy(kmer_int_type_t kmer, unsigned int kmer_length);

std::string replace_nonGATC_chars_with_A(std::string& input_seq);

std::vector<std::string> potential_variations(const std::string& kmer);


// Print progress message 
void progress(char *format, ...)
     __attribute__((format(printf, 1, 2)));
// Print error message but do not exit
void err(char *format, ...)
     __attribute__((format(printf, 1, 2)));
// Print error message to stderr and exit
void errAbort(char *format, ...)
     __attribute__((noreturn, format(printf, 1, 2)));


#endif
