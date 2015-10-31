// This file was modified from file KmerCounter.cpp
// in Inchworm modules of Trintiy.
// The original copyright info is listed below
// Copyright (c) 2010, The Broad Institute, Inc. 
// Distributed under the  Distributed under Trinity Software LICENSE
// (See accompanying file LICENSE)

#ifndef KMERHASH_H
#define KMERHASH_H

/*
 *  kmer_hash.h
 */

#include "utility.h"
#include "common.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <boost/unordered_map.hpp>

using namespace std;

typedef pair<kmer_int_type_t, size_t> kmer_occurence_pair_t;

class info_base_t {
public:
  void virtual add(size_t i, size_t j) = 0;
  size_t virtual get_count() const = 0;
  void virtual clear() = 0;
  virtual ~info_base_t() {}
};

class info_t : info_base_t {
public:

  info_t() : count(0) {}

  inline void add(size_t i, size_t j) {
    count++;
  }

  inline void add(size_t i, size_t j, size_t interval) {
    count++;
  }

  inline size_t get_count() const {return count;}

  inline void clear() {count = 0;}

private:
  size_t count;
};

class info_with_read_set_t : info_base_t {
public:

  info_with_read_set_t() : count(0) {}

  void add(size_t i, size_t j) {
    count++;
    if ( j % g_interval == 0)
      readset.push_back(i);  // record read index
  }

  void add(size_t i, size_t j, size_t interval) {
    count++;
    if ( j % interval == 0)
      readset.push_back(i);  // record read index
  }

  inline std::vector<size_t> get_readset() {
    return readset;
  }

  inline size_t get_count() const {
    return count;
  } 

  inline void clear() {
    count = 0; 
    //readset.clear();
  }

private:

  size_t count;
  std::vector<size_t> readset;
};



template<typename T>
class KmerMap {
private:

  typedef typename boost::unordered_map<kmer_int_type_t, T> kmer_map_base_t;
  typedef typename boost::unordered_map<kmer_int_type_t, T>::iterator kmer_map_base_iterator_t;
  typedef typename boost::unordered_map<kmer_int_type_t, T>::const_iterator kmer_map_const_iterator_t;

  /* sorter class */
  class kmer_sorter_by_count_desc_t {
  public:
    kmer_sorter_by_count_desc_t(KmerMap<T>& k) : kmer_map_(k) {};

    bool operator() (const kmer_int_type_t& i,
      const kmer_int_type_t& j) {
        return ( (kmer_map_[i].get_count() > kmer_map_[j].get_count())
          || (kmer_map_[i].get_count() == kmer_map_[j].get_count() && i > j) );
    }

    bool operator() (const kmer_occurence_pair_t& i,
      const kmer_occurence_pair_t& j) {
        return ( (i.second > j.second)
          || (i.second == j.second && i.first > j.first) );
    }

    bool operator() (const kmer_map_base_iterator_t& i,
      const kmer_map_base_iterator_t& j) {
	size_t count_i = i->second.get_count();
	size_t count_j = j->second.get_count();
        return ( (count_i > count_j)
		|| ( count_i == count_j && i->first > j->first));
    }

    bool operator() (const std::string& i, const std::string& j) {
      kmer_int_type_t val_i = kmer_to_intval(i);
      kmer_int_type_t val_j = kmer_to_intval(j);
      return ( (kmer_map_[val_i].get_count() > kmer_map_[val_j].get_count())
          || (kmer_map_[val_i].get_count() == kmer_map_[val_j].get_count() && val_i > val_j) );
    }
    
  private:
    KmerMap& kmer_map_;
  };

  kmer_map_base_t kmer_map;
  bool DS_MODE;
  int kmer_length;  

public:
  // constructor
  KmerMap () { }
  KmerMap (int kmer_length = g_kmer_length, bool is_ds = false) {
    this->kmer_length = kmer_length;
    this->DS_MODE = is_ds;
  }
  // define []
  T & operator[](kmer_int_type_t kmer_val) {
    if (DS_MODE) {
      kmer_val = get_DS_kmer_val(kmer_val, kmer_length);
    }
    return kmer_map[kmer_val];
  }

  size_t get_size() {
    return kmer_map.size();
  }

  bool empty() {
    return kmer_map.empty();
  }

  // find a kmer from kmer 
  kmer_map_base_iterator_t find_kmer(kmer_int_type_t kmer_val) {

    if (DS_MODE) 
      kmer_val = get_DS_kmer_val(kmer_val, kmer_length);
   
    return kmer_map.find(kmer_val);
  }

  bool reuse(const kmer_int_type_t kmer_val) {
    kmer_map_base_iterator_t it = find_kmer(kmer_val);
    return (it != kmer_map.end());
  }
  
  // get abundance of one kmer
  size_t get_kmer_count(kmer_int_type_t kmer_val) {

    kmer_map_base_iterator_t it = find_kmer(kmer_val);

    if (it != kmer_map.end())
      return((it->second).get_count());
    else
      return 0;
  }

  size_t get_kmer_count(const string& kmer) {

    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    return (get_kmer_count(kmer_val));
  }

  // check a kmer exist or not
  bool exists(const kmer_int_type_t kmer_val) {

    return (get_kmer_count(kmer_val) > 0);
  }

  bool exists(const std::string& kmer) {

    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    return (exists(kmer_val));
  }

  // get kmer hash from reads
  void get_hash(std::vector<std::string>& data, bool flag = true) {

    // warning : calculating read_length and data_size is not a wise choice !
    size_t data_size = data.size();
    if (flag) 
      std::cerr << "Beginning kmer hash ..." << std::endl;
    time_t beg = time(NULL);

    if (data.empty())
      return;

    if ( (int)data[0].length() < 2 * g_kmer_length) { // each read id occur at least twice
       g_interval = data[0].length() - g_kmer_length;
    } else {
       g_interval = g_kmer_length;
    }

    for (size_t i = 0; i < data_size; ++i) {

      const std::string& sequence = data[i];
      // check all kmer for each read
      for (size_t j = 0; j <= sequence.length()-kmer_length; ++j) {

        const std::string& kmer = sequence.substr(j, kmer_length);
        if (contains_non_gatc(kmer, kmer_length)) 
          continue;
        
        kmer_int_type_t kmer_val = kmer_to_intval(kmer, kmer_length);
        if (DS_MODE) 
          kmer_val = get_DS_kmer_val(kmer_val, kmer_length);
        
	// update it if exist
        if (exists(kmer_val)) {
          kmer_map[kmer_val].add(i,j);          
          continue;
        }

	// add a new kmer into kmer_map	
	T info;
        if (flag)
          info.add(i,j);
        else
          info.add(i,j,5);
        kmer_map[kmer_val] = info;
      }
    }

    time_t end = time(NULL);
    if (flag)
      std::cerr << "Kmer hash finished, get " << kmer_map.size() 
        << " kmers! (elapsed time: " << (end-beg) << " s)" << std::endl;
  }

/*
  void get_hash_from_kmers(std::string filename, unsigned int kmer_length) {
    std::fstream infile;
    infile.open(filename.c_str(), std::fstream::in);
    if (!infile.is_open()) {
      std::cerr << "error,the file " << filename << " can't be opened." << std::endl;
      exit(1);
    }
    std::cerr << "Begin loading kmers ..." << std::endl;
    string header, kmer;
    while(!infile.eof()) {
      getline(infile, header);
      if (header[0] == '>'){	
	getline(infile, kmer);
	if (kmer.length() == kmer_length) {
	  kmer_int_type_t intval = kmer_to_intval(kmer);
	  header.erase(0,1);
          size_t count = atoi(header.c_str());
	  kmer_map[intval] = count;
	}
      }
    }
    cerr << kmer_map.size() << " kmers loaded." << std::endl;
  }
*/

  // remove it by set its count to be 0
  bool remove (kmer_int_type_t kmer_val) {

    kmer_map_base_iterator_t it = find_kmer(kmer_val);

    if (it != kmer_map.end()) {
       (it->second).clear();
      //(it->second).count = 0;
      return(true);
    } else {
      return(false);
    }
  }

  // delete a kmer from hash table
  bool delete_kmer (kmer_int_type_t kmer_val) {

    kmer_map_base_iterator_t it = find_kmer(kmer_val);

    if (it != kmer_map.end()) {
      kmer_map.erase(it);
      return(true);
    } else {
      return(false);
    }
  }


  // exclude the singleton and low complexity kmers hash table 
  bool remove_erroneous_kmers(float min_ratio_non_error) {

    time_t beg = time(NULL);

    kmer_map_const_iterator_t it;
    std::vector<kmer_int_type_t> deletion_list;
    for (it = kmer_map.begin(); it != kmer_map.end(); ++it) {

      kmer_int_type_t kmer_val = it->first;
      std::vector<kmer_occurence_pair_t> candidates = get_forward_candidates(kmer_val);

      int dominant_count = 0;
      for (unsigned int i = 0; i < candidates.size(); ++i) {

        if (candidates[i].second) {
          int candidate_count = candidates[i].second;
          if (dominant_count == 0) {
            dominant_count = candidate_count;
          } else if ( (float) candidate_count/dominant_count < min_ratio_non_error ) {
            kmer_map_base_iterator_t candidate = find_kmer(candidates[i].first);
            deletion_list.push_back(candidate->first);
            candidate->second.clear(); // disable when encountered next time
          }
        }
      }

    }

    cerr << "remove erroneous kmers..." << endl;

    // remove the erroneous kmers
    if (!deletion_list.empty()) {

      for (unsigned int i = 0; i < deletion_list.size(); ++i) 
        delete_kmer(deletion_list[i]);
      
      time_t end = time(NULL);
      std::cerr << "Done! (elapsed time: " << (end-beg) << " s)" << std::endl;
      std::cerr << "kmer count after errors deletion: " << kmer_map.size() << std::endl;

      return (true);
    } else {

      return (false);
    }
  }

  bool prune_hash(int min_kmer_coverage, float min_kmer_entropy) {

    if (min_kmer_coverage == 1 && min_kmer_entropy < 1e-6) {
      return (true);
    }

    time_t beg = time(NULL);
    std::cerr << "Pruning kmer hash ..." << std::endl;

    kmer_map_base_iterator_t pos;
    for (pos = kmer_map.begin(); pos != kmer_map.end(); ) {
      std::string kmer = intval_to_kmer(pos->first, kmer_length);
      
      if ( static_cast<int>((pos->second).get_count()) < min_kmer_coverage
          || compute_entropy(pos->first, kmer_length) < min_kmer_entropy ) {
          //  Why error if using hash_map ??
          kmer_map.erase(pos++);
      } else {
        ++pos;
      }
    }

    time_t end = time(NULL);
    std::cerr << "Pruning kmer hash has been Done ! (elapsed time: " 
        << (end-beg) << " s)" << std::endl;

    std::cerr << "kmer count after prune: " << kmer_map.size() << std::endl;

    return (true);
  }


  // get sorted seeds
  std::vector<kmer_int_type_t> get_seed_sort_descending_counts() {

    std::vector<kmer_int_type_t> kmer_list;
    kmer_map_base_iterator_t it;
    for (it = kmer_map.begin(); it != kmer_map.end(); ++it) {
      if ( static_cast<int>((it->second).get_count()) < g_min_seed_coverage
        || compute_entropy(it->first, kmer_length) < g_min_seed_entropy ) {
          //  Why error if using hash_map ??
          continue;
      }
      kmer_list.push_back(it->first);
    }

    std::cerr << "Sorting kmers ..." << std::endl;

    kmer_sorter_by_count_desc_t sorter (*this);
    sort(kmer_list.begin(), kmer_list.end(), sorter);

    std::cerr << "Done sorting." << std::endl;
    return (kmer_list);
  }

  void get_seed_sort_descending_counts(std::vector<kmer_int_type_t> & kmer_list) {

    kmer_map_base_iterator_t it;
    for (it = kmer_map.begin(); it != kmer_map.end(); ++it) {
      if ( static_cast<int>((it->second).get_count()) < g_min_seed_coverage
        || compute_entropy(it->first, kmer_length) < g_min_seed_entropy ) {
          continue;
      }
      kmer_list.push_back(it->first);
    }

    time_t beg = time(NULL);
    std::cerr << "Sorting kmers ..." << std::endl;

    kmer_sorter_by_count_desc_t sorter(*this);
    sort(kmer_list.begin(), kmer_list.end(), sorter);

    time_t end = time(NULL);
    std::cerr << "Done sorting. (elapsed time: " << (end-beg) << " s)" << std::endl;
  }

  // get forward candidate kmers
  std::vector<kmer_occurence_pair_t> get_forward_candidates(kmer_int_type_t seed_kmer, bool getZero = false) {

    std::vector<kmer_occurence_pair_t> candidates;
    kmer_int_type_t forward_prefix 
      = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;

    for (kmer_int_type_t i = 0; i < 4; ++i) {

      kmer_occurence_pair_t candidate;
      candidate.first = forward_prefix | i;
      candidate.second = get_kmer_count(candidate.first);

      if (candidate.second) {
        candidates.push_back(candidate);      
      } else {
        if (getZero && reuse(candidate.first))
          candidates.push_back(candidate);
      }
    }

    kmer_sorter_by_count_desc_t sorter(*this);
    sort (candidates.begin(), candidates.end(), sorter);

    return (candidates);
  }

  void get_forward_candidates(kmer_int_type_t seed_kmer, std::vector<kmer_occurence_pair_t>& candidates, bool getZero = false) {

    candidates.clear(); // clear the vector
    kmer_int_type_t forward_prefix
      = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;

    for (kmer_int_type_t i = 0; i < 4; ++i) {

      kmer_occurence_pair_t candidate;
      candidate.first = forward_prefix | i;
      candidate.second = get_kmer_count(candidate.first);

      if (candidate.second) {
        candidates.push_back(candidate);
      } else {
        if (getZero && reuse(candidate.first))
          candidates.push_back(candidate);
      }

    }

    kmer_sorter_by_count_desc_t sorter(*this);
    sort (candidates.begin(), candidates.end(), sorter);

  }


  // get reverse candidate kmers
  std::vector<kmer_occurence_pair_t> get_reverse_candidates(kmer_int_type_t seed_kmer, bool getZero = false) {
    
    kmer_int_type_t reverse_suffix = seed_kmer >> 2;
    std::vector<kmer_occurence_pair_t> candidates;

    for (kmer_int_type_t i = 0; i < 4; ++i) {

      kmer_occurence_pair_t candidate;
      candidate.first = (i << (kmer_length*2-2)) | reverse_suffix;
      candidate.second = get_kmer_count(candidate.first);

      if (candidate.second) {
        candidates.push_back(candidate);
      } else {
        if (getZero && reuse(candidate.first))
          candidates.push_back(candidate);
      }

    }

    kmer_sorter_by_count_desc_t sorter(*this);
    sort (candidates.begin(), candidates.end(), sorter);

    return (candidates);
  }

  void get_reverse_candidates (kmer_int_type_t seed_kmer, std::vector<kmer_occurence_pair_t>& candidates, bool getZero = false) {

    candidates.clear(); // clear the vector    
    kmer_int_type_t reverse_suffix = seed_kmer >> 2;

    for (kmer_int_type_t i = 0; i < 4; ++i) {

      kmer_occurence_pair_t candidate;
      candidate.first = (i << (kmer_length*2-2)) | reverse_suffix;
      candidate.second = get_kmer_count(candidate.first);

      if (candidate.second) {
        candidates.push_back(candidate);
      } else {
        if (getZero && reuse(candidate.first))
          candidates.push_back(candidate);
      }

    }

    kmer_sorter_by_count_desc_t sorter(*this);
    sort (candidates.begin(), candidates.end(), sorter);
  }

}; // end of class KmerMap

#endif
