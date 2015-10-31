//// This file was modified from file KmerCounter.cpp
// in Inchworm modules of Trintiy.
// The original copyright info is listed below
// Copyright (c) 2010, The Broad Institute, Inc. 
// Distributed under the  Distributed under Trinity Software LICENSE
// (See accompanying file LICENSE)


#ifndef SPLICING_GRAPH_H
#define SPLICING_GRAPH_H

/*
 * splicing_graph.h
 */

#include "utility.h"
#include "kmerhash.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <iomanip>
#include <numeric>
#include <algorithm>

// serialize/restore data as a text stream
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// headers privating serialize() function for many c++ classes
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>


typedef std::pair<int,int> pair_t;


// functions
bool is_similar(const std::string& str1, const std::string& str2, char mode);
bool is_aligned(const std::string& str1, const std::string& str2);
bool compatible(const std::string& str1, const std::string& str2);


template<class T>
class SplicingGraph {

private:  //SplicingGraph

  typedef int node_idx_t;  // it can not be unsigned

  // define Node as a private class
  class Node {

    friend class boost::serialization::access;
    /*
     * When the class Archive corresponds to an output
     * Archive, the & operator is defined similar to <<.
     * Likewise, when the class Archive is a type of input
     * Archive, the & operator is defined similar to >>.
     */
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & sequence;
      ar & parents;
      ar & children;
    }

  public:
    Node(): sequence("") {};
    Node(const std::string& mysequence) {
      sequence = mysequence;
    }

    Node(const Node& node) {
      sequence = node.sequence;
      children = node.children;
      parents = node.parents;
    }

    void set_sequence(const std::string& myseq) {
      sequence = myseq;
    }

    std::string get_sequence() {
      return sequence;
    }

    bool add_child(node_idx_t child) {
      if (child < 0) 
        return false;
      if (!children.empty()) {
        for (size_t i = 0; i < children.size(); ++i) {
          if (children[i] == child ) // if exist already
            return false;
        }
      }
      this->children.push_back(child);
      return true;
    }

    bool add_parent(node_idx_t parent) {
      if (parent < 0) 
        return false; 
      if (! parents.empty()) {
        for (size_t i = 0; i < parents.size(); ++i) {
          if (parents[i] == parent )
            return false;
        }
      }
      this->parents.push_back(parent);
      return true;
    }

    bool is_child(node_idx_t child) {
      if (child < 0)
	return false;
      std::vector<node_idx_t>::iterator it = children.begin();
      for ( ; it != children.end(); ++it ) {
        if (*it == child) 
	  return true;
      }
      return false;
    }

    bool is_parent(node_idx_t parent) {
      if (parent < 0) 
	return false;
      std::vector<node_idx_t>::iterator it = parents.begin();
      for ( ; it != parents.end(); ++it ) {
        if (*it == parent)
          return true;
      }
      return false;
    }
 
    bool delete_child(node_idx_t child) {
      if (child < 0)
        return false;
      std::vector<node_idx_t>::iterator it = children.begin();
      for ( ; it != children.end(); ++it) {
        if (*it == child) 
          break;
      }
      if (it != children.end()) {
	children.erase(it);
        return true;
      } else {
        return false;
      }
    }
  
    bool delete_parent(node_idx_t parent) {
      if (parent < 0)
        return false;
      std::vector<node_idx_t>::iterator it = parents.begin();
      for ( ; it != parents.end(); ++it) {
        if (*it == parent) 
	  break;
      }
      if (it != parents.end()) {
        parents.erase(it);
        return true;
      } else {
        return false;
      }
    }
 
    void clear_children() {
      children.clear();
    }

    void clear_parents() {
      parents.clear();
    }

    void clear() {
      sequence.clear();
      children.clear();
      parents.clear();
    }

  public:
    std::string sequence;
    std::vector<node_idx_t> parents;
    std::vector<node_idx_t> children;

  };

  class node_sorter_by_seq_length_t {

  public:
    node_sorter_by_seq_length_t(std::vector<Node>& ns) : node_s(ns) {};
    bool operator() (const node_idx_t i, const node_idx_t j) {
      return (node_s[i].sequence.length() < node_s[j].sequence.length() ||
	       (node_s[i].sequence.length() == node_s[j].sequence.length() && i < j));
    }
  private:
    std::vector<Node>& node_s;

  };

public:  //SplicingGraph

  SplicingGraph() {
    size_ = 0;
  }

  size_t get_size() {
    return size_;
  } 
  
  double compute_coverage(const std::string& seq) {

    int cov_size = static_cast<int>(seq.length())-static_cast<int>(g_kmer_length)+1;
    if (cov_size <= 0) {
      //std::cerr << "[Warning] coverage of seq " << seq << " is 0." << std::endl;
      return 0.0f;
    }

    std::vector<int> cov_v(cov_size);
    for (int i = 0; i < cov_size; ++i) {
      const std::string& kmer = seq.substr(i, g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);
      cov_v[i] = get_kmer_count(intval);
    }

    sort(cov_v.begin(), cov_v.end());

    // compute average coverage based on some "confidence interval", 
    // considering non-uniform frequency distribution of kmers.
    int quantile = static_cast<int>(0.05*cov_size + 0.5);
    std::vector<int>::iterator first = cov_v.begin();
    std::vector<int>::iterator last = cov_v.end();
    if (quantile > 0) {
      first = cov_v.begin() + quantile;
      last = cov_v.end() - quantile;
    }
    int sum = 0;
    for (; first != last; ++first) {
      sum += *first;
    }

    return static_cast<double>(sum*1.0/(cov_size-2*quantile));
  } 

  std::string show_coverage_detail(const std::string& sequence) {

    std::stringstream cov;
    if ((int)sequence.length() < g_kmer_length)
      return cov.str();
    
    for (int i = 0; i <= (int)sequence.length()-g_kmer_length; ++i) {
      const std::string& kmer = sequence.substr(i, g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);
      cov << get_kmer_count(intval) << ",";
    }
    cov << std::endl;

    return cov.str();
  }

  bool is_used(kmer_int_type_t intval) {
    if (g_double_stranded_mode)
       intval = get_DS_kmer_val(intval, g_kmer_length);
    return (used_kmers_.find(intval) != used_kmers_.end());
  }

  // return 0 if intval does not exist in used_kmers_
  size_t get_kmer_count(kmer_int_type_t intval) {

    if (g_double_stranded_mode)
       intval = get_DS_kmer_val(intval, g_kmer_length);
 
    if (is_used(intval))
      return used_kmers_[intval];
    else
      return 0;
  }

  size_t get_kmer_count(const std::string& kmer) {

    kmer_int_type_t intval = kmer_to_intval(kmer);
    return get_kmer_count(intval);
  }


  bool add_used_kmer(kmer_int_type_t intval, size_t cov) {

    if (g_double_stranded_mode)
       intval = get_DS_kmer_val(intval, g_kmer_length);
    used_kmers_[intval] = cov;

    return true;
  }

  bool add_used_kmers(KmerMap<T>& kmer_map, const std::string str) {

    int str_len = str.length();
    if (str_len >= g_kmer_length) {
      for (int i = 0; i <= str_len - g_kmer_length; ++i) {
        kmer_int_type_t intval = kmer_to_intval(str.substr(i, g_kmer_length));
        size_t cov = kmer_map.get_kmer_count(intval);
        add_used_kmer(intval, cov);
      }
    }

    return true;
  }


  bool delete_kmer(kmer_int_type_t intval) {

    if (g_double_stranded_mode)
       intval = get_DS_kmer_val(intval, g_kmer_length);

    if (used_kmers_.find(intval) != used_kmers_.end()) 
      used_kmers_.erase(used_kmers_.find(intval));

    return true;
  }

  std::string describe_graph () {
   
    std::stringstream graph;

    graph << "** Edges **" << std::endl;

    for (size_t i = 0; i < size_; ++i) {
      if (node_set_[i].children.empty())
        continue;
      for (size_t j = 0; j < node_set_[i].children.size(); ++j) {
        graph << i << "->" << node_set_[i].children[j];
	if (j < node_set_[i].children.size()-1) 
	  graph << ", ";
        else 
	  graph << ";" << std::endl;
      }
    }

    graph << "** Nodes **" << std::endl;

    for (size_t i = 0; i < size_; ++i) {

      graph << "node id = " << i << "\tlength = "
        << node_set_[i].sequence.length();

      std::string sequence = node_set_[i].sequence;
      // extend it if sequence is too short
      if ((int)sequence.length() < g_kmer_length) {
	// only one parent
	if (node_set_[i].parents.size() == 1) {
	  node_idx_t j = node_set_[i].parents[0];
	  const std::string& parent = node_set_[j].sequence;
	  if ((int)parent.length() < g_kmer_length) {
	    sequence = parent + sequence;
	  } else {
	    // only last (k-1) bases
	    sequence = parent.substr(parent.length()-g_kmer_length+1) + sequence;
	  }
	}
	// only one child
        if (node_set_[i].children.size() == 1) {
	  node_idx_t j = node_set_[i].children[0];
	  const std::string& child = node_set_[j].sequence;
	  sequence += child.substr(0, g_kmer_length-1);
	}
      }

      graph << "\tcov = " << setiosflags(ios::fixed) 
	<< setprecision(2) << compute_coverage(sequence)
        << "\tsequence : " << std::endl
        <<  node_set_[i].sequence << std::endl;
      graph << show_coverage_detail(sequence);      
    }

    return graph.str();
  }

  std::string forward_extend(KmerMap<T>& kmer_map, kmer_int_type_t kmer_val, 
   std::vector<kmer_int_type_t>& bifurcation_points) {

    //cerr << "forward_extend... " << endl;
    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    size_t sum = 0;
    std::vector<kmer_occurence_pair_t> candidates;

    while (1) {
      //candidates = kmer_map.get_forward_candidates(intval);
      kmer_map.get_forward_candidates(intval, candidates);
      if (candidates.empty()) break; // no candidates
      kmer_int_type_t candidate;
      bool flag = false;  // indicate if there exist a unused candidate
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (!is_used(candidates[i].first)) {
          flag = true;
          candidate = candidates[i].first;
          // record the kmer as a bifurcation point if it has alternative candidate
          if (i < candidates.size()-1) 
            bifurcation_points.push_back(intval);
          break;
        }
      }
      if (!flag) break;  // all candidates have been used before
      size_t cov = kmer_map[candidate].get_count(); 
      add_used_kmer(candidate, cov);
      sum = sum + cov;
      int base_num = candidate & 3ll;
      char base = int_to_base(base_num);
      str += base;
      // update the kmer to be extended
      intval = candidate;
    }

    // check average coverage of the contig
    if ((int)str.length() > g_kmer_length) {
      int avg_cov = static_cast<int>(sum*1.0/(str.length()-g_kmer_length) + 0.2);
      if (avg_cov < g_min_average_coverage) {
	//str = intval_to_kmer(kmer_val, g_kmer_length);
	str.erase(str.begin()+g_kmer_length, str.end());
      }
    }

    return str;
  }  

  std::string forward_extend(KmerMap<T>& kmer_map, kmer_int_type_t kmer_val, bool getZero = false) {

    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    size_t sum = 0;

    std::vector<kmer_occurence_pair_t> candidates;
    while (1) {

      //candidates = kmer_map.get_forward_candidates(intval);
      kmer_map.get_forward_candidates(intval, candidates, getZero);
      if (candidates.empty()) break; // no candidates
      kmer_int_type_t candidate;
      bool flag = false;  // indicate if there exist a unused candidate
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (!is_used(candidates[i].first)) {
          flag = true;
          candidate = candidates[i].first;
          break;
        }
      }

      if (!flag) break;  // all candidates have been used before
      size_t cov = kmer_map[candidate].get_count();
      add_used_kmer(candidate, cov);
      sum = sum + cov;
      int base_num = candidate & 3ll;
      char base = int_to_base(base_num);
      str += base;

      // update the kmer to be extended
      intval = candidate;
    }
 
    // check average coverage of the contig
    if ((int)str.length() > g_kmer_length) {
      int avg_cov = static_cast<int>(sum*1.0/(str.length()-g_kmer_length) + 0.2);
      if (avg_cov < g_min_average_coverage) {
        //str = intval_to_kmer(kmer_val, g_kmer_length);
        str.erase(str.begin()+g_kmer_length, str.end());
      }
    }

    return str;
  }


  std::string reverse_extend(KmerMap<T>& kmer_map, kmer_int_type_t kmer_val, bool getZero = false) {

    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    size_t sum = 0;

    std::vector<kmer_occurence_pair_t> candidates;
    while (true) {
      
      kmer_map.get_reverse_candidates(intval, candidates, getZero);
      if (candidates.empty()) break;  
      kmer_int_type_t candidate;
      bool flag = false;  // if exist a candidate that has not been used
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (!is_used(candidates[i].first)) {
          flag = true;
          candidate = candidates[i].first;
          break;
        }
      }

      if (!flag) break;  // all candidates have been used before
      size_t cov = kmer_map[candidate].get_count();
      add_used_kmer(candidate, cov);
      sum = sum + cov;
      // get the extended base
      int base_num = (candidate >> (g_kmer_length*2-2)) & 3ll;
      char base = int_to_base(base_num);
      str = base + str;

      intval = candidate;
    }


    // check average coverage 
    if ((int)str.length() > g_kmer_length) {
      int avg_cov = static_cast<int>(sum*1.0/(str.length()-g_kmer_length) + 0.2);
      if (avg_cov < g_min_average_coverage) {
        str.erase(str.begin()+g_kmer_length, str.end());
      }
    }

    return str;
  }

  std::string forward_extend_for_completeness(KmerMap<T>& kmer_map, kmer_int_type_t kmer_val) {

    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    std::map<kmer_int_type_t, bool> kmers;
    std::vector<kmer_occurence_pair_t> candidates;
    int repeat = 0;
    int start = g_kmer_length;

    while(1) {
      kmer_map.get_forward_candidates(intval, candidates);
      if (candidates.empty()) break; 
      kmer_int_type_t candidate;
      bool flag = false; 
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (!is_used(candidates[i].first)) {
          flag = true;
          candidate = candidates[i].first;
          break;
        } else if (kmers.find(candidates[i].first) == kmers.end()) {
          flag = true;
          candidate = candidates[i].first;
          if (repeat == 0) 
            start = str.length();
          repeat++;
          break;
        }
      }

      if (!flag) break;  // all candidates have been used
      size_t cov = kmer_map[candidate].get_count();
      add_used_kmer(candidate, cov);
      kmers[candidate] = true;
      int base_num = candidate & 3ll;
      char base = int_to_base(base_num);
      str += base;

      // update the kmer to be extended
      intval = candidate;
    }

    if (repeat > 0.5*str.length())
      return str.substr(0, start);

    return str;
  }

  std::string reverse_extend_for_completeness(KmerMap<T>& kmer_map, kmer_int_type_t kmer_val) {

    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    std::vector<kmer_occurence_pair_t> candidates;
    std::map<kmer_int_type_t, bool> kmers;

    int repeat = 0;
    int start = g_kmer_length;
    while (true) {

      kmer_map.get_reverse_candidates(intval, candidates);
      if (candidates.empty()) break;
      kmer_int_type_t candidate;
      bool flag = false;  // if exist a candidate that has not been used
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (!is_used(candidates[i].first)) {
          flag = true;
          candidate = candidates[i].first;
          break;
        } else if (kmers.find(candidates[i].first) == kmers.end()) {
          flag = true;
          candidate = candidates[i].first;
          if (repeat == 0)
            start = str.length();
          repeat++;
          break;
        }
      }

      if (!flag) break;  // all candidates have been used before
      size_t cov = kmer_map[candidate].get_count();
      add_used_kmer(candidate, cov);
      kmers[candidate] = true;
      int base_num = (candidate >> (g_kmer_length*2-2)) & 3ll;
      char base = int_to_base(base_num);
      str = base + str;

      intval = candidate;
    }

    if (repeat > 0.5*str.length())
      return str.substr(str.length()-start);

    return str;
  }



  kmer_int_type_t get_reverse_end_kmer(KmerMap<T>& kmer_map, kmer_int_type_t seed_val,
   kmer_int_type_t stop_kmer = 0ULL, bool getZero = true) {

    kmer_int_type_t intval = seed_val;
    std::vector<kmer_occurence_pair_t> candidates;
    std::map<kmer_int_type_t, bool> kmers;

    while (1) {

      kmer_map.get_reverse_candidates(intval, candidates, getZero);
      if (candidates.empty()) 
        break;
      kmer_int_type_t candidate;
      bool flag = false;  // if there exists a candidate or not
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (kmers.find(candidates[i].first) == kmers.end()) {
          flag = true;
          candidate = candidates[i].first;
          break;
        }
      }
      
      if ((!flag) || candidate == stop_kmer)
        break; 
      kmers[candidate] = true;
      intval = candidate;
    }

    return intval;
  }

  kmer_int_type_t get_forward_end_kmer(KmerMap<T>& kmer_map, kmer_int_type_t seed_val,
   kmer_int_type_t stop_kmer = 0ULL, bool getZero = true) {

    kmer_int_type_t intval = seed_val;
    std::map<kmer_int_type_t, bool> kmers; 
    std::vector<kmer_occurence_pair_t> candidates;

    while (1) {

      kmer_map.get_forward_candidates(intval, candidates, getZero);
      if (candidates.empty()) 
        break; // no candidates
      kmer_int_type_t candidate;
      bool flag = false;  // indicate if there exist a unused candidate
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (kmers.find(candidates[i].first) == kmers.end()) {
          flag = true;
          candidate = candidates[i].first;
          break;
        }
      }

      if (!flag || candidate == stop_kmer) 
        break;
      kmers[candidate] = true;
      intval = candidate;
    }

    return intval;
  }

  kmer_int_type_t get_reverse_end_kmer(KmerMap<T>& kmer_map, kmer_int_type_t seed_val, std::string& str,
   kmer_int_type_t stop_kmer = 0ULL, bool getZero = true) {

    kmer_int_type_t intval = seed_val;
    str = intval_to_kmer(intval, g_kmer_length);
    std::vector<kmer_occurence_pair_t> candidates;
    std::map<kmer_int_type_t, bool> kmers;

    while (1) {

      kmer_map.get_reverse_candidates(intval, candidates, getZero);
      if (candidates.empty())
        break;

      kmer_int_type_t candidate;
      bool flag = false;  // if there exists a candidate or not
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (kmers.find(candidates[i].first) == kmers.end()) {
          flag = true;
          candidate = candidates[i].first;
          break;
        }
      }

      if ((!flag) || candidate == stop_kmer)
        break;

      kmers[candidate] = true;
      int base_num = (candidate >> (g_kmer_length*2-2)) & 3ll;
      char base = int_to_base(base_num);
      str = base + str;

      intval = candidate;
    }

    return intval;
  }

  kmer_int_type_t get_forward_end_kmer(KmerMap<T>& kmer_map, kmer_int_type_t seed_val, std::string & str,
   kmer_int_type_t stop_kmer = 0ULL, bool getZero = true) {

    kmer_int_type_t intval = seed_val;
    str = intval_to_kmer(intval, g_kmer_length);
    std::map<kmer_int_type_t, bool> kmers;
    std::vector<kmer_occurence_pair_t> candidates;

    while (1) {

      kmer_map.get_forward_candidates(intval, candidates, getZero);
      if (candidates.empty())
        break; // no candidates

      kmer_int_type_t candidate;
      bool flag = false;  // indicate if there exist a unused candidate
      for (size_t i = 0; i < candidates.size(); ++i) {
        if (kmers.find(candidates[i].first) == kmers.end()) {
          flag = true;
          candidate = candidates[i].first;
          break;
        }
      }

      if (!flag || candidate == stop_kmer)
        break;

      kmers[candidate] = true;
      int base_num = candidate & 3ll;
      char base = int_to_base(base_num);
      str += base;

      intval = candidate;
    }

    return intval;
  }


  bool get_trunk(KmerMap<T>& kmer_map, std::vector<std::string>& data, kmer_int_type_t seed) {

    add_used_kmer(seed, kmer_map.get_kmer_count(seed));

    //extend seed kmer in two direction
    const std::string& left = reverse_extend(kmer_map, seed);
    const std::string& right = forward_extend(kmer_map, seed);
    std::string trunk = left + right.substr(g_kmer_length);

    if (static_cast<int>(trunk.length()) < g_min_trunk_length)
      return false;

    //cerr << "trunk :\n" << trunk  << endl << show_coverage_detail(trunk) << endl;
    Node node(trunk);
    int p = add_node(node);  //add trunk

    refine_trunk(kmer_map);

    if (g_is_paired_end) {

      // extend by pair info
      forward_extend_by_pair_info(kmer_map, data, p);
      //cerr << "After forward extend :\n" << endl << node_set_[0].sequence << endl;
      reverse_extend_by_pair_info(kmer_map, data, p);
      //cerr << "After reverse extend :\n" << endl << node_set_[0].sequence << endl;
    } else {

       //allows used kmers
       const std::string& start_kmer = node_set_[p].sequence.substr(0, g_kmer_length);
       kmer_int_type_t reverse_seed = kmer_to_intval(start_kmer);
       const std::string& str1 = reverse_extend(kmer_map, reverse_seed, true);

       const std::string& end_kmer = node_set_[p].sequence.substr(node_set_[p].sequence.length() - g_kmer_length, g_kmer_length);
       kmer_int_type_t forward_seed = kmer_to_intval(end_kmer);
       const std::string& str2 = forward_extend(kmer_map, forward_seed, true);

       node_set_[p].sequence = str1 + node_set_[p].sequence.substr(g_kmer_length) + str2.substr(g_kmer_length);
    }

    return true;
  }
 
  void refine_trunk(KmerMap<T>& kmer_map) {

    // reverse
    std::vector<kmer_occurence_pair_t> candidates;
    for (int i = 0; i < 2*g_kmer_length; ++i) {
      const std::string& kmer = node_set_[0].sequence.substr(i,g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);
      kmer_map.get_reverse_candidates(intval, candidates);
      if (candidates.size() <= 1)
        continue;
      bool hit_point = false;
      for (size_t j = 0; j < candidates.size(); ++j) {
        if (!is_used(candidates[j].first)) {
          hit_point = true;
          break;
        }
      }
      if (hit_point) {
        const std::string & sequence = reverse_extend(kmer_map, intval);
        if ((int)sequence.length() > g_kmer_length + g_min_exon_length) {
	  node_set_[0].sequence = sequence + node_set_[0].sequence.substr(i+g_kmer_length);
          break;
        }
      }
    }

   // forward
   int length = node_set_[0].sequence.length();
   assert(length >= 3*g_kmer_length);
   for (int i = length - 3*g_kmer_length; i < length-g_kmer_length; ++i) {
      const std::string& kmer = node_set_[0].sequence.substr(i,g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);
      kmer_map.get_forward_candidates(intval, candidates);
      if (candidates.size() <= 1)
        continue;
      bool hit_point = false;
      for (size_t j = 0; j < candidates.size(); ++j) {
        if (!is_used(candidates[j].first)) {
          hit_point = true;
          break;
        }
      }
      if (hit_point) {
        const std::string & sequence = forward_extend(kmer_map, intval);
        if ((int)sequence.length() > g_kmer_length + g_min_exon_length) {
          node_set_[0].sequence = node_set_[0].sequence.substr(0,i) + sequence;
          break;
        }
      }
    }


  }

  node_idx_t add_node(Node& node) {

    node_set_.push_back(node);
    return (size_++);
  }

  // return the index of node containing the anchor 
  node_idx_t find_node_index (const std::string& anchor) {

    int idx = -1;

    for (size_t i = 0; i < size_ ; ++i) {
      if (node_set_[i].sequence.find(anchor) != std::string::npos) {
        idx = i;
        break;
      }
    }

    return idx;
  }

  /*
     Is a circle generated if add q to be child of p ?
     Acutally, we only need to check if p is the offspring of q.
     return true if there exists path q -> ... -> p
  */
  bool is_circle(const node_idx_t p, const node_idx_t q, std::set<node_idx_t>& checked) {

   // checked should be non-const reference because it need iterate many time
    bool flag = false;

    if (node_set_[q].children.empty() || 
	checked.find(q) != checked.end()) {
      return false;
    } else {
      checked.insert(q);
    }

    std::vector<node_idx_t>::iterator it;
    for (it = node_set_[q].children.begin(); 
     it != node_set_[q].children.end(); ++it) {
      if ( *it == p ) {
        flag = true;
        break;
      } else {
        flag = is_circle(p, *it, checked);
        if (flag) 
          break;
      }
    }

    return flag;
  }

  // return true if there exists path p -> ... -> q 
  bool has_path(const node_idx_t p, const node_idx_t q) {

    std::set<node_idx_t> checked;
    return is_circle(q, p, checked);
  }

  bool erase_node(const node_idx_t p) {

    if (p >= size_)
      return (false); 
    node_set_.erase(node_set_.begin()+p);
    size_--;

    // update each node
    for (size_t i = 0; i < size_ ; ++i) {

      // update its children if any
      for (size_t j = 0; j < node_set_[i].children.size(); ) {
        if (node_set_[i].children[j] < p) {
	  ++j;
        } else if (node_set_[i].children[j] == p) {
          node_set_[i].children.erase(node_set_[i].children.begin()+j);
	  if (j == node_set_[i].children.size()) 
            break;
        } else if (node_set_[i].children[j] > p) {
          node_set_[i].children[j]--;
	  ++j;
        }
      }

      // update its parents if any
      for (size_t j = 0; j < node_set_[i].parents.size();) {
	if (node_set_[i].parents[j] < p) {
 	  ++j;
        } else if (node_set_[i].parents[j] == p) {
          node_set_[i].parents.erase(node_set_[i].parents.begin()+j);
	  if (j == node_set_[i].parents.size()) 
            break;
        } else if (node_set_[i].parents[j] > p) {
          node_set_[i].parents[j]--;
	  ++j;
        }
      }
    }

    return (true);
  }

  // restore some used kmers to kmer_hash
  void restore_kmers(const std::string& sequence) {

    if ((int)sequence.length() < g_kmer_length)
      return;

    for (int i = 0; i <= (int)sequence.length() - g_kmer_length; ++i) {

      const std::string& kmer = sequence.substr(i, g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);

      if (g_double_stranded_mode)
	intval = get_DS_kmer_val(intval, g_kmer_length);

      restore_kmers_.insert(intval);
    }

  }

  void grow_and_branch(KmerMap<T>& kmer_map, std::vector<std::string>& data, node_idx_t p,
    std::vector<kmer_int_type_t>& bifurcation_points) {

      //std::cerr << "grow_and_branch" << std::endl;

      while (bifurcation_points.size() > 0) {

        kmer_int_type_t intval = bifurcation_points.back();
        bifurcation_points.pop_back();

        std::vector<kmer_int_type_t> bifurcation_points1;
        std::string sequence = forward_extend(kmer_map, intval, bifurcation_points1);  // can not be const
	if (bifurcation_points1.size() > 0 && bifurcation_points1[0] == intval ) {
          bifurcation_points.push_back(intval);
          bifurcation_points1.erase(bifurcation_points1.begin());
        }

	// the end kmer
        const std::string& endkmer = sequence.substr(sequence.length()-g_kmer_length);
	kmer_int_type_t end_val = kmer_to_intval (endkmer);
        std::vector<kmer_occurence_pair_t> candidates;
        kmer_map.get_forward_candidates(end_val, candidates);

        if (candidates.size() > 0) {  // add bubble

	  node_idx_t r = add_bubble(p, sequence, kmer_map);
          if (r > 0) {   
            forward_branches.insert(r);
          } else if (r == -2) { // return -2 if branch
	    r = add_branch(p, sequence, kmer_map, data);
	    if (r > 0) {
              forward_branches.insert(r);
	      node_set_[r].sequence = node_set_[r].sequence.substr(g_kmer_length); 
            }
          }
	} else { // add a branch

          node_idx_t r = add_branch(p, sequence, kmer_map, data);
          if (r > 0) {  
            forward_branches.insert(r);
            // change sequence of node r
            node_set_[r].sequence = node_set_[r].sequence.substr(g_kmer_length);
          }
        }  
 
      } // while
  }

  // make sure this function is only used for paired data
  bool check_forward_branch_with_pair_info(KmerMap<T>& kmer_map, std::vector<std::string>& data, const std::string& branch, size_t count) {

    bool is_branch = false;
    const string& check_seq = branch.substr(g_kmer_length, 2* g_kmer_length);

    size_t max_read_id = data.size() / 2;
    std::set<size_t> reads;
    set_reads(kmer_map, check_seq, reads);

    std::set<kmer_int_type_t> kmers_in_branch;
    if (g_double_stranded_mode) {
      for (int i = 0; i <= (int)branch.length()-g_kmer_length; ++i) {
        kmer_int_type_t intval = kmer_to_intval(branch.substr(i,g_kmer_length));
        kmer_int_type_t revcomp_intval = revcomp_val(intval, g_kmer_length);
        kmers_in_branch.insert(intval);
        kmers_in_branch.insert(revcomp_intval);
      }
    }

    int support = 0;
    int unsupport = 0;
    std::set<size_t>::const_iterator its;
    for (its = reads.begin(); its != reads.end(); ++its) {

      if (g_double_stranded_mode) {
        if (compatible(branch,data[*its])) {
          size_t mate_id ;
          if (*its < max_read_id)
            mate_id = *its + max_read_id;
          else
            mate_id = *its - max_read_id;
          const std::string& mate_read = data[mate_id];
          if (mate_read.find("N") != std::string::npos)
            continue;
         
          unsigned int used = 0;
          for (int i = 0; i <= (int)mate_read.length()-g_kmer_length; ++i) {
            kmer_int_type_t intval = kmer_to_intval(mate_read.substr(i,g_kmer_length));
            // used, but not used in this branch
            if (is_used(intval) && kmers_in_branch.find(intval) == kmers_in_branch.end())
              used++;
          }
          if (used == 0) {
            unsupport++;
            if ((unsupport >= g_min_ratio_welds * count) && (unsupport >= 0.5 * reads.size() + 2 ))
              break;
          } else if (used + 5 >= mate_read.length() - g_kmer_length) {
            support++;
            if ((support >= g_min_ratio_welds * count) && (support >= g_min_reads_span_junction)) {
              is_branch = true;
              break;
            }
          }
        }
      } else { // strand-specfic mode

        if (g_fr_strand == 1) {   // fr-firststrand

          if (*its < max_read_id) { 
            const std::string& mate_left = data[*its];
            if (compatible(branch, mate_left)) {
              const std::string& mate_right = data[*its+max_read_id];
              if (mate_right.find("N") != std::string::npos)
                continue;
              
              unsigned int used = 0;
              for (unsigned int i = 0; i <= mate_right.length()-g_kmer_length; ++i) {
                kmer_int_type_t intval = kmer_to_intval(mate_right.substr(i,g_kmer_length));
                if (is_used(intval))
                  used++;
              }
              if (used == 0) {
                unsupport++;
                if ((unsupport >= g_min_ratio_welds * count) && (unsupport >= g_min_reads_span_junction))
                  break;
              } else if (used + 5 >= mate_right.length() - g_kmer_length) {
                support++;
                if ((support >= g_min_ratio_welds * count) && (support >= g_min_reads_span_junction)) {
                  is_branch = true;
                  break;
                }
              }
            }
          }
        } else { // fr-secondstrand

           if (*its >= max_read_id) {  
            const std::string& mate_right = data[*its];
            if (compatible(branch, mate_right)) {
              const std::string& mate_left = data[*its-max_read_id];
              //cerr<< "left " << mate_left << endl;
              //cerr << "right" << mate_right << endl;
              if (mate_left.find("N") != std::string::npos)
                continue;
              int used = 0;
              for (int i = 0; i <= (int)mate_left.length() - g_kmer_length; ++i) {
                kmer_int_type_t intval = kmer_to_intval(mate_left.substr(i,g_kmer_length));
                if (is_used(intval))
                  used++;
              }
              if (used == 0) {
                unsupport++;
                if ((unsupport >= g_min_ratio_welds * count) && (unsupport >= g_min_reads_span_junction))
                  break;
              } else if (used + 5 >= (int)mate_left.length() - g_kmer_length) {
                support++;
                if ((support >= g_min_ratio_welds * count) && (support >= g_min_reads_span_junction)) {
                  is_branch = true;
                  break;
                }
              }
            }
          }
        }
      } // else
    }

    return is_branch;
  }

  // make sure this function is only used for paired data
  bool check_reverse_branch_with_pair_info(KmerMap<T>& kmer_map, std::vector<std::string>& data, const std::string& str, size_t count) {

    bool is_branch = false;
    size_t max_read_id = data.size() / 2;
    int str_length = static_cast<int>(str.length())-g_kmer_length;
    std::set<size_t> reads;
    int pos = str_length > 2*g_kmer_length ? (str_length - 2*g_kmer_length) : 0;
    set_reads(kmer_map, str.substr(pos, 2*g_kmer_length), reads);

    std::set<kmer_int_type_t> kmers_in_branch;
    if (g_double_stranded_mode) {
      for (int i = 0; i <= (int)str.length()-g_kmer_length; ++i) {
        kmer_int_type_t intval = kmer_to_intval(str.substr(i,g_kmer_length));
        kmer_int_type_t revcomp_intval = revcomp_val(intval, g_kmer_length);
        kmers_in_branch.insert(intval);
        kmers_in_branch.insert(revcomp_intval);
      }
    }

    int support = 0;
    int unsupport = 0;
    std::set<size_t>::const_iterator its;
    for (its = reads.begin(); its != reads.end(); ++its) {

      if (g_double_stranded_mode) { // non-strand specific
        if (compatible(str, data[*its])) {
          size_t mate_id ;
          if (*its < max_read_id)
            mate_id = *its + max_read_id;
          else
            mate_id = *its - max_read_id;
          const std::string& mate_read = data[mate_id];
          if (mate_read.find("N") != std::string::npos)
            continue;
     
          int used = 0;
          for (int i = 0; i <= (int)mate_read.length()-g_kmer_length; ++i) {
            kmer_int_type_t intval = kmer_to_intval(mate_read.substr(i,g_kmer_length));
            // used, but not used in this branch
            if (is_used(intval) && kmers_in_branch.find(intval) == kmers_in_branch.end())
              used++;
          }
          if (used == 0) {
            unsupport++;
            if ((unsupport >= g_min_ratio_welds * count) && (unsupport >= 0.5 * reads.size() + 2))
              break;
          } else if (used + 5 >= (int)mate_read.length() - g_kmer_length) {
            support++;
            if ((support >= g_min_ratio_welds * count) && (support >= g_min_reads_span_junction)) {
              is_branch = true;
              break;
            }
          }
        }
      } else { // strand-specfic mode

        if (g_fr_strand == 1) { // fr-firststrand
          if (*its >= max_read_id) { // --2--> ....... <--1--
            const std::string& mate_right = data[*its];
            if (compatible(str, mate_right)) {
              const std::string& mate_left = data[*its-max_read_id];
              if (mate_left.find("N") != std::string::npos)
                continue;
              int used = 0;
              for (int i = 0; i <= (int)mate_left.length() - g_kmer_length; ++i) {
                kmer_int_type_t intval = kmer_to_intval(mate_left.substr(i,g_kmer_length));
                if (is_used(intval))
                  used++;
              }

              if (used == 0) {
                unsupport++;
                if ((unsupport > g_min_ratio_welds * count) && (unsupport >= g_min_reads_span_junction))
                  break;
              } else if (used + 5 >= (int)mate_left.length()-g_kmer_length) {
                support++;
                if ((support >= g_min_ratio_welds * count) && (support >= g_min_reads_span_junction)) {
                  is_branch = true;
                  break;
                }
              }
            }
          }
        } else if (g_fr_strand == 2) { //--1--> ....... <--2--
          if (*its < max_read_id) {
            const std::string& mate_left = data[*its];
            if (compatible(str, mate_left)) {
              const std::string& mate_right = data[*its+max_read_id];
              if (mate_right.find("N") != std::string::npos)
                continue;
              unsigned int used = 0;
              for (unsigned int i = 0; i <= mate_right.length()-g_kmer_length; ++i) {
                kmer_int_type_t intval = kmer_to_intval(mate_right.substr(i,g_kmer_length));
                if (is_used(intval))
                  used++;
              }
              if (used == 0) {
                unsupport++;
                if ((unsupport >= g_min_ratio_welds * count) && (unsupport >= g_min_reads_span_junction))
                  break;
              } else if (used + 5 >= mate_right.length() - g_kmer_length) {
                support++;
                if ((support >= g_min_ratio_welds * count) && (support >= g_min_reads_span_junction)) {
                  is_branch = true;
                  break;
                }
              }
            }
          }
        }

      } // else (strand-spefic) 
    }

    return is_branch;
  }

  node_idx_t add_branch(const node_idx_t p, const std::string& branch, KmerMap<T>& kmer_map, std::vector<std::string>& data) {

    if ((int)node_set_[p].sequence.length() <= g_kmer_length)
      return -1;

    if ((int)branch.length() > g_kmer_length + g_min_exon_length) {

      if (g_debug) 
        std::cerr << "add_branch:" << branch << std::endl;

      const std::string& startkmer = branch.substr(0, g_kmer_length);
      std::string::size_type start = node_set_[p].sequence.find(startkmer);
      // the anchor of branch can be less than g_kmer_length
      if (start == std::string::npos) {
        restore_kmers(branch.substr(1));
        return -1; 
      }

      //if (start == 0) return -1;  // necessary if grow on a branch 
      if (p == 0 && (int)start < g_pair_gap_length) { // empirical trick
        restore_kmers(branch.substr(1));
        return -1;
      }

      if (g_is_paired_end) { 

        size_t count = get_kmer_count(startkmer);
        bool is_branch = check_forward_branch_with_pair_info(kmer_map, data, branch, count);
        
        if (!is_branch) {
          restore_kmers(branch.substr(1));
          return -1;
        }
      }


      Node node1, node2;
      
      if (start+g_kmer_length < node_set_[p].sequence.length())
        node1.sequence = node_set_[p].sequence.substr(start+g_kmer_length);
      
      node2.sequence = branch;
 
      if (!g_is_paired_end) {
        // check the junction
        const std::string& junc_branch_kmer = branch.substr(static_cast<int>(0.5*g_kmer_length),g_kmer_length);            
        if ((int)get_kmer_count(junc_branch_kmer) < g_min_junction_coverage) {
  	  if ((int)branch.length() > g_min_transcript_length) 
            restore_kmers(branch.substr(1));	  
          return -1;
        }
      }

      if (is_similar(branch.substr(g_kmer_length), node1.sequence, 'F')) {
        if (node_set_[p].children.empty() && (int)node1.sequence.length() < g_min_exon_length) {
	  node_set_[p].sequence = node_set_[p].sequence.substr(0, start) + branch;
          return -2;
        } else if ((int)node1.sequence.length() + g_kmer_length > (int)branch.length()) {
	  return -1;
        }
      } 

      if (size_ == 1 && (int)node1.sequence.length() < 2*g_kmer_length) {
        node_set_[p].sequence = node_set_[p].sequence.substr(0, start) + branch;
        return -2;
      }      

      if (!node1.sequence.empty()) {
        int q1 = add_node(node1);
        node_set_[q1].children = node_set_[p].children;
        node_set_[p].children.clear();
        node_set_[p].add_child(q1);
      }

      int q2 = add_node(node2);
      node_set_[p].add_child(q2);
      // update sequence of node p
      node_set_[p].sequence = node_set_[p].sequence.substr(0, start+g_kmer_length);

      //cerr << describe_graph() << endl;
      return q2;  // return new branch

    } else { 
      // ingore too short branch
      return -1;
    }

  }

  node_idx_t add_bubble(const node_idx_t p, const std::string& sequence, KmerMap<T>& kmer_map) {

    if ((int)sequence.length() < 2 * g_min_anchor_length)
      return -1;

    std::string anchor_left = sequence.substr(0, g_kmer_length);
    std::string::size_type start = node_set_[p].sequence.find(anchor_left);
    if (start == std::string::npos) {
      while ((int)anchor_left.length() > g_min_anchor_length) {
        anchor_left = anchor_left.substr(0, anchor_left.length()-1);
        start = node_set_[p].sequence.find(anchor_left);
        if (start != std::string::npos) 
          break;
      }
      if (start == std::string::npos) {
        if (g_debug)
          std::cerr << "Can not find start kmer(" << anchor_left << ") when add bubble!" << std::endl;
        return -1;
      }
    }

    std::string anchor_right = sequence.substr(sequence.length()-g_kmer_length+1);
    node_idx_t q = -1;
    while ((int)anchor_right.length() > g_min_anchor_length) {
      q = find_node_index(anchor_right);
      if (q >= 0) break;
      anchor_right = anchor_right.substr(1);
    }
    if (q == -1) {
      restore_kmers(sequence);
      return -2; 
    }

    if (g_debug)
       cerr << "add_bubble..., p = " << p << ", seq = " << sequence << endl;

    int anchor_left_length = anchor_left.length(); // maybe less than g_kmer_length
    int anchor_right_length = anchor_right.length();
    
    std::string::size_type end =  node_set_[q].sequence.find(anchor_right);

    
    // check coverage of the junctions
    int junc_start1 = static_cast<int>(0.5*anchor_left_length);
    const std::string& junc_kmer1 = sequence.substr(junc_start1, g_kmer_length);	
    int junc_start2 = static_cast<int>(sequence.length())-g_kmer_length-static_cast<int>(0.5*anchor_right_length);
    const std::string& junc_kmer2 = sequence.substr(junc_start2, g_kmer_length);
    int junc_cov1 = get_kmer_count(junc_kmer1);
    int junc_cov2 = get_kmer_count(junc_kmer2);
    if (junc_cov1 < g_min_junction_coverage || junc_cov2 < g_min_junction_coverage) {
      if (g_debug)
        std::cerr << "not enough coverage to be junction." << std::endl;
      return -1;
    }
    

    int length = static_cast<int>(sequence.length()) - anchor_left_length - anchor_right_length;
   
    //cerr << "left=" << anchor_left_length << "\t" << "right =" << anchor_right_length << endl;

    if (p == q) {
      
      if (end <= start ) return -1; 
      int distance = static_cast<int>(end-start)-anchor_left_length;
      if (length + distance <= 4)
         return -1;
      if ((distance == 0 && length < g_kmer_length) || (length == 0 && distance < g_kmer_length))
         return -1;
      if ( distance > 0 && length == distance
        && is_similar(node_set_[p].sequence.substr(start+anchor_left_length, distance),
        sequence.substr(anchor_left_length,length), 'F')
        ) {
          return -1;
      }
      if ( length <= 0 && distance <= 0) {
	return -1;
      } else if ( length < 0) {
	anchor_left_length = anchor_left_length + length;
      } else if (distance < 0) {
	anchor_left_length = anchor_left_length + distance;
      }

      if ((int)start+anchor_left_length <= g_kmer_length) // necessary
	return -1;

      Node node1,node2;
      node1.sequence = node_set_[p].sequence.substr(end);
 
      node_idx_t q1 = add_node(node1);
      node_set_[q1].children = node_set_[p].children;
      node_set_[p].children.clear();
      node_idx_t q2 = -1;
      if ( distance <= 0) {
	node_set_[p].sequence = node_set_[p].sequence.substr(0, start+anchor_left_length);
	node_set_[p].add_child(q1);
      } else {
        // divide the node into 3 nodes
        Node node3;
	//  distance + g_kmer_length - anchor_left_length);
        if (length < 0) 
	  node3.sequence = node_set_[p].sequence.substr(start+anchor_left_length, distance-length);
	else 
	  node3.sequence = node_set_[p].sequence.substr(start+anchor_left_length, distance);

        node_idx_t q3 = add_node(node3);
        node_set_[q3].add_child(q1);

        node_set_[p].sequence = node_set_[p].sequence.substr(0, start+anchor_left_length);
        node_set_[p].add_child(q3);
      }

      if (length > 0) {
        if (distance < 0)  // two exons have same first/last several nucleotides 
          node2.sequence = sequence.substr(anchor_left_length, length-distance);
        else 
          node2.sequence = sequence.substr(anchor_left_length, length);
        q2 = add_node(node2);
        node_set_[q2].add_child(q1);
        node_set_[p].add_child(q2);
      } else {
        node_set_[p].add_child(q1);
      }
      
      //if (g_debug)
      //  std::cerr << describe_graph() << std::endl;
      return q2;

    }  // end of case p = q
    else {    // the case p != q

      //cerr << "p != q, q = " << q << ", q size =" <<  node_set_[q].sequence.length() <<endl;
      std::set<node_idx_t> checked;
      if (q == 0 || is_circle(p, q, checked)) 
	return -1; 
      if (node_set_[p].is_child(q) && length > 0 
	&& static_cast<int>(node_set_[p].sequence.length()-start+end) < g_kmer_length) {
        // notice end and start is not in the same node
	return -1;
      }
      Node node1,node2;
      node_idx_t q1 = -1;
      if (length < 0) 
	anchor_left_length = anchor_left_length + length;
      
      if ((int)start+anchor_left_length <= g_kmer_length) //necessary
        return -1;

      //cerr << "anchor_left_length :" << anchor_left_length << endl;
      if (start+anchor_left_length < node_set_[p].sequence.length()) { // does not need split p if equal
        node1.sequence = node_set_[p].sequence.substr(start+anchor_left_length);
        node_set_[p].sequence = node_set_[p].sequence.substr(0,start+anchor_left_length);
        q1 = add_node(node1);
        node_set_[q1].children = node_set_[p].children;
        node_set_[p].children.clear();
        node_set_[p].add_child(q1);
      }
      node_idx_t q2 = -1;
      if (length > 0) {
        node2.sequence = sequence.substr(anchor_left_length, length);
        q2 = add_node(node2);
        node_set_[p].add_child(q2);
      } else {
        q2 = p;
      }
      if (end == 0) {
        node_set_[q2].add_child(q);
      } else {
        // split node q into two nodes
        Node node3;
        assert(end <= node_set_[q].sequence.length());
        node3.sequence = node_set_[q].sequence.substr(end);
        node_idx_t q3 = add_node(node3);
        node_set_[q3].children = node_set_[q].children;
        node_set_[q].children.clear();
        node_set_[q].add_child(q3);
        node_set_[q2].add_child(q3);
        node_set_[q].sequence = node_set_[q].sequence.substr(0, end);
      }
      // return the new bubble node 
      if (q2 == p) {
        return -1;
      } else {
        return q2;
      }
    } // end-of-else
  }


  node_idx_t add_reverse_branch(const node_idx_t p, const std::string& str, KmerMap<T>& kmer_map, std::vector<std::string>& data) {

    int str_length = static_cast<int>(str.length())-g_kmer_length;
    if (str_length >= g_min_exon_length) {
 
      const std::string& end_kmer = str.substr(str.length()-g_kmer_length, g_kmer_length);
      std::string::size_type start = node_set_[p].sequence.find(end_kmer);
      if (start == std::string::npos) {
        restore_kmers(str.substr(0, str.length()-1));
        return -1;
      }
      if (g_debug)
        cerr << "add_reverse_branch : " << str << endl;
      // use pair info to solve short repeat
      if (g_is_paired_end) {
        size_t count = get_kmer_count(end_kmer);
        bool is_branch = check_reverse_branch_with_pair_info(kmer_map, data, str, count);
        if (!is_branch) {
          restore_kmers(str.substr(0, str.length()-1));
          return -1;
        }
      }

      Node node1,node2;
      node1.sequence = node_set_[p].sequence.substr(start);
      node2.sequence = str.substr(0, str_length);

      if (start == 0) {  // possible if not node 0
        node_idx_t q2 = add_node(node2);
        node_set_[q2].add_child(p);  // node1 is p 
        reverse_branches.insert(q2);
        return q2;
      }

      if (is_similar(node_set_[p].sequence.substr(0,start), node2.sequence,'R')) { 
        if (str_length > (int)start && p == 0) {
          node_set_[p].sequence = str.substr(0, str_length) + node_set_[p].sequence.substr(start);
          return -2;
        } else if ((int)start > g_kmer_length) {
          return -1;
        }
      }

      if (p == 0 && (int)start < 2 * g_kmer_length) {
        node_set_[p].sequence = str.substr(0, str_length) + node_set_[p].sequence.substr(start);
        return -2;
      }

      // check junction
      if (!g_is_paired_end) {
        const std::string& junc_kmer = str.substr(str.length()-static_cast<int>(1.5*g_kmer_length),g_kmer_length);
        if ((int)get_kmer_count(junc_kmer) < g_min_junction_coverage) {
           if (str_length > g_min_transcript_length)
             restore_kmers(str.substr(0, str.length()-1));
           return -1;
        }
      }
      // add a reverse branch
      node_set_[p].sequence = node_set_[p].sequence.substr(0, start);
      node_idx_t q1 = add_node(node1);
      node_set_[q1].children = node_set_[p].children;
      node_set_[p].children.clear();
      node_set_[p].add_child(q1);
      node_idx_t q2 = add_node(node2);
      node_set_[q2].add_child(q1);
      reverse_branches.insert(q2);

      //cerr << describe_graph() << endl;
      return q2;

    } else {
      return -1;
    }
  }

  node_idx_t add_reverse_bubble(const node_idx_t p, const std::string& str, KmerMap<T>& kmer_map) {

    if ((int)str.length() < 2 * g_min_anchor_length)
       return -1;

    const std::string& anchor_right = str.substr(str.length()-g_kmer_length, g_kmer_length);
    std::string::size_type start = node_set_[p].sequence.find(anchor_right);
    if (start == std::string::npos)  // possible if overlap 
      return -1;

    std::string anchor_left = str.substr(0, g_kmer_length-1);
    node_idx_t q = -1;
    while ((int)anchor_left.length() > g_min_anchor_length) {
      q = find_node_index(anchor_left);
      if (q >= 0) break;
      anchor_left = anchor_left.substr(1);
    }

    if (q < 0 || p == q || has_path(p, q))
      return -1;
    
    int anchor_length = anchor_left.length();

    std::string::size_type r_end = node_set_[q].sequence.find(anchor_left);
    if (anchor_length + (int)r_end == (int)node_set_[q].sequence.length()) {
      //cerr << "add bubble in reverse check!" << endl;  
      int length = (int)str.length() - anchor_length - g_kmer_length;

      if (length > 0 && start > 0 && length + start < 10)
         return -1;

      Node node1;
      node_idx_t q1 = -1;
      if (length > 0) {
        node1.sequence = str.substr(anchor_length,length);
        q1 = add_node(node1);
        reverse_branches.insert(q1);
        node_set_[q].add_child(q1);
      } else {
        q1 = q;
      }

      if (start == 0) {
        node_set_[q1].add_child(p);
      } else {
        Node node2;
        node2.sequence = node_set_[p].sequence.substr(start);
        node_idx_t q2 = add_node(node2);
        node_set_[q1].add_child(q2);
        node_set_[q2].children = node_set_[p].children;
        node_set_[p].children.clear();
        node_set_[p].add_child(q2);
        node_set_[p].sequence = node_set_[p].sequence.substr(0, start);
      }

      if (q1 != q)
        return q1;
    }
   
    return -1;
  }

  void reverse_check_and_extend(KmerMap<T>& kmer_map, std::vector<std::string>& data, const node_idx_t p) {

    //std::cerr << "reverse_check_and_extend:" << p << std::endl;
    std::vector<kmer_int_type_t> bifurcation_points;
    std::vector<kmer_occurence_pair_t> candidates;
    int check_size = static_cast<int>(node_set_[p].sequence.length()) - g_kmer_length;
    if (check_size < 0)
      return ; 
    for (int i = 0; i < check_size; ++i) {
      const std::string& kmer = node_set_[p].sequence.substr(i,g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);
      kmer_map.get_reverse_candidates(intval, candidates);
      if (candidates.size() <= 1)
        continue;
      for (int j = 0; j < (int)candidates.size(); ++j) {
        if (!is_used(candidates[j].first)) {   
	  bifurcation_points.push_back(intval);
          break;
        }
      }
    }

    if (bifurcation_points.empty())
      return;
    
    for (int i = bifurcation_points.size(); i > 0; --i) {

      kmer_int_type_t intval = bifurcation_points[i-1];
      const std::string& str = reverse_extend(kmer_map, intval);

      const std::string& r_endkmer = str.substr(0, g_kmer_length);
      kmer_int_type_t r_end_val = kmer_to_intval(r_endkmer);
      kmer_map.get_reverse_candidates(r_end_val, candidates);

      if (candidates.size() > 0) {  // add bubble (possible)

        node_idx_t r = add_reverse_bubble(p, str, kmer_map); 
        if (r > 0)
          reverse_branches.insert(r); 
      } else {    // add branch
      
        node_idx_t r = add_reverse_branch(p, str, kmer_map, data);
        if (r > 0)
          reverse_branches.insert(r);
      } //else
    } // for
  }

  void reverse_check_and_extend(KmerMap<T>& kmer_map, std::vector<std::string>& data) {

     int cur_size = size_;

     for (int i = 0; i < cur_size; ++i) {

       if (forward_branches.find(i) != forward_branches.end())
         continue;

       reverse_check_and_extend(kmer_map, data, i);
     }
  }

  void forward_check_and_extend(KmerMap<T>& kmer_map, std::vector<std::string>& data, const node_idx_t p) {

    //cerr << "forward_check_and_extend..." << endl;
    std::vector<kmer_int_type_t> bifurcation_points;
    std::vector<kmer_occurence_pair_t> candidates;
 
    int check_size = static_cast<int>(node_set_[p].sequence.length())-g_kmer_length;
    if (check_size < 0) 
      return; 

    for (int i = 0; i < check_size; ++i) {
      const std::string& kmer = node_set_[p].sequence.substr(i,g_kmer_length);
      kmer_int_type_t intval = kmer_to_intval(kmer);

      kmer_map.get_forward_candidates(intval, candidates);
      if (candidates.size() <= 1) 
        continue;

      //if (compute_entropy(intval, g_kmer_length) < g_min_seed_entropy)
      // continue;

      for (size_t j = 0; j < candidates.size(); ++j) {
        if (!is_used(candidates[j].first)) {
          bifurcation_points.push_back(intval);
          break;
        }
      }
    }

    grow_and_branch(kmer_map, data, p, bifurcation_points);
  }


  void set_parents() {

    for (size_t i = 0; i < size_; ++i) {
      if (!node_set_[i].parents.empty())
        node_set_[i].parents.clear();
    }

    // reset parents of each node 
    for (size_t i = 0; i < size_; ++i) {
      std::vector<node_idx_t>::const_iterator it;
      for (it = node_set_[i].children.begin(); 
	it != node_set_[i].children.end(); ++it) {
          node_set_[*it].add_parent(i);
      }
    }
  }

  void set_reads(KmerMap<info_t>& kmer_map, std::vector<string>& data) {
    ;
  }

  void set_reads(KmerMap<info_with_read_set_t>& kmer_map, std::vector<string>& data) {

    if (kmer_map.empty() || data.empty())
      return;

    size_t max_read_id = data.size();
    if (g_is_paired_end)   // paired
      max_read_id = data.size()/2;

    std::set<size_t> reads;    
    std::map<kmer_int_type_t, size_t>::iterator it;
    for (it = used_kmers_.begin(); it != used_kmers_.end(); ++it) {
      kmer_int_type_t kmer_val = it->first;
      // get the reads 
      const std::vector<size_t>& tempset = (kmer_map[kmer_val]).get_readset();
      if (tempset.empty())
        continue;
      for (size_t i = 0; i < tempset.size(); ++i) {
        if (tempset[i] < max_read_id) // it always meet this condiction for single reads
          reads.insert(tempset[i]);
        else
          reads.insert(tempset[i]-max_read_id);
      }
    } 

    std::set<size_t>::iterator its;
    for(its = reads.begin(); its != reads.end(); ++its) {
      reads_.push_back(data[*its]);
      if (g_is_paired_end)
	reads_.push_back(data[*its+max_read_id]);
    }

  }

  void set_reads(KmerMap<T>& kmer_map, const std::string& seq, std::set<size_t>& reads) {

    if ((int)seq.length() < g_kmer_length) 
      return;

    for (int j = 0; j <= (int)seq.length() - g_kmer_length; ++j) {

        const std::string& kmer = seq.substr(j, g_kmer_length);
        kmer_int_type_t intval = kmer_to_intval(kmer);

        // get reads containing this kmer
        const std::vector<size_t>& tempset = (kmer_map[intval]).get_readset();
        if (tempset.empty()) 
          continue;

        for(size_t k = 0; k < tempset.size(); ++k) 
          reads.insert(tempset[k]);
    }
  }

  void get_reads_from_graph() {
    
    if (reads_.empty()) {
      std::cerr << "Error, No reads!" << std::endl;
      return;
    }

    if (g_is_paired_end) {

      size_t max_read_id = reads_.size()/2;

      for (size_t i = 0; i < max_read_id; ++i) 
        std::cout << ">r" << i << "_1" << std::endl << reads_[2*i] << std::endl;

      for (size_t i = 0; i < max_read_id; ++i)
        std::cout << ">r" << i << "_2" << std::endl << reads_[2*i+1] << std::endl;
      
    } else {

      for (size_t i = 0; i < reads_.size(); ++i) 
        std::cout << ">r" << i << std::endl << reads_[i] << std::endl;
    }
    
  }

 
  int reads_support_contig(KmerMap<T>& kmer_map, std::vector<std::string>& data, const std::string& seq) {

    std::set<size_t> reads;
    int num_aligned = 0;

    // make sure neither kmer_map or data are empty 
    if (seq.length() > g_kmer_length) { 
      for (size_t i = 0; i <= seq.length()-g_kmer_length; ++i) {
        const std::string& kmer = seq.substr(i, g_kmer_length);
        kmer_int_type_t intval = kmer_to_intval(kmer);
        // get reads containing this kmer
        const std::vector<size_t>& tempset = (kmer_map[intval]).get_readset();
        for(size_t j = 0; j < tempset.size(); ++j)
          reads.insert(tempset[j]);
      }
      // align each read to reference contig
      std::set<size_t>::iterator its = reads.begin();
      for (; its != reads.end(); ++its) {
	const std::string& read = data[*its];
        // they are comparible if :
	//   	read :        -------
        //  	contig:  ----------
        // or
        //	read :	   -----
        //      contig:  ------------- 
        // or
        //      read :  -----
        //      contig:    ------------
        if (compatible(seq,read)) 
	  num_aligned++;
      }
    }

    return num_aligned;
  }

  bool has_support(KmerMap<T>& kmer_map, std::vector<std::string>& data, const std::string& seq) {

    if (kmer_map.empty() || data.empty()) 
      return true;

    std::set<size_t> reads;
    set_reads(kmer_map, seq, reads);

    if (reads.empty()) {  // exclude the reason of short seq
      if ((int)seq.length() < g_kmer_length + 5);
        return true;
    }

    int num_aligned = 0;
    // align each read to reference contig
    std::set<size_t>::iterator its = reads.begin();
    for (; its != reads.end(); ++its) {
        const std::string& read = data[*its];
        // they are comparible if :
        //      read :        -------
        //      contig:  ----------
        // or
        //      read :     -----
        //      contig:  ------------- 
        // or
        //      read :  -----
        //      contig:    ------------
        if (compatible(seq,read))
          num_aligned++;
        if (num_aligned >= g_min_reads_span_junction)
	  break;
    }

    return (num_aligned >= g_min_reads_span_junction);
    
  }


  void map_reads_to_graph(KmerMap<T>& kmer_map, std::vector<std::set<size_t> >& reads_mapped_to_node) {
    
    if (size_ <= 1 || kmer_map.empty()) 
      return;

    reads_mapped_to_node.resize(size_);
    for (size_t i = 0; i < size_; ++i) {

      std::string seq = node_set_[i].sequence;
      //if (seq.length() < g_kmer_length) {  // try best to use pair information 
	if (node_set_[i].parents.size() == 1) {
	  node_idx_t p = node_set_[i].parents[0]; 
	  const std::string & parent = node_set_[p].sequence;
          if ((int)parent.length() < g_kmer_length) {
	    seq = parent + seq;
	  } else{
	    seq = parent.substr(parent.length()-g_kmer_length+1) + seq;
          }
        }

        if (node_set_[i].children.size() == 1) {
	  node_idx_t c = node_set_[i].children[0];
          const std::string& child = node_set_[c].sequence;
          seq += child.substr(0, g_kmer_length-1);
        }

        // if seq is a long branch, mapping reads to a part is enough.
        if (node_set_[i].children.empty()) {
	  if ((int)seq.length() > g_max_pair_gap_length) 
	    seq = seq.substr(0, g_max_pair_gap_length);
        } else if (node_set_[i].parents.empty()) {
          if ((int)seq.length() > g_max_pair_gap_length)
            seq = seq.substr(seq.length()-g_max_pair_gap_length, g_max_pair_gap_length);
        }
      //}

      set_reads(kmer_map, seq, reads_mapped_to_node[i]);
      //cerr << "node " << i << " :" << reads_mapped_to_node[i].size() << endl;
    } 

  }

  // given some reads, get their mate distribution 
  // mode :  0 -> only find right mate(/2) ; 1-> only find left mate (/1)
  void get_mate_distribution(std::set<size_t>& reads, 
			     std::map<node_idx_t, size_t>& distribution,
  			     std::vector<std::set<size_t> >& reads_mapped_to_node, 
                             int mode = 2) {

    // only for paired data 
    if (!g_is_paired_end || reads.empty())
      return;

    if (g_double_stranded_mode) 
      mode = 2;

    if (mode == 0) {

       // only find right mate distribution
        std::set<size_t>::iterator it = reads.begin();
        for ( ; it != reads.end(); ++it) {
	  size_t left_id, right_id;
          if ( (*it) % 2 == 0) {
            left_id = *it;
            right_id = left_id + 1;
          } else {
            continue;
          }

          for (size_t i = 0; i < size_; ++i) {
            if (reads_mapped_to_node[i].find(right_id) != reads_mapped_to_node[i].end()) {
              if (distribution.find(i) == distribution.end()) 
                distribution[i] = 1;
              else
		distribution[i]++;
            }
          }
        }
     } else if (mode == 1) { 

        // only find left mate distribution
        std::set<size_t>::iterator it = reads.begin();
        for ( ; it != reads.end(); ++it) {
          size_t left_id, right_id;
          if ((*it) % 2 == 1) {
            right_id = *it;
            left_id = right_id - 1;
          } else {
            continue;
          }

          for (size_t i = 0; i < size_; ++i) {
            if (reads_mapped_to_node[i].find(left_id) != reads_mapped_to_node[i].end()) {
              if (distribution.find(i) == distribution.end())
                distribution[i] = 1;
              else
                distribution[i]++;
            }
          }
        }
     } else if (mode == 2) {

        // both
        std::set<size_t>::iterator it = reads.begin();
        for ( ; it != reads.end(); ++it) {
          size_t read_id = *it;
          size_t mate_id;
          if (read_id % 2 == 0) {
            mate_id = read_id+1;
          } else {
            mate_id = read_id-1;
          }

	  for (size_t i = 0; i < size_; ++i) {
            if (reads_mapped_to_node[i].find(mate_id) != reads_mapped_to_node[i].end()) {
              if (distribution.find(i) == distribution.end())
                distribution[i] = 1;
              else
                distribution[i]++;
            }
          }
        }
    }

    // prune error
    float max = 0;
    std::map<node_idx_t, size_t>::iterator it;
    for (it = distribution.begin(); it != distribution.end(); ++it) {
      if (static_cast<float>(it->second) > max)
	 max = it->second;
    }

    for (it = distribution.begin(); it != distribution.end(); ) {
      if (static_cast<float>(it->second)/max < g_min_ratio_non_error || it->second < 2) 
	distribution.erase(it++);
      else
	++it;
    }

  }

  bool build(KmerMap<T>& kmer_map, kmer_int_type_t seed_val, std::vector<std::string>& data) {

    if (!get_trunk(kmer_map, data, seed_val))
      return false;
   
    // forward check and extend the trunk (its id is 0)
    forward_check_and_extend(kmer_map, data, 0);

    //reverse check and extend the trunk 
    // now, the trunk may have been splited into several nodes
    reverse_check_and_extend(kmer_map, data);
    
    check_edges(kmer_map, data);

    //check if graph is big enough 
    if ( (size_ == 1 && (int)node_set_[0].sequence.length() > g_min_transcript_length)
       || (size_ > 1 && (int)get_total_amount_of_kmers() > g_min_kmers_per_graph)) {
    //if ((int) get_total_amount_of_kmers() > g_min_kmers_per_graph) {
 
      // set parents
      set_parents();

      // extend using pair info
      check_completeness_of_graph(kmer_map, data);

      // collect reads  
      if (size_ > 1) 
        set_reads(kmer_map, data);

      //std::cerr << "finish a graph, size = " << size_ << std::endl;
      //if (!reads_.empty()) 
      //  std::cerr << "involved reads: " << reads_.size() << std::endl;
      return true;
    } else {
      return false;
    }
  }

  size_t get_total_amount_of_kmers() {

    size_t kmer_count = 0;
    for(unsigned int i = 0; i < size_; ++i) 
      kmer_count += node_set_[i].sequence.length();

    return (kmer_count-g_kmer_length);
  }
 

  // compute sequence of edge from source to target
  std::string get_edge_sequence(node_idx_t source, node_idx_t target) {

    // get left anchor sequence
    int length = g_kmer_length-1;
    int start = (int)node_set_[source].sequence.length() > length ?
       static_cast<int>(node_set_[source].sequence.length())-length : 0;
    std::string edge = node_set_[source].sequence.substr(start); 
    while ((int)edge.length() < length) {
      if (node_set_[source].parents.size() >= 1) {
	int parent = node_set_[source].parents[0];
	int remain = length - edge.length();
        int start_p = (int)node_set_[parent].sequence.length() > remain ? 
	  node_set_[parent].sequence.length()-remain : 0;
        edge = node_set_[parent].sequence.substr(start_p) + edge;
	source = parent;
      } else {
	break;
      }
    }
    
    // get right anchor sequence
    if ((int)node_set_[target].sequence.length() < length) {
      edge = edge + node_set_[target].sequence;
    } else {
      edge += node_set_[target].sequence.substr(0,length);
    }

    while ((int)edge.length() < g_kmer_length) {
      if (node_set_[target].children.size() >= 1) {
        int child = node_set_[target].children[0];
        int remain = length - node_set_[target].sequence.length();
        // no problem even if remain > node_set_[child].sequence.length() 
        edge += node_set_[child].sequence.substr(0,remain);
        target = child;
      } else {
	break;
      }
    }
    // refine
    int edge_len = edge.length();
    if (edge_len > g_kmer_length + 8) {
       edge = edge.substr(2, edge_len-4);   
    } else if (edge_len > g_kmer_length + 4) {
      edge = edge.substr(1, edge_len-2);
    }

    return edge;
  }


  // find all edges and compute their average coverages
  void get_coverage_of_edges(std::map<pair_t, double>& edge_cov_map, std::vector<pair_t>& edges) {

    //std::vector<pair_t> edges;
    for (int i = 0; i < (int)get_size(); ++i) {
      if (node_set_[i].children.empty())
        continue;
      for (int j = 0; j < (int)node_set_[i].children.size(); ++j)
        edges.push_back(pair_t(i,node_set_[i].children[j]));
    }

    // coverage for each edge
    for (int i = 0; i < (int)edges.size(); ++i) {
      int source = edges[i].first;
      int target = edges[i].second;
      // compute sequence of edge from source to target
      const std::string& edge = get_edge_sequence(source, target);
      // compute coverage
      edge_cov_map[edges[i]] = compute_coverage(edge);
      if (g_debug) {
        cerr << "edge(" << source << "-> " << target << ") :" << edge << "\t" << show_coverage_detail(edge) << endl;
      }
    }
  }

  
  void get_coverage_of_edges(std::map<pair_t, double>& edge_cov_map) {

    std::vector<pair_t> edges;
    for (int i = 0; i < (int)get_size(); ++i) {
      if (node_set_[i].children.empty())
        continue;
      for (int j = 0; j < (int)node_set_[i].children.size(); ++j)
        edges.push_back(pair_t(i,node_set_[i].children[j]));
    }

    // coverage for each edge
    for (int i = 0; i < (int)edges.size(); ++i) {
      int source = edges[i].first;
      int target = edges[i].second;
      // compute sequence of edge from source to target
      const std::string& edge = get_edge_sequence(source, target);
      // compute coverage
      edge_cov_map[edges[i]] = compute_coverage(edge);
    }
  }
    

  void get_coverage_of_nodes(std::map<node_idx_t, double>& node_cov_map) {

    for (size_t i = 0; i < get_size(); ++i) {

      const std::string & seq = node_set_[i].sequence;

      // extend it if sequence is too short
      if ((int)seq.length() < g_kmer_length + 10) {
        // only one parent
        std::string sequence;
        if (node_set_[i].parents.size() == 1) {
          node_idx_t j = node_set_[i].parents[0];
          const std::string& parent = node_set_[j].sequence;
          if ((int)parent.length() < g_kmer_length) {
            sequence = parent + seq;
          } else {
            // only last (k-1) bases
            sequence = parent.substr(parent.length()-g_kmer_length+1) + seq;
          }
        }

        // only one child
        if (node_set_[i].children.size() == 1) {
          node_idx_t j = node_set_[i].children[0];
          const std::string& child = node_set_[j].sequence;
          sequence = seq + child.substr(0, g_kmer_length-1);
        }
        node_cov_map[i] = compute_coverage(sequence);
      } else {
         node_cov_map[i] = compute_coverage(seq);
      }
    }

  }

  void compact_graph() { 

     if (size_ <= 100) return;
     std::vector<node_idx_t> nodes;
     for (int i = 0; i < (int)get_size(); ++i) {
       if (node_set_[i].sequence.length() < 5) 
	 nodes.push_back(i);  
     }

     // compact
     for (size_t i = 0; i < nodes.size(); ++i) {

       int n = nodes[i];

       // forward extention
       if (node_set_[n].parents.size() <= 1 && node_set_[n].children.size() > 1) {

         bool compact = true;
         for (size_t j = 0; j < node_set_[n].children.size(); ++j) {
	   if (node_set_[node_set_[n].children[j]].parents.size() != 1) 
	     compact = false;
	 }

         if (compact) {
           std::vector<node_idx_t> children = node_set_[n].children; // need this copy
	   for (size_t j = 0; j < children.size(); ++j) {
             int c = children[j];
	     node_set_[c].sequence = node_set_[n].sequence + node_set_[c].sequence;
             node_set_[c].delete_parent(n);
	     node_set_[n].delete_child(c);
           }
           if (node_set_[n].parents.size() == 1) {
	     int p = node_set_[n].parents[0];
             node_set_[n].delete_parent(p);
             node_set_[p].delete_child(n);
             for (size_t j = 0; j < children.size(); ++j) {
               int c = children[j];
               node_set_[c].add_parent(p);
               node_set_[p].add_child(c);
             }
           }
         }
       } else if (node_set_[n].parents.size() > 1 && node_set_[n].children.size() <= 1) {

         bool compact = true;
         for (size_t j = 0; j < node_set_[n].parents.size(); ++j) {
           if (node_set_[node_set_[n].parents[j]].children.size() != 1)
             compact = false;
         }

         if (compact) {
           std::vector<node_idx_t> parents = node_set_[n].parents;
           for (size_t j = 0; j < parents.size(); ++j) {
             int p = parents[j];
             node_set_[p].sequence = node_set_[p].sequence + node_set_[n].sequence;
             node_set_[p].delete_child(n);
             node_set_[n].delete_parent(p);
           }
           if (node_set_[n].children.size() == 1) {
             int c = node_set_[n].children[0];
	     node_set_[c].delete_parent(n);
             node_set_[n].delete_child(c);
             for (size_t j = 0; j < parents.size(); ++j) {
               int p = parents[j];
               node_set_[c].add_parent(p);
               node_set_[p].add_child(c);
             }
           }
         }
       } //else if

     }
  }

  void forward_extend_by_pair_info(KmerMap<T>& kmer_map, std::vector<std::string>& data, node_idx_t c) {

    size_t max_read_id = data.size() / 2;

    while (1) {

      int length = node_set_[c].sequence.length();
      int start_pos = length > 2 * g_kmer_length ? (length - 2 * g_kmer_length) : 0;

      const std::string& check_seq = node_set_[c].sequence.substr(start_pos);
      std::set<size_t> reads;
      set_reads(kmer_map, check_seq, reads);

      if (reads.empty()) return;

      std::set<size_t>::iterator its;
      bool extend_flag = false;
      kmer_int_type_t extend_val;
      std::string extend_str;
      for (its = reads.begin(); its != reads.end(); ++its) {

        size_t mate_id = 0;
        bool need_rev_comp = false;
        if (g_double_stranded_mode) {
          if (compatible(check_seq, data[*its])) {
            need_rev_comp = true;
          } else {
            const std::string & rev_comp = revcomp(data[*its]);
            if (!compatible(check_seq, rev_comp))
              continue;
          }
          if (*its < max_read_id) {
            mate_id = *its + max_read_id;
          } else {
            mate_id = *its - max_read_id;
          }
        } else {
          if (g_fr_strand == 1) {   // --2--> .......... <--1--
            if (*its < max_read_id || (!compatible(check_seq, data[*its])))
              continue;
            mate_id = *its - max_read_id;
          } else if (g_fr_strand == 2) {  // --1--> .......... <--2--
             if (*its >= max_read_id || (!compatible(check_seq, data[*its])))
              continue;
            mate_id = *its + max_read_id;
          } 
        }
        const std::string& mate_read = data[mate_id];

        // not contain N
        if (mate_read.find("N") != std::string::npos)
          continue;

        kmer_int_type_t max_kmer = 0ll;
        int max = 0;
        for (int j = 0; j <= (int)mate_read.length()-g_kmer_length; ++j) {
          kmer_int_type_t intval = kmer_to_intval(mate_read.substr(j, g_kmer_length));
          if (is_used(intval)) {
            max = 0;
            break;
          } else if ((int)kmer_map.get_kmer_count(intval) > max) {
            max_kmer = intval;
            max = kmer_map.get_kmer_count(max_kmer);
          }
        }

        if (max >= g_min_seed_coverage && compute_entropy(max_kmer,g_kmer_length) >= g_min_seed_entropy) {
          kmer_int_type_t stop_kmer = kmer_to_intval(node_set_[c].sequence.substr(length-g_kmer_length));
          if (g_double_stranded_mode && need_rev_comp) {
            max_kmer = revcomp_val(max_kmer, g_kmer_length);
          }
          std::string str;
          extend_val = get_reverse_end_kmer(kmer_map, max_kmer, str, stop_kmer);
          const std::string & kmer = intval_to_kmer(extend_val, g_kmer_length);
          const std::string & anchor = kmer.substr(0, 5);
          std::string::size_type start = check_seq.find(anchor);
          if (start != std::string::npos) {
            if (is_similar(check_seq.substr(start), kmer, 'F')) {
              node_set_[c].sequence = node_set_[c].sequence.substr(0, start_pos+start) + str;
              extend_flag = true;
            }
          } else {
            if (!g_double_stranded_mode) {
              if (((int)str.length() > g_pair_gap_length - 80) && (str.length() > extend_str.length())) 
                extend_str = str;
            }
          }
          if (extend_flag) {
            add_used_kmers(kmer_map, str);
            break;
          }
        }
      }

      if (extend_flag) {

        const std::string & extend_kmer = node_set_[c].sequence.substr(node_set_[c].sequence.length()-g_kmer_length);
        extend_val =  kmer_to_intval(extend_kmer);
        const std::string& seq = forward_extend(kmer_map, extend_val, true);
        //const std::string& seq = forward_extend_for_completeness(kmer_map, extend_val);
        //cerr << "completeness / forward extend :" << seq.length() << endl;
        node_set_[c].sequence = node_set_[c].sequence + seq.substr(g_kmer_length);
        if ((int)seq.length() < 2 * g_kmer_length)
          return;
      } else {
         if (!g_double_stranded_mode) {
           if (((int)extend_str.length() > g_pair_gap_length - 80) && ((int)extend_str.length() < g_max_pair_gap_length)) {
             node_set_[c].sequence = node_set_[c].sequence + extend_str;
             add_used_kmers(kmer_map, extend_str);
             //cerr <<"add withoust overlap:" << endl<< node_set_[c].sequence << endl << show_coverage_detail(node_set_[c].sequence) << endl;
           }
         }
         return;
      }

    }
  }

  // only used for strand specific pair_end data
  void reverse_extend_by_pair_info(KmerMap<T>& kmer_map, std::vector<std::string>& data, node_idx_t p) {

    size_t max_read_id = data.size() / 2;

    while (1) {

      const std::string& check_seq = node_set_[p].sequence.substr(0, 2*g_kmer_length);
      std::set<size_t> reads;
      set_reads(kmer_map, check_seq, reads);
      if (reads.empty())
        return;

      std::set<size_t> ::iterator its;
      bool extend_flag = false;
      kmer_int_type_t extend_val;
      std::string extend_str;
      for (its = reads.begin(); its != reads.end(); ++its) {

        size_t mate_id = 0;
        bool need_rev_comp = false;
        if (g_double_stranded_mode) {
          if (compatible(check_seq, data[*its])) {
            need_rev_comp = true;
          } else {
            const std::string & rev_comp = revcomp(data[*its]);
            if (!compatible(check_seq, rev_comp))
              continue;
          }
          if (*its < max_read_id) {
            mate_id = *its + max_read_id;
          } else {
            mate_id = *its - max_read_id;
          }
        } else {
          if (g_fr_strand == 1) {   // --2--> .......... <--1--
            if (*its >= max_read_id || (!compatible(check_seq, data[*its])))
              continue;
            mate_id = *its + max_read_id;
          } else if (g_fr_strand == 2) {  // --1--> .......... <--2--
             if (*its < max_read_id || (!compatible(check_seq, data[*its])))
              continue;
            mate_id = *its - max_read_id;
          }
        }
        const std::string& mate_read = data[mate_id];
        // not contain N
        if (mate_read.find("N") != std::string::npos)
          continue;

        kmer_int_type_t max_kmer = 0ll;
        int max = 0;
        for (int j = 0; j <= (int)mate_read.length()-g_kmer_length; ++j) {
          kmer_int_type_t intval = kmer_to_intval(mate_read.substr(j,g_kmer_length));
          if (is_used(intval)) {
            max = 0;
            break;
          } else if ((int)kmer_map.get_kmer_count(intval) > max) {
            max_kmer = intval;
            max = kmer_map.get_kmer_count(intval);
          }
        }

        if (max > g_min_seed_coverage && compute_entropy(max_kmer, g_kmer_length) >= g_min_seed_entropy) {

          kmer_int_type_t stop_kmer = kmer_to_intval(node_set_[p].sequence.substr(0, g_kmer_length));
          if (g_double_stranded_mode && need_rev_comp) {
            max_kmer = revcomp_val(max_kmer, g_kmer_length);
          }
          std::string str;
          extend_val = get_forward_end_kmer(kmer_map, max_kmer, str, stop_kmer);
          const std::string & kmer = intval_to_kmer(extend_val, g_kmer_length);
          const std::string & anchor = kmer.substr(g_kmer_length-5); // last five bases
          std::string::size_type start = check_seq.find(anchor);
          if (start != std::string::npos) {
            if (is_similar(check_seq.substr(0,start+5), kmer, 'R')) {
              node_set_[p].sequence = str + node_set_[p].sequence.substr(start+5);
              extend_flag = true;
            }
          } else {
            if (!g_double_stranded_mode) {
              if (((int)str.length() > g_pair_gap_length - 80) && (str.length() > extend_str.length())) 
                extend_str = str;
            }
          }

          if (extend_flag) {
            add_used_kmers(kmer_map, str);
            break;
          }
        }
      }

      if (extend_flag) {

        const std::string& extend_kmer = node_set_[p].sequence.substr(0, g_kmer_length);
        extend_val = kmer_to_intval(extend_kmer);
        // reverse extend
        const std::string& seq = reverse_extend(kmer_map, extend_val, true);
        //const std::string& seq = reverse_extend_for_completeness(kmer_map, extend_val);
        //cerr << "completeness / reverse extend :" << seq.length() << endl;
        node_set_[p].sequence = seq.substr(0, seq.length()-g_kmer_length) + node_set_[p].sequence;
        if ((int)seq.length() < 2 * g_kmer_length)
          return; // extend fail
      } else {
        if (!g_double_stranded_mode) {
          if (((int)extend_str.length() > g_pair_gap_length - 80) && ((int)extend_str.length() < g_max_pair_gap_length)) {
            node_set_[p].sequence = extend_str + node_set_[p].sequence;
            add_used_kmers(kmer_map, extend_str);
          }
        }
        return;  // can not extend
      }

    }
  }


  void check_completeness_of_graph(KmerMap<T>& kmer_map, std::vector<std::string>& data) {

    if (!g_is_paired_end)  // use pair info
      return;

    std::vector<node_idx_t> nodes_without_parents;  
    std::vector<node_idx_t> nodes_without_children; 
    for (int i = 0; i < (int)size_; ++i) {
      if (node_set_[i].children.empty()) 
	nodes_without_children.push_back(i);
      if (node_set_[i].parents.empty())
        nodes_without_parents.push_back(i);
    }

    // initial exons
    for (size_t i = 0; i < nodes_without_parents.size(); ++i) {

      node_idx_t p = nodes_without_parents[i];
      // reverse extend by pair 
      reverse_extend_by_pair_info(kmer_map, data, p);
    } 

 
    // terminal exons
    for (size_t i = 0; i < nodes_without_children.size(); ++i) {

      node_idx_t c = nodes_without_children[i];
      // forward extend by pair 
      forward_extend_by_pair_info(kmer_map, data, c);
    } 

  }

  void trim_graph(KmerMap<T>& kmer_map, std::vector<std::string>& data) {

    if (size_ <= 1) return;

    time_t start = time(NULL);
    bool has_trimed = false;

    std::vector<pair_t> edges;
    std::map<pair_t, double> edge_cov_map;
    get_coverage_of_edges(edge_cov_map, edges);

    std::map<node_idx_t, double> node_cov_map;
    get_coverage_of_nodes(node_cov_map);
    
    std::map<node_idx_t, double> total_out_cov;
    std::map<node_idx_t, double> total_in_cov;
    for (size_t i = 0; i < edges.size(); ++i) {

      int source = edges[i].first;
      if (total_out_cov.find(source) == total_out_cov.end())
        total_out_cov[source] = edge_cov_map[edges[i]];
      else 
        total_out_cov[source] += edge_cov_map[edges[i]];

      int target = edges[i].second;
      if (total_in_cov.find(target) == total_in_cov.end())
        total_in_cov[target] = edge_cov_map[edges[i]];
      else
        total_in_cov[target] += edge_cov_map[edges[i]];
    }

    for (size_t i = 0; i < edges.size(); ++i) {
      
      int source = edges[i].first;
      int target = edges[i].second;

      //if (node_set_[source].parents.empty() || node_set_[target].children.empty()) {
      if ((node_set_[source].parents.empty() && node_set_[source].children.size() == 1) ||
          (node_set_[target].children.empty() && node_set_[target].parents.size() == 1) ||
          (!g_is_paired_end) ) {  // check every edge for single reads

        const std::string& check_edge = get_edge_sequence(source,target);
        if ((int)check_edge.length() < g_kmer_length)
          continue;

        double e_cov = edge_cov_map[edges[i]];
        if (e_cov < 1e-8) { 
          //std::cerr << "graph " << rg_file << ": edge cov = 0" << std::endl;
          if (g_is_paired_end)
            continue;
        }

        double flanking_node_cov = node_cov_map[source] > node_cov_map[target] ?
          node_cov_map[source] : node_cov_map[target];

        if ( (!has_support(kmer_map, data, check_edge)) ||
           e_cov < g_min_junction_coverage ||
           e_cov < g_min_ratio_welds * flanking_node_cov ||
           e_cov < g_min_ratio_branch * total_out_cov[source] ||
           e_cov < g_min_ratio_branch * total_in_cov[target] ||
           (total_in_cov.find(source) != total_in_cov.end() && e_cov < g_min_ratio_in_out * total_in_cov[source]) ||
           (total_out_cov.find(target) != total_out_cov.end() && e_cov < g_min_ratio_in_out * total_out_cov[target])) {
          if (g_debug)
            std::cerr << "delete edge : " << source << "->" << target << std::endl;
          has_trimed = true;
          node_set_[source].delete_child(target);
          node_set_[target].delete_parent(source);
        }

      }
    }

    if (has_trimed)
      recompact_graph();

    time_t end = time(NULL);
    if (g_debug)
      std::cerr << "trim graph : " << (end-start) << " s." << endl;
  }


  void recompact_graph() {

    for (size_t i = 0; i < size_; ++i) {
      node_idx_t p = i;
      if (node_set_[p].children.size() == 1) {
        node_idx_t c = node_set_[p].children[0];
        if (node_set_[c].parents.size() == 1 && node_set_[c].parents[0] == p) {
          node_set_[p].children.clear();
          node_set_[p].children = node_set_[c].children;
          node_set_[p].sequence = node_set_[p].sequence + node_set_[c].sequence;
          // children of node c
          if (!node_set_[c].children.empty()) {
            for (size_t j = 0; j < node_set_[c].children.size(); ++j) {
              node_idx_t cc = node_set_[c].children[j];
              node_set_[cc].delete_parent(c);
              node_set_[cc].add_parent(p);
            }
          }
          node_set_[c].clear();
          i--;  // iteratively check node p
        }

      }
    }
  }

  void process_tips() { // before refine_tips

    std::vector<node_idx_t> nodes_without_parents;
    std::vector<node_idx_t> nodes_without_children;
    for (size_t i = 0; i < size_; ++i) {
      if ((int) node_set_[i].sequence.length() < g_kmer_length) {
        if (node_set_[i].parents.empty() && (!node_set_[i].children.empty()))
          nodes_without_parents.push_back(i);
        else if (node_set_[i].children.empty() && (!node_set_[i].parents.empty()))
          nodes_without_children.push_back(i);
      }
    }

    if (!nodes_without_parents.empty()) {
      for (size_t i = 0; i < nodes_without_parents.size(); ++i) {
        int n = nodes_without_parents[i];
        bool compact = true;
        for (int j = 0; j <  node_set_[n].children.size(); ++j) {
          if (node_set_[node_set_[n].children[j]].parents.size() != 1)
            compact = false;
        }
        if (compact) {
          std::vector<node_idx_t> children = node_set_[n].children;
          for (size_t j = 0; j < children.size(); ++j) {
            int c = children[j];
            node_set_[c].sequence = node_set_[n].sequence + node_set_[c].sequence;
            node_set_[c].delete_parent(n);
            node_set_[n].delete_child(c);
          }
        }
      }
    }

    if (!nodes_without_children.empty()) {
      for (int i = 0; i < nodes_without_children.size(); ++i) {
        int n = nodes_without_children[i];
        bool compact = true;
        for (int j = 0; j <  node_set_[n].parents.size(); ++j) {
          if (node_set_[node_set_[n].parents[j]].parents.size() != 1)
            compact = false;
        }
        if(compact) {
          std::vector<node_idx_t> parents = node_set_[n].parents;
          for (int j = 0; j < parents.size(); ++j) {
            int p = parents[j];
            node_set_[p].sequence = node_set_[p].sequence + node_set_[n].sequence;
            node_set_[n].delete_parent(p);
            node_set_[p].delete_child(n);
          }
        }
      }
    }

  }



  void refine_tips(KmerMap<T>& kmer_map, std::vector<std::string>& data) {

    std::vector<node_idx_t> nodes_without_parents;  
    std::vector<node_idx_t> nodes_without_children; 
    for (size_t i = 0; i < size_; ++i) {
      if (node_set_[i].children.empty()) {
	if (node_set_[i].parents.empty() && (int)node_set_[i].sequence.length() < g_min_transcript_length)
          continue;
        nodes_without_children.push_back(i);
      }
      if (node_set_[i].parents.empty())
        nodes_without_parents.push_back(i);
    }

    int extend = 0;
    for (int i = 0; i < (int)nodes_without_parents.size(); ++i) {

      node_idx_t p = nodes_without_parents[i];
      const string& str = node_set_[p].sequence.substr(0, 2*g_kmer_length);
      std::set<size_t> reads;
      set_reads(kmer_map, str, reads);
      if (reads.empty()) // it maybe empty if str.length() < g_kmer_length
	 continue;
      size_t read_id = 0;
      int extend_len = 0;
      int read_start = -1;
      int node_start = -1;

      std::set<size_t>::iterator it;
      for (it = reads.begin(); it != reads.end(); ++it) {
        const string& read = data[*it];
        if (!compatible(read, str))
          continue;
        for (int j = read.length()-g_kmer_length; j > 0; --j) {
          const std::string& kmer = read.substr(j, g_kmer_length);
          std::string::size_type start = node_set_[p].sequence.find(kmer);
          if (start != std::string::npos) { // get alignment
            if ((int)start < j) {
              int len = j - start;
              if (len > extend_len) {
                 read_id = *it;
                 extend_len = len;
                 read_start = j;
                 node_start = start;
              }
            }
            break;
          } 
        } // for each j
      }

      if (extend_len > 0) {
        node_set_[p].sequence = data[read_id].substr(0,read_start) + node_set_[p].sequence.substr(node_start);
        extend += extend_len;
        if (extend < g_pair_gap_length) 
          i--;   
        else
	  extend = 0;  // clear
      } else {
        extend = 0; 
      }
    }    

    extend = 0;  // reset it to be zero
    for (int i = 0; i < (int)nodes_without_children.size(); ++i) {

      node_idx_t c = nodes_without_children[i];
      int start = (int)node_set_[c].sequence.length() > 2*g_kmer_length ? 
                    (node_set_[c].sequence.length()-2*g_kmer_length) : 0;
      const string& str = node_set_[c].sequence.substr(start, 2*g_kmer_length);
      std::set<size_t> reads;
      set_reads(kmer_map, str, reads);
      if (reads.empty()) continue;
      std::set<size_t>::iterator it;
      size_t read_id = 0;
      int extend_len = 0;
      int read_start = -1;
      int node_start = -1;

      for (it = reads.begin(); it != reads.end(); ++it) {
        const string& read = data[*it];
        if (!compatible(read, str))
          continue;
        for (int j = 0; j < (int)read.length()-g_kmer_length; ++j) {
          const std::string& kmer = read.substr(j, g_kmer_length);
          std::string::size_type start = node_set_[c].sequence.find(kmer);
          if (start != std::string::npos) { // get alignment
            int seq_left = node_set_[c].sequence.length() - start;
            int read_left = read.length()-j;
            if ( seq_left < read_left) {
              int len = read_left - seq_left;
              if (len > extend_len) {
                read_id = *it;
                extend_len = len;
                read_start = j;
                node_start = start;
              }
            }
            break;
          } 
        } // for each j
      }

      if (extend_len > 0) {
        node_set_[c].sequence = node_set_[c].sequence.substr(0, node_start) + data[read_id].substr(read_start);
        extend += extend_len;
        if (extend < g_pair_gap_length && (!g_is_paired_end))
          i--;
        else
          extend = 0;
      } else {
        extend = 0;
      }
    }

  }


  void check_edges(KmerMap<T>& kmer_map, std::vector<std::string>& data) {

    //cerr << "check_edges..." << endl;
    typedef std::pair<int,int> pair_t;
    std::vector<pair_t> edges;
    edges.reserve(size_);

    for (unsigned int i = 0; i < size_; ++i) {
      if (node_set_[i].children.empty()) 
	continue;
      for (unsigned int j = 0; j < node_set_[i].children.size(); ++j) { 
        node_idx_t c = node_set_[i].children[j];
        edges.push_back(pair_t(i,c));
      }
    }

    for (size_t i = 0; i < edges.size(); ++i) {

      int source = edges[i].first;
      int target = edges[i].second;
      if ( (int) node_set_[source].sequence.length() <= g_kmer_length || 
	(int) node_set_[target].sequence.length() <= g_kmer_length )
	continue;

      int source_len = node_set_[source].sequence.length();
      std::string edge = node_set_[source].sequence.substr(source_len+1-g_kmer_length) +
        node_set_[target].sequence.substr(0, g_kmer_length-1);

      for (unsigned int j = 0; j <= edge.length()-g_kmer_length; ++j) {
	  const std::string& kmer = edge.substr(j, g_kmer_length);
          kmer_int_type_t intval = kmer_to_intval(kmer);
          // forward check 
	  std::vector<kmer_occurence_pair_t> candidates;
          kmer_map.get_forward_candidates(intval, candidates);
          if (candidates.size() > 1) {
            std::vector<kmer_int_type_t> bifurcation_points;
	    std::string branch = forward_extend(kmer_map,intval,bifurcation_points);
            if ((int)branch.length() < 2 * g_min_anchor_length)
              continue;

            // Be careful: node_set_[source].sequence.length() can be less than g_kmre_length
	    const std::string& endkmer = branch.substr(branch.length()-g_kmer_length);
            kmer_int_type_t endval = kmer_to_intval(endkmer);
            kmer_map.get_forward_candidates(endval, candidates);
            if (candidates.size() > 0) {
              if ((int)node_set_[source].sequence.length() > g_kmer_length) {
		branch = node_set_[source].sequence.substr(source_len-g_kmer_length) +
                         branch.substr(g_kmer_length-1-j);
                add_bubble(source, branch, kmer_map);
              } 
            } else {
	      if ((int)branch.length() < g_kmer_length + g_min_exon_length)
		continue;
              
              if (g_is_paired_end) {
                size_t count = get_kmer_count(intval);
                bool is_branch = check_forward_branch_with_pair_info(kmer_map, data, branch, count);

        	if (!is_branch) {
          	  restore_kmers(branch.substr(1));
	          continue;
		}
              }

	      Node node1;
	      node1.sequence = branch.substr(g_kmer_length-1-j);
	      if (!is_similar(node1.sequence, node_set_[target].sequence, 'F')) {
                node_idx_t q = add_node(node1);
                node_set_[source].add_child(q);
                if (j + g_kmer_length == edge.length() && node_set_[target].parents.size() > 1) {
                   for (unsigned int k = 0; k < node_set_[target].parents.size(); ++k) {
                     node_idx_t p = node_set_[target].parents[k];
                     if (p != source)   // wrong edges will be deleted by split_graph or trim_graph
                       node_set_[p].add_child(q);
                   }
                }  
              }	      
            }
            if ((int)candidates.size() > 2) { // three junction derive from this exon
              bool redo = false;
	      for (unsigned int k = 0; k < candidates.size(); ++k) {
                if (!is_used(candidates[k].first)) {
		  redo = true; 
		  break;
                }
              }
	      if (redo) {
		--j;
		continue;
              }
            }
	  }

	  // reverse check
	  std::vector<kmer_occurence_pair_t> candidates_reverse;
	  kmer_map.get_reverse_candidates(intval, candidates_reverse);
	  if ((int)candidates_reverse.size() > 1) {
	    const std::string& branch_rev = reverse_extend(kmer_map, intval);
	    if ((int)branch_rev.length() < g_kmer_length + g_min_exon_length)
	      continue; 

            if (g_is_paired_end) {
              size_t count = get_kmer_count(intval);
              bool is_branch = check_reverse_branch_with_pair_info(kmer_map, data, branch_rev, count);
	      if (!is_branch) {
                restore_kmers(branch_rev.substr(0, branch_rev.length()-1));
	        continue;
	      }
            }

	    Node node2;
            node2.sequence = branch_rev.substr(0, branch_rev.length()-j-1);
            if (!is_similar(node2.sequence, node_set_[source].sequence, 'R')) {
              node_idx_t q = add_node(node2);
              if (j == 0 && (int)node_set_[source].children.size() > 1) {
                bool add_flag = false;
                std::vector<size_t> temp = (kmer_map[intval]).get_readset();
                if (!temp.empty()) {
                  std::string anchor;
                  for (unsigned int k = 0; k < temp.size(); ++k) {
                    std::string::size_type start = data[temp[k]].find(intval_to_kmer(intval, g_kmer_length)) ;
                    if (start != std::string::npos && (start+g_kmer_length+5) <= data[temp[k]].length() ) {
                      anchor = data[temp[k]].substr(start+g_kmer_length, 5);
                      break;
                    }
                  }
                  if ((int)anchor.length() > 0) {
                    for (unsigned int k = 0; k < node_set_[source].children.size(); ++k) {
                      if (node_set_[node_set_[source].children[k]].sequence.find(anchor) != std::string::npos) {
                        node_set_[q].add_child(node_set_[source].children[k]);
                        add_flag = true;
                      }
                    }
                  }
                }
                if (!add_flag) {
                  for (unsigned int k = 0; k < node_set_[source].children.size(); ++k)
                    node_set_[q].add_child(node_set_[source].children[k]);
                }
              } else {
                node_set_[q].add_child(target);
              }
            }

	  }          
      } // for j

    } // for i
  }

  bool remove_used_kmers_from_kmer_map(KmerMap<T>& kmer_map) {

    // clear used kmer from hash table
    map<kmer_int_type_t,size_t>::iterator it;
    for (it = used_kmers_.begin(); it != used_kmers_.end(); ++it) {
      
      if (restore_kmers_.find(it->first) != restore_kmers_.end())
	continue;

      kmer_map.remove(it->first);
    }

    return true;
  }

  bool clear(KmerMap<T>& kmer_map) {

    if (!used_kmers_.empty()) 
      remove_used_kmers_from_kmer_map(kmer_map);

    node_set_.clear();
    reads_.clear();
    size_ = 0;
    restore_kmers_.clear();

    return true;
  }

  bool save(const std::string& filename) {

    // create and open an archive for output
    std::ofstream ofs(filename.c_str());
    boost::archive::text_oarchive oa(ofs);

    // write to archive
    oa << node_set_;
    oa << used_kmers_;
    oa << reads_;  
    // archive and stream closed when destructors are called

    return true;
  }
  
  bool save_debug(const std::string& filename){

    std::ofstream ofs(filename.c_str());    
    const std::string graph = describe_graph(); //a readable format
    ofs << graph;

    return true;
  }

  bool load(const std::string& filename) {

    // create and open an archive for input
    std::ifstream ifs(filename.c_str());
    boost::archive::text_iarchive ia(ifs);

    // load from archive
    ia >> node_set_;
    ia >> used_kmers_;
    ia >> reads_; 
    // set size
    size_ = node_set_.size();

    return true;
  }

public:

  std::vector<Node> node_set_;  // node set
  size_t size_;  // the number of nodes in this graph
  std::vector<std::string> reads_;  // reads (single or pairs)
  std::map<kmer_int_type_t, size_t> used_kmers_; // kmers encounted when building graph
  // auxiliary container
  std::set<kmer_int_type_t> restore_kmers_;
  std::set<node_idx_t> forward_branches;
  std::set<node_idx_t> reverse_branches;
};


#endif
