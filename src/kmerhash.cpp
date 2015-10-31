#include "kmerhash.h"
#include <iostream>
#include <algorithm>

// kmerhash.cpp


/*
std::vector<kmer_int_type_t> get_root_kmers(Kmer_hash_map& kmer_hash) {
  std::vector<kmer_int_type_t> roots;
  Kmer_hash_map_iterator it;
  for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
    if(!(it->second).has_prekmer) {
      roots.push_back(it->first);
    }
  }
  return (roots);
}
*/

/*
bool purne_kmer_extensions (Kmer_hash_map& kmer_hash,
  std::vector<Kmer_Occurence_Pair>& candidates,
  int dominant_count,
  float min_ratio_non_error) {
  if (candidates.size() <= 1) {
    return (true);
  } else if (dominant_count < candidates[0].second){
    dominant_count = candidates[0].second;
  }
  for ( unsigned int i = 1; i < candidates.size(); ++i) {
    int candidate_count =  candidates[i].second;
    if ((float)candidate_count/dominant_count < min_ratio_non_error) {
      // prune error kmer
      // Kmer_hash_map_iterator pos = kmer_hash.find(candidates[i].first);
      // kmer_hash.erase(pos);
      delete_kmer(kmer_hash, candidates[i].first);
      candidates.erase(candidates.begin()+i);
    }
  }
  return(true);
}

*/





