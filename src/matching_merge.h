// This file was modified from file matching_merge.cpp in Cufflinks
// The original copyright info is listed below
// Copyright 2009 Cole Trapnell.
// Distributed under the  Distributed under Cufflinks LICENSE
// (See accompanying file LICENSE)

#ifndef MATCHING_MERGE_H
#define MATCHING_MERGE_H

/*
 *  matching_merge.h
 *
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <sstream>

//#include <lemon/topology.h>
#include <lemon/smart_graph.h>
#include <lemon/bipartite_matching.h>

#include "compatibility_graph.h"

using namespace std;
using namespace boost;

typedef lemon::SmartBpUGraph reachability_bp_graph_t;

template<class Matcher>
void make_chains_from_matching(const reachability_bp_graph_t& bp,const Matcher& matcher,
                               std::vector<std::vector<directed_acyclic_graph_node_t> >& chains) {

  // Make chains out of the matching
  reachability_bp_graph_t::ANodeMap<reachability_bp_graph_t::UEdge> matched_a_nodes(bp);
  matcher.aMatching(matched_a_nodes);
  reachability_bp_graph_t::BNodeMap<reachability_bp_graph_t::UEdge> matched_b_nodes(bp);
  matcher.bMatching(matched_b_nodes);

  set<reachability_bp_graph_t::ANode> chain_heads;
  for (reachability_bp_graph_t::ANodeIt i(bp); i!=lemon::INVALID; ++i) {
    int a_id = bp.aNodeId(i);
    reachability_bp_graph_t::ANode a = bp.nodeFromANodeId(a_id);
    if (matched_a_nodes[a] == lemon::INVALID)
      chain_heads.insert(bp.nodeFromANodeId(bp.aNodeId(i)));
  }

  for (set<reachability_bp_graph_t::ANode>::iterator i = chain_heads.begin();
       i != chain_heads.end();
       ++i) {

    std::vector<directed_acyclic_graph_node_t> chain;
    reachability_bp_graph_t::ANode n = *i;
    chain.push_back(bp.aNodeId(*i));

    while (true) {
      //int a_id = bp.aNodeId(n);
      int b_id = bp.bNodeId(n);
      reachability_bp_graph_t::BNode b = bp.nodeFromBNodeId(b_id);
      if (matched_b_nodes[b] == lemon::INVALID) {
        break;
      } else {
        reachability_bp_graph_t::ANode a_match_to_b = bp.source(matched_b_nodes[b]);
        chain.push_back(bp.aNodeId(a_match_to_b));
        n = a_match_to_b;
      }
    }

    chains.push_back(chain);
  }

  assert (chains.size() == chain_heads.size());
}


template <class T>
void find_path(const directed_acyclic_graph_t& bundle_dag,
               const adjacency_list<>& TC,
               const directed_acyclic_graph_node_t& source,
               const directed_acyclic_graph_node_t& target,
               std::vector<directed_acyclic_graph_node_t>& path,
               std::vector<pair_t >& splicing_graph_edges,
               SplicingGraph<T>& splicing_graph,
               std::map<pair_t,double> & edge_cov_map,
               bool find_long_path = false) {

  if (source == target) return;
  bool done = false;
  directed_acyclic_graph_node_t curr = source;

  while(!done) {
    graph_traits<directed_acyclic_graph_t>::adjacency_iterator i, iend;
    for (tie(i,iend) = adjacent_vertices(curr, bundle_dag); i != iend; ++i) {
      directed_acyclic_graph_node_t I = *i;
      pair<adjacency_list<>::edge_descriptor, bool> p;
      p = edge(I, target, TC);
      // p.second will be true if there exists an edge between I and target
      if (p.second) {
        /* add by myself from here */
        if ( (i+1) != iend ) {
          directed_acyclic_graph_node_t J = *(i+1);
          pair<adjacency_list<>::edge_descriptor, bool> q;
          q = edge(J, target, TC);
          if (q.second) {
            // target node of two edge
            int tp = splicing_graph_edges[*i].second;
            int tq = splicing_graph_edges[*(i+1)].second;

            // choose a longer path       
            if (find_long_path && splicing_graph.has_path(tq, tp)) {
              path.push_back(*(i+1));
              curr = *(i+1);
              break;
            }
            
            // choose the most supported path
            if (edge_cov_map[splicing_graph_edges[J]] > edge_cov_map[splicing_graph_edges[I]]) {
              path.push_back(*(i+1));
              curr = *(i+1);
              break;
            }
          }
        } /* to here */
        path.push_back(*i);
        curr = *i;
        break;
      }
      if (*i == target) {
        path.push_back(*i);
        done = true;
        break;
      }
    } // for
  } // while
}


// output the fasta sequence into multiple lines
std::string add_fasta_seq_line_breaks(const std::string& sequence, int interval);
// oupput the path for debugging
std::string describe_path(const std::vector<directed_acyclic_graph_node_t>& path);


template <class T>
void extend_chains_to_paths(const directed_acyclic_graph_t& bundle_dag,
                            std::vector<std::vector<directed_acyclic_graph_node_t> >& chains,
                            adjacency_list<>& transitive_closure,
                            directed_acyclic_graph_node_t source_node,
                            directed_acyclic_graph_node_t sink_node,
                            std::vector<std::vector<directed_acyclic_graph_node_t> >& paths,
                            std::vector<pair_t >& splicing_graph_edges,
                            SplicingGraph<T>& splicing_graph,
                            std::map<pair_t, double>& edge_cov_map) {

  paths.reserve(chains.size());

  // extend each chain to a path
  for(size_t c = 0; c < chains.size(); ++c) {
    std::vector<directed_acyclic_graph_node_t>& chain = chains[c];
    assert(!chain.empty());
    if (source_node == chain[0] || sink_node == chain[0])
      continue;
    // reverse the chain
    reverse(chain.begin(), chain.end());
    std::vector<directed_acyclic_graph_node_t> path;
    // find path source -> chain[0]
    find_path(bundle_dag, transitive_closure, source_node, chain[0], path, splicing_graph_edges, splicing_graph, edge_cov_map);

    for (size_t n = 1; n < chain.size(); ++n) {
      assert (path.back() == chain[n-1]);
      directed_acyclic_graph_node_t last = chain[n-1];
      directed_acyclic_graph_node_t next = chain[n];
      // find path between two adjacent nodes last -> next
      find_path(bundle_dag, transitive_closure, last, next, path, splicing_graph_edges, splicing_graph, edge_cov_map);
    }
    // find path chain[n] -> sink
    find_path(bundle_dag, transitive_closure, chain.back(), sink_node, path, splicing_graph_edges, splicing_graph, edge_cov_map);
    assert(path.back() == sink_node);
    path.pop_back();  // delete sink
    if (g_debug)
      std::cerr << "chain " << c << ": " << describe_path(chain) << std::endl;
    paths.push_back(path);
  }


  // record alternative extentions of chains 
  std::vector<std::vector<directed_acyclic_graph_node_t> > long_paths;
  long_paths.reserve(chains.size());
  for(size_t c = 0; c < chains.size(); ++c) {
    std::vector<directed_acyclic_graph_node_t>& chain = chains[c];
    if (source_node == chain[0] || sink_node == chain[0])
      continue;
    std::vector<directed_acyclic_graph_node_t> path;
    // find path source -> chain[0]
    find_path(bundle_dag, transitive_closure, source_node, chain[0], path, splicing_graph_edges, splicing_graph, edge_cov_map, true);
    for (size_t n = 1; n < chain.size(); ++n) {
      directed_acyclic_graph_node_t last = chain[n-1];
      directed_acyclic_graph_node_t next = chain[n];
      // find path between two adjacent nodes last -> next
      find_path(bundle_dag, transitive_closure, last, next, path, splicing_graph_edges, splicing_graph, edge_cov_map, true);
    }
    // find path chain[n] -> sink
    find_path(bundle_dag, transitive_closure, chain.back(), sink_node, path, splicing_graph_edges, splicing_graph, edge_cov_map, true);
    path.pop_back();  // delete sink
    long_paths.push_back(path);
  }

  for (size_t i = 0; i < long_paths.size(); ++i) {
    if (long_paths[i] == paths[i]) {
      continue;
    } else {
      paths.push_back(long_paths[i]);
    }
  }

  if (g_debug) {
    for (size_t i = 0; i < paths.size(); ++i) 
      std::cerr << "path " << i << ": " << describe_path(paths[i]) << std::endl;
  }

}



// the following function maybe need rewrite
template <class T>
void get_isoforms_from_paths(SplicingGraph<T>& splicing_graph,
                             const std::vector<std::vector<directed_acyclic_graph_node_t> >& paths,
			     std::vector<pair_t >& splicing_graph_edges) {

  // change from compatibility graph path to splicing graph path
  std::vector<std::vector<std::size_t> > node_paths;
  for (size_t i = 0; i < paths.size(); ++i) {
    std::vector<std::size_t> node_path;
    for (size_t j = 0; j < paths[i].size(); ++j) {
      pair_t node_pair = splicing_graph_edges[paths[i][j]];
      if (j == 0)
        node_path.push_back(node_pair.first);
      node_path.push_back(node_pair.second); 
    }
    node_paths.push_back(node_path);
  }

  size_t index = 0; // index of isoforms
  // output the sequences of the isoforms one by one
  for (size_t i = 0; i < node_paths.size(); ++i) {
    std::string isoform;
    std::stringstream path;
    for (size_t j = 0; j < node_paths[i].size(); ++j) {
      path << node_paths[i][j] << ":" 
	<< splicing_graph.node_set_[node_paths[i][j]].sequence.length();
      if (j != node_paths[i].size()-1)
	path << ";";
      isoform += splicing_graph.node_set_[node_paths[i][j]].sequence;
    }
   
    if ((int)isoform.length() > g_min_transcript_length) {
      std::cout << ">" << rg_file << "_seq" << index 
        << " len=" << isoform.length() 
        << " path=[" << path.str() << "]" << std::endl;
      const std::string& seq = add_fasta_seq_line_breaks(isoform, 60);
      std::cout << seq << std::endl;
      ++index;
    }
  }

  // check isolated node in splicing graph
  for (size_t i = 0; i < splicing_graph.get_size(); ++i) {
    if (splicing_graph.node_set_[i].children.empty() && splicing_graph.node_set_[i].parents.empty()) {
      if ((int)splicing_graph.node_set_[i].sequence.length() > g_min_transcript_length) {
        const std::string& isoform = splicing_graph.node_set_[i].sequence;
	std::cout << ">" << rg_file << "_seq" << index << " len=" << isoform.length() 
	  << " path=[" << i << ":" << isoform.length() << "]"<< std::endl;
        const std::string& seq = add_fasta_seq_line_breaks(isoform, 60);
        std::cout << seq << std::endl;
        ++index;
      }
    }
  }

}


template <class T>
bool recover_path(SplicingGraph<T>& splicing_graph,
 		  KmerMap<T>& kmer_map,
		  std::vector<pair_t >& splicing_graph_edges,
		  std::map<pair_t, int>& edge_index,
		  std::vector<std::set<size_t> >& reads_mapped_to_node,
		  std::vector<directed_acyclic_graph_node_t>& path,
                  int inhibit_edge,
                  std::vector<directed_acyclic_graph_node_t>& ref_path,
		  std::map<int, int>& edge_frequency,
	          int mode = 0) {

  // recover it only if it contain nodes used only once
  bool need_recover = false;
  for (size_t i = 0; i < path.size(); ++i) {
    if (edge_frequency[path[i]] <= 1) {
      need_recover = true;
      break;
    }
  }
  //cerr << "inhibit_edge: " << inhibit_edge << endl;
  if ((!need_recover) || path.empty())
    return false; 

  std::set<int> reference;
  for (size_t i = 0; i < ref_path.size(); ++i) { 
    reference.insert(ref_path[i]);
  }

  if (mode == 0) { // find right part

    while (1) {
      int end = path[path.size()-1];
      int new_end = -1;
      //int s = splicing_graph_edges[end].first;
      int t = splicing_graph_edges[end].second;
      //cerr << "end edge: " << end <<"(" << s <<"->" << t <<")"<<endl;
      if (splicing_graph.node_set_[t].children.empty()) 
        break;
      if (splicing_graph.node_set_[t].children.size() == 1) {
	int child = splicing_graph.node_set_[t].children[0];
        std::pair<int,int> edge(t, child);
        new_end = edge_index[edge];
        path.push_back(new_end);
        edge_frequency[new_end]++;
        continue;
      }
      /*
      const string& edge = splicing_graph.get_edge_sequence(s,t);
      std::set<size_t> reads;
      splicing_graph.set_reads(kmer_map, edge, reads);
      std::map<int, size_t> distribution;
      //splicing_graph.get_mate_distribution(reads, distribution, reads_mapped_to_node, 1);
      splicing_graph.get_mate_distribution(reads, distribution, reads_mapped_to_node, 2); // need change
      
      std::map<int, size_t>::iterator its = distribution.begin(); 
      cerr << "distribution:" << endl;
      for ( ; its != distribution.end(); ++its) {
        cerr << its->first << " : " << its->second << endl;
      }     
      */
      std::vector<int> candidates;
      for (int i = 0; i < (int)splicing_graph.node_set_[t].children.size(); ++i) {
        int child = splicing_graph.node_set_[t].children[i];
        pair<int,int> edge(t, child);
        int temp = edge_index[edge];
        if (temp != inhibit_edge)
          candidates.push_back(temp);
      }

      assert(!candidates.empty());

      new_end = candidates[0];
      for (size_t i = 0; i < candidates.size(); ++i) {
        // use as many edges of original path as possible
        if (reference.find(candidates[i]) != reference.end()) {
	  new_end = candidates[i];
          edge_frequency[new_end]++;
	  break;
	}
      } 
     
      assert(new_end >= 0);	
      path.push_back(new_end); 
      edge_frequency[new_end]++;
      //cerr << "current path: " << describe_path(path) << endl;
    }    

  } else if (mode == 1) { // find left part

    while (1) {
      int start = path[0];
      int new_start = -1;
      int s = splicing_graph_edges[start].first;
      //int t = splicing_graph_edges[start].second;
      //cerr << "start edge:" << start << "(" << s <<","<< t <<")" << endl;
      if (splicing_graph.node_set_[s].parents.empty())
        break;
      if (splicing_graph.node_set_[s].parents.size() == 1) {
	int parent = splicing_graph.node_set_[s].parents[0];
        std::pair<int,int> edge(parent, s);
        new_start = edge_index[edge];
        path.insert(path.begin(),new_start);
        edge_frequency[new_start]++;
        continue; 
      }
      /*
      const string& edge = splicing_graph.get_edge_sequence(s,t);
      std::set<size_t> reads;
      splicing_graph.set_reads(kmer_map, edge, reads);
      std::map<int, size_t> distribution;
      //splicing_graph.get_mate_distribution(reads, distribution, reads_mapped_to_node, 0);
      splicing_graph.get_mate_distribution(reads, distribution, reads_mapped_to_node, 2); // need change

      cerr << "distribution:" << endl;
      std::map<int, size_t>::iterator its = distribution.begin();
      for ( ; its != distribution.end(); ++its) {
        cerr << its->first << " : " << its->second << endl;
      }
      */
      std::vector<int> candidates;
      for (size_t i = 0; i < splicing_graph.node_set_[s].parents.size(); ++i) {
        int parent = splicing_graph.node_set_[s].parents[i];
        pair<int,int> edge(parent,s);
        int temp = edge_index[edge];
        if (temp != inhibit_edge) 
          candidates.push_back(temp);
      }
      assert(candidates.size() >= 1); //debug

      new_start = candidates[0];
      for (size_t i = 0; i< candidates.size(); ++i) {
        // use edge of original path
        if (reference.find(candidates[i]) != reference.end()) {
          new_start = candidates[i];
          break;
        }
      }
      // add new node to path
      assert(new_start >= 0);  // debug

      path.insert(path.begin(), new_start);
      edge_frequency[new_start]++;
      //cerr << "current path: " << describe_path(path) << endl;
    } // while
  }
  return true;
}

/*
template <class T>
bool check_node(SplicingGraph<T>& splicing_graph, 
                std::map<int, size_t>& distribution,
                std::vector<int>& node_path, int index) {

  assert(index >= 0);
  int check_length = splicing_graph.node_set_[node_path[index]].sequence.length();
 
  // first, there must exist one parent of node index 
  bool has_parent = false;
  std::vector<int>::iterator it = splicing_graph.node_set_[node_path[index]].parents.begin();
  for ( ; it != splicing_graph.node_set_[node_path[index]].parents.end(); ++it) {
    if ( distribution.find(*it) != distribution.end()) {
      has_parent = true;
      break;
    }
  }
  if (!has_parent) {  // too long or too short
   if (distribution.find(index) != distribution.end())
      return true;
  }
  //bool has_support = false;
  for (int j = index-1; j >= 0; --j) {   //check the ancestors of node index
    if (distribution.find(node_path[j]) != distribution.end()) {
      return true;
    } else {
      check_length += splicing_graph.node_set_[node_path[j]].sequence.length();
      if (check_length > g_max_pair_gap_length)
        break;
    }
  }

  //if (check_length < g_max_pair_gap_length) 
  //  return true; // too short for finding paired read
  //else
  return false;
}
*/

template <class T>
bool check_node(SplicingGraph<T>& splicing_graph, 
                std::map<int, size_t>& distribution,
                std::vector<int>& node_path, int index) {

  int check_length = splicing_graph.node_set_[node_path[index]].sequence.length();

  for (int j = index-1; j >= 0; --j) {   //check the ancestors of node index
    if (distribution.find(node_path[j]) != distribution.end()) {
      return true;
    } else {
      check_length += splicing_graph.node_set_[node_path[j]].sequence.length();
      if (check_length > g_max_pair_gap_length)
        break;
    }
  }
  
  return (distribution.find(node_path[index]) != distribution.end());
}


template <class T>
std::pair<int, int> check_path(SplicingGraph<T>& splicing_graph,
	 			KmerMap<T>& kmer_map, 
			        std::vector<pair_t >& splicing_graph_edges,
				std::set<int>& check_nodes,
                                std::map<int,int>& edge_frequency,
                                std::map<pair_t, double> edge_cov_map,
 				std::vector<std::set<size_t> >& reads_mapped_to_node,
			  	std::vector<directed_acyclic_graph_node_t>& path,
				std::vector<std::vector<directed_acyclic_graph_node_t> >& paths_checked) {

  std::pair<int, int> flag(-1,-1); // return false and index if it has no support

  std::vector<int> node_path;
  for (size_t i = 0; i < path.size(); ++i) {
    pair_t node_pair = splicing_graph_edges[path[i]];
    if (i == 0)
      node_path.push_back(node_pair.first);
    node_path.push_back(node_pair.second);
  }

  /*
    relationship between chain and node_path here:
   
    chain:   0 --> 11 --> 10 --> 12   (node of compatibility graph)
           /  \  /  \   /  \   /  \
    path : 0 --> 4 --> 10 --> 7 --> 1 (node of splicing graph)
  */

  if (g_debug)
    std::cerr << "check path: " << describe_path(path) << std::endl;

  // for path containing only one node
  if (path.size() == 1) {
    /*
    int s = node_path[0];
    int t = node_path[1];
    int s_length = splicing_graph.node_set_[s].sequence.length();
    int t_length = splicing_graph.node_set_[t].sequence.length();
    if (s_length < 150 || t_length < 150) return flag;  // enough long to use pair
    int start = splicing_graph.node_set_[s].sequence.length() > 150 ?
        static_cast<int>(splicing_graph.node_set_[s].sequence.length())-150 : 0;
    const string& seq_s = splicing_graph.node_set_[s].sequence.substr(start, 100);
    const string& seq_t = splicing_graph.node_set_[t].sequence.substr(50, 100);
    std::set<size_t> reads_s;
    splicing_graph.set_reads(kmer_map, seq_s, reads_s);
    std::map<int,size_t> d_s;
    splicing_graph.get_mate_distribution(reads_s, d_s, reads_mapped_to_node, 1);
    std::set<size_t> reads_t;
    splicing_graph.set_reads(kmer_map, seq_t, reads_t);
    std::map<int,size_t> d_t;
    splicing_graph.get_mate_distribution(reads_t, d_t, reads_mapped_to_node, 0);
    bool del_edge = false;
    if (d_s.find(t) == d_s.end()) {
      if (t_length > g_max_pair_gap_length) 
        del_edge = true;
    } else {
      eds_span_junctionf (d_s[t] < g_min_reads_span_junction)
        del_edge = true;
    }
    if (!del_edge) {
      if (d_t.find(s) == d_t.end()) {
        if (s_length > g_max_pair_gap_length) 
          del_edge = true;
      } else {
        if (d_t[t] < g_min_reads_span_junction)
          del_edge = true;
      }
    }
    if (del_edge) {     // delete this edge
      if (g_debug)
        std::cerr << "delete edge: " << s << "->" << t << std::endl;
      splicing_graph.node_set_[s].delete_child(t);
      splicing_graph.node_set_[t].delete_parent(s);
      flag.second = -2; // indicates a isolated node in splicing graph
    }
    */
    return flag;
    
  }


  // for path containing multiple nodes
  std::vector<int> unique_edges;
  for (int i = 1; i < (int)path.size(); ++i) {
    pair_t node_pair = splicing_graph_edges[path[i]];
    int s = node_pair.first;
    int t = node_pair.second;
    if (edge_frequency[path[i]] == 1) {  
      if (splicing_graph.node_set_[s].children.size() >= 2 &&
	  splicing_graph.node_set_[t].parents.size() >= 2)   
	unique_edges.push_back(path[i]);
    }
    if (check_nodes.find(s) != check_nodes.end()) {
      if ((int)splicing_graph.node_set_[s].sequence.length() > g_max_pair_gap_length)
        continue;
      const string& edge = splicing_graph.get_edge_sequence(s,t);
      /* 
      std::string edge;
      if (splicing_graph.node_set_[s].sequence.length() > g_min_exon_length) 
        edge = splicing_graph.node_set_[s].sequence.substr(0, 100);
      else
        edge = splicing_graph.node_set_[s].sequence + splicing_graph.node_set_[t].sequence.substr(0,g_kmer_length-1);
      */
      
      std::set<size_t> reads;
      splicing_graph.set_reads(kmer_map, edge, reads);
      std::map<int, size_t> distribution;
      //splicing_graph.get_mate_distribution(reads, distribution, reads_mapped_to_node, 0);
      splicing_graph.get_mate_distribution(reads, distribution, reads_mapped_to_node, 2);  // fr-unstrand
      if (distribution.size() < 1) continue;    
      
      // find support in a region less than 500 nucleotides
      bool has_support = check_node(splicing_graph, distribution, node_path, i);
      if (!has_support) {
        flag.first = i-1;
        flag.second = i; 
        break;
      } 
    }
  }
  
  //cerr << "unique edges..." << endl;
  //  check deeply if there exist several unique edges in one path
  if (flag.first && unique_edges.size() >= 2) {  
    int e1 = unique_edges[0];
    for (int i = 1; i < (int)unique_edges.size(); ++i) {
      int e2 = unique_edges[i];
      int e1_t = splicing_graph_edges[e1].second;
      int e2_s = splicing_graph_edges[e2].first;
      if (splicing_graph.node_set_[e1_t].parents.size() <= 1 ||
          splicing_graph.node_set_[e2_s].children.size() <= 1)
        continue;
      // if their coverage is different
      double d1 = edge_cov_map[splicing_graph_edges[e1]];
      double d2 = edge_cov_map[splicing_graph_edges[e2]];
      if (d1/d2 > 5 || d2/d1 > 5) {
        for (int j = 0; j < (int)path.size(); ++j) {
	  if ((int)path[j] == e1) flag.second = j+1;
	  if ((int)path[j] == e2) flag.first = j-1;
        }
        if (flag.first >= 0 && flag.second >= 0) {
          cerr << rg_file << ": edge " << e1 << "(cov=" << d1 << ") and edge" <<  e2 << "(cov=" << d2 << ") are in one path." << endl;
          break;
        }
      }
      e1 = e2;
    }
  }

  return flag;

}

template <class T>
void check_paths(SplicingGraph<T>& splicing_graph,
                KmerMap<T>& kmer_map,
                std::vector<pair_t >& splicing_graph_edges,
                std::map<pair_t, double> edge_cov_map,
                std::vector<std::set<size_t> >& reads_mapped_to_node,
                std::vector<std::vector<directed_acyclic_graph_node_t> >& paths,
                std::vector<std::vector<directed_acyclic_graph_node_t> >& paths_checked) {

  time_t beg = time(NULL);
  std::vector<std::vector<directed_acyclic_graph_node_t> > new_paths;

  // check only for paired reads
  if (!g_is_paired_end) {
    paths_checked = paths;
    return;
  }
  // frequency of each edge
  std::map<int, int> edge_frequency;
  for (size_t i = 0; i < paths.size(); ++i) {
    for (size_t j = 0; j < paths[i].size(); ++j) {
      if (edge_frequency.find(paths[i][j]) == edge_frequency.end())
        edge_frequency[paths[i][j]] = 1;
      else 
        edge_frequency[paths[i][j]] += 1;
    }
  }

  std::map<pair_t, int> edge_index;
  for (int i = 0; i < (int)splicing_graph_edges.size(); ++i) {
    edge_index[splicing_graph_edges[i]] = i;
  }

  std::set<int> check_nodes;
  for (int i = 0; i < (int)splicing_graph.get_size(); ++i) {
    if (splicing_graph.node_set_[i].children.size() <= 1
       || splicing_graph.node_set_[i].parents.size() <= 1) {
      continue;
    } else if ((int)splicing_graph.node_set_[i].sequence.length() < g_max_pair_gap_length) {
      check_nodes.insert(i);
    }
  }

  for (unsigned int i = 0; i < paths.size(); ++i) {
    std::pair<int,int> flag = check_path(splicing_graph, kmer_map, splicing_graph_edges,
      check_nodes, edge_frequency, edge_cov_map, reads_mapped_to_node, paths[i], paths_checked);
    if (flag.first < 0) {
      if (flag.second == -2) continue;
      paths_checked.push_back(paths[i]);
    } else {
      int left = flag.first;
      int right = flag.second;
      //cerr << "left = " << left << ", right=" << right << endl;
      std::vector<directed_acyclic_graph_node_t> path_left(paths[i].begin(), paths[i].begin()+left+1);
      std::vector<directed_acyclic_graph_node_t> path_right(paths[i].begin()+right, paths[i].end());
      if (g_debug) {
         std::cerr << "this path is splited!" << std::endl;
         std::cerr << "path_left: " << describe_path(path_left) << std::endl;
         std::cerr << "path_right: " << describe_path(path_right) << std::endl;
      }
      // recover the right part based on left
      bool recover_left = recover_path(splicing_graph, kmer_map, splicing_graph_edges, edge_index,
	reads_mapped_to_node, path_left, paths[i][left+1], path_right, edge_frequency, 0);
      if (recover_left) {
        new_paths.push_back(path_left);
        //cerr << "new left path : " << describe_path(path_left) << endl;
      }
      // recover the left part based on right 
      bool recover_right = recover_path(splicing_graph, kmer_map, splicing_graph_edges, edge_index,
	reads_mapped_to_node, path_right, paths[i][right-1], path_left, edge_frequency, 1);
      if (recover_right) {
	new_paths.push_back(path_right);
	//cerr << "new right path : " << describe_path(path_right) << endl; 
      }
      if (!(recover_left || recover_right)) {
	paths_checked.push_back(paths[i]);
      }
    }
  }
 
  // add new paths to path_checked
  for (unsigned int i = 0; i < new_paths.size(); ++i) {
    // check redundancy
    bool no_redundancy = true;
    for (unsigned int j = 0; j < paths_checked.size(); ++j) {
      if (new_paths[i].size() == paths_checked[j].size()) {
	if (new_paths[i] == paths_checked[j]) {
          no_redundancy = false;
          break;
        }
      }
    }
    if (no_redundancy)
      paths_checked.push_back(new_paths[i]);
  }

  time_t end = time(NULL);
  if (g_debug)
    std::cerr << "check paths : " << (end-beg) << " s" << std::endl;
}


#endif


