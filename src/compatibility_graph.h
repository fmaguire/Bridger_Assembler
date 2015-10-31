#ifndef COMPATIBILITY_GRAPH_H
#define COMPATIBILITY_GRAPH_H

/*
 *  compatibility_graph.h
 *
 */

#include "utility.h"
#include "kmerhash.h"
#include "splicing_graph.h"
#include <vector>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/version.hpp>

#if (BOOST_VERSION < 103800)
#include <boost/vector_property_map.hpp>
#else
#include <boost/property_map/vector_property_map.hpp>
#endif


struct vertex_percent_spliced_in_t {
  typedef boost::vertex_property_tag kind;
};

typedef boost::adjacency_list<boost::vecS,
  boost::vecS,
  boost::bidirectionalS,
  boost::property<boost::vertex_name_t, std::string, 
  boost::property<vertex_percent_spliced_in_t, double> > > directed_acyclic_graph_t;

typedef boost::graph_traits<directed_acyclic_graph_t>::vertex_descriptor directed_acyclic_graph_node_t;

typedef boost::property_map<directed_acyclic_graph_t, vertex_percent_spliced_in_t>::type 
  PercentSplicedInMap;


typedef std::pair<int,int> pair_t;

template<class T>
void check_incompatible_edges(SplicingGraph<T> & splicing_graph,
                              KmerMap<T>& kmer_map,
                              std::vector<std::set<size_t> >& reads_mapped_to_node,
                              std::map<pair_t,int>& edge_index,
                              int node_id,
			      std::set<pair_t>& inhibit_edges) {

  assert(splicing_graph.node_set_[node_id].parents.size() >= 2);
  assert(splicing_graph.node_set_[node_id].children.size() >= 2);
  int p1 = splicing_graph.node_set_[node_id].parents[0];
  int p2 = splicing_graph.node_set_[node_id].parents[1];
  int c1 = splicing_graph.node_set_[node_id].children[0];
  int c2 = splicing_graph.node_set_[node_id].children[1];
  int p1_length = splicing_graph.node_set_[p1].sequence.length();
  int p2_length = splicing_graph.node_set_[p2].sequence.length();
  int c1_length = splicing_graph.node_set_[c1].sequence.length();
  int c2_length = splicing_graph.node_set_[c2].sequence.length();
  int length = splicing_graph.node_set_[node_id].sequence.length();
  int e1 = edge_index[pair_t(p1,node_id)];
  int e2 = edge_index[pair_t(p2,node_id)];
  int e3 = edge_index[pair_t(node_id,c1)];
  int e4 = edge_index[pair_t(node_id,c2)];

  /* draft 

    p1              c1
        \1       3/
          node_id 
        /2       4\
    p2              c2

  */

  bool p1_to_c1 = true;
  bool p1_to_c2 = true;
  bool p2_to_c1 = true;
  bool p2_to_c2 = true;
  const std::string& edge1 = splicing_graph.get_edge_sequence(p1, node_id);
  std::set<size_t> reads1;
  splicing_graph.set_reads(kmer_map, edge1, reads1);
  std::map<int,size_t> d1;
  //splicing_graph.get_mate_distribution(reads1,d1,reads_mapped_to_node, 1);
  splicing_graph.get_mate_distribution(reads1,d1,reads_mapped_to_node, 2);

  const std::string& edge2 = splicing_graph.get_edge_sequence(p2, node_id);
  std::set<size_t> reads2;
   splicing_graph.set_reads(kmer_map, edge2, reads2);
  std::map<int,size_t> d2;
  //splicing_graph.get_mate_distribution(reads2,d2,reads_mapped_to_node, 1);
  splicing_graph.get_mate_distribution(reads2,d2,reads_mapped_to_node, 2);

  const std::string& edge3 = splicing_graph.get_edge_sequence(node_id, c1);
  std::set<size_t> reads3;
  splicing_graph.set_reads(kmer_map, edge3, reads3);
  std::map<int,size_t> d3;
  //splicing_graph.get_mate_distribution(reads3,d3,reads_mapped_to_node, 0);
  splicing_graph.get_mate_distribution(reads3,d3,reads_mapped_to_node, 2);

  const std::string& edge4 = splicing_graph.get_edge_sequence(node_id, c2);
  std::set<size_t> reads4;
   splicing_graph.set_reads(kmer_map, edge4, reads4);
  std::map<int,size_t> d4;
  //splicing_graph.get_mate_distribution(reads4,d4,reads_mapped_to_node, 0);
  splicing_graph.get_mate_distribution(reads4,d4,reads_mapped_to_node, 2);


  int tolerance = 80;
  if ( ((!d1.empty() && d1.find(c1) == d1.end() && (length + c1_length >= g_pair_gap_length + tolerance)) ||
        (!d3.empty() && d3.find(p1) == d3.end() && (length + p1_length >= g_pair_gap_length + tolerance))) &&   // p1 -| c1
      // (d1.find(c2) != d1.end() || d4.find(p1) != d4.end()) &&     // p1 -> c2
       (d2.find(c1) != d2.end() || d3.find(p2) != d3.end()) ) {     // p2 -> c1
    if (d1.find(c2) != d1.end() || d4.find(p1) != d4.end()) {  // p1 -> c2
      p1_to_c1 = false;
    } else if (length + p1_length < g_pair_gap_length + tolerance && 
      length + c2_length < g_pair_gap_length + tolerance) {
       int parent = -1;
       int child = -1;
       if (splicing_graph.node_set_[p1].parents.size() == 1) 
	 parent = splicing_graph.node_set_[p1].parents[0];
       if (splicing_graph.node_set_[c2].children.size() == 1)
	 child = splicing_graph.node_set_[c2].children[0];
       if ( (parent > 0 && d4.find(parent) != d4.end()) ||
	    (child > 0 && d1.find(child) != d1.end()) )
	 p1_to_c1 = false; 
    }
  }
  if ( ((!d1.empty() && d1.find(c2) == d1.end() && (length + c2_length >= g_pair_gap_length + tolerance)) ||
        (!d4.empty() && d4.find(p1) == d4.end() && (length + p1_length >= g_pair_gap_length + tolerance))) &&   // p1 -| c2
       //(d1.find(c1) != d1.end() || d3.find(p1) != d3.end()) &&     // p1 -> c1
       (d2.find(c2) != d2.end() || d4.find(p2) != d4.end()) ) {      // p2 -> c2
    if (d1.find(c1) != d1.end() || d3.find(p1) != d3.end()) {        // p1 -> c1
      p1_to_c2 = false;
    } else if (length + c1_length <  g_pair_gap_length + tolerance &&
      length + p1_length < g_pair_gap_length + tolerance) {
      int parent = -1;
      int child = -1;
      if (splicing_graph.node_set_[p1].parents.size() == 1)
        parent = splicing_graph.node_set_[p1].parents[0];
      if (splicing_graph.node_set_[c1].children.size() == 1)
        child = splicing_graph.node_set_[c1].children[0];
      if ( (parent > 0 && d3.find(parent) != d3.end()) ||
           (child > 0 && d1.find(child) != d1.end()) )
        p1_to_c2 = false;
    }
  }


  if ( ((!d2.empty() && d2.find(c1) == d2.end() && (length + c1_length >= g_pair_gap_length + tolerance)) ||
        (!d3.empty() && d3.find(p2) == d3.end() && (length + p2_length >= g_pair_gap_length + tolerance))) &&   // p2 -| c1
       //(d2.find(c2) != d2.end() || d4.find(p2) != d4.end()) &&      // p2 -> c2
       (d1.find(c1) != d1.end() || d3.find(p1) != d3.end()) ) {      // p1 -> c1
    if (d2.find(c2) != d2.end() || d4.find(p2) != d4.end()) {        // p2 -> c2
      p2_to_c1 = false;
    } else if (length + c2_length <  g_pair_gap_length + tolerance &&
      length + p2_length < g_pair_gap_length + tolerance) {
      int parent = -1;
      int child = -1;
      if (splicing_graph.node_set_[p2].parents.size() == 1)
        parent = splicing_graph.node_set_[p2].parents[0];
      if (splicing_graph.node_set_[c2].children.size() == 1)
        child = splicing_graph.node_set_[c2].children[0];
      if ( (parent > 0 && d4.find(parent) != d4.end()) ||
           (child > 0 && d2.find(child) != d2.end()) )
         p2_to_c1 = false;
    }
  }

  if ( ((!d2.empty() && d2.find(c2) == d2.end() && (length + c2_length >= g_pair_gap_length + tolerance)) ||
        (!d4.empty() && d4.find(p2) == d4.end() && (length + p2_length >= g_pair_gap_length + tolerance))) &&   // p2 -| c2
      // (d2.find(c1) != d2.end() || d3.find(p2) != d3.end()) &&      // p2 -> c1
       (d1.find(c2) != d1.end() || d4.find(p1) != d4.end()) ) {      // p1 -> c2
    if (d2.find(c1) != d2.end() || d3.find(p2) != d3.end()) {        // p2 -> c1
      p2_to_c2 = false;
    } else if (length + c1_length <  g_pair_gap_length + tolerance &&
      length + p2_length < g_pair_gap_length + tolerance) {
      int parent = -1;
      int child = -1;
      if (splicing_graph.node_set_[p2].parents.size() == 1)
        parent = splicing_graph.node_set_[p2].parents[0];
      if (splicing_graph.node_set_[c1].children.size() == 1)
        child = splicing_graph.node_set_[c1].children[0];
      if ( (parent > 0 && d3.find(parent) != d3.end()) ||
           (child > 0 && d2.find(child) != d2.end()) )
         p2_to_c2 = false;
    }
  }

  if (!p1_to_c1) inhibit_edges.insert(pair_t(e1,e3));
  if (!p1_to_c2) inhibit_edges.insert(pair_t(e1,e4));
  if (!p2_to_c1) inhibit_edges.insert(pair_t(e2,e3));
  if (!p2_to_c2) inhibit_edges.insert(pair_t(e2,e4));

}

template<class T>
bool create_compatibility_graph(directed_acyclic_graph_t & bundle_dag, 
				SplicingGraph<T> & splicing_graph, 
				std::vector<pair_t> & splicing_graph_edges, 
				std::map<pair_t, double> & edge_cov_map,
				KmerMap<T>& kmer_map,
			        std::vector<std::set<size_t> >& reads_mapped_to_node) {

  if (g_debug)
    std::cerr << "Creating an compatibility graph..." << std::endl;

  time_t beg = time(NULL);
  if (!splicing_graph_edges.empty())
    splicing_graph_edges.clear();

  splicing_graph_edges.reserve(edge_cov_map.size());
  for (size_t i = 0; i < splicing_graph.node_set_.size(); ++i) {
    if (splicing_graph.node_set_[i].children.empty())
      continue;
    for (size_t j = 0; j < splicing_graph.node_set_[i].children.size(); ++j)
      splicing_graph_edges.push_back(pair_t(i,splicing_graph.node_set_[i].children[j]));
  }

  std::vector<double> coverage_of_edges(splicing_graph_edges.size());
  std::vector<double> total_coverage_from_sources(splicing_graph.get_size(), 0.001f);
  std::vector<double> total_coverage_into_targets(splicing_graph.get_size(), 0.001f);

  for (int i = 0; i < (int)splicing_graph_edges.size(); ++i) {
    int source = splicing_graph_edges[i].first;
    int target = splicing_graph_edges[i].second;
    coverage_of_edges[i] = edge_cov_map[splicing_graph_edges[i]];
    if (g_debug) 
      progress(const_cast<char *>("[%d] %d->%d\tcov=%.3f\n"),i, source, target, coverage_of_edges[i]);
    total_coverage_from_sources[source] += coverage_of_edges[i];
    total_coverage_into_targets[target] += coverage_of_edges[i];
  }

  // find some inhibit edges for compatibility graph
  std::set<pair_t> inhibit_edges;
  if (g_is_paired_end && (!reads_mapped_to_node.empty())) {
    std::map<pair_t, int> edge_index;
    for (int i = 0; i < (int)splicing_graph_edges.size(); ++i) {
      edge_index[splicing_graph_edges[i]] = i;
    }
    std::vector<int> check_nodes;
    for (int i = 0; i < (int)splicing_graph.get_size(); ++i) {
      if (splicing_graph.node_set_[i].children.size() > 1 &&
        splicing_graph.node_set_[i].parents.size() > 1) { 
        check_nodes.push_back(i);
      }
    }
    if (!check_nodes.empty()) {
      for (int i = 0; i < (int)check_nodes.size(); ++i) {
        int id = check_nodes[i];
        if ((int)splicing_graph.node_set_[id].sequence.length() < g_kmer_length ||
         (int)splicing_graph.node_set_[id].sequence.length() > g_max_pair_gap_length)
          continue; 
        check_incompatible_edges(splicing_graph, kmer_map, reads_mapped_to_node, edge_index, id, inhibit_edges);
      }
    }
  }
 
  if (g_debug) { 
    std::cerr << "inhibit edges..." << inhibit_edges.size() << std::endl;  
    std::set<pair_t>::iterator it;
    for (it = inhibit_edges.begin() ; it != inhibit_edges.end();++it) {
      std::cerr << (*it).first << " ->" << (*it).second<< std::endl;
    }
  }

  // build compatibility graph
  bundle_dag = directed_acyclic_graph_t(splicing_graph_edges.size());
  // add edges to dag
  for (int i = 0; i < (int)splicing_graph_edges.size(); ++i) {
    for (int j = 0; j < (int)splicing_graph_edges.size(); ++j) {
      if (inhibit_edges.find(pair_t(i,j)) != inhibit_edges.end()) 
        continue; // inhibit some edges
      if (splicing_graph_edges[i].second == splicing_graph_edges[j].first)
        boost::add_edge(i, j, bundle_dag);
    }
  }

  // compute percent-spliced-in metric for each node of dag
  PercentSplicedInMap
    psi_map = boost::get(vertex_percent_spliced_in_t(), bundle_dag);

  for (int i = 0; i < (int)splicing_graph_edges.size(); ++i) {

    int source = splicing_graph_edges[i].first;
    int target = splicing_graph_edges[i].second;

    //assert(total_coverage_from_sources[source] > 0);
    psi_map[i] = coverage_of_edges[i]/total_coverage_from_sources[source];
    if ( splicing_graph.node_set_[source].children.size() == 1 
      && psi_map[i] > coverage_of_edges[i]/total_coverage_into_targets[target]) {
        psi_map[i] = coverage_of_edges[i]/total_coverage_into_targets[target];
    }

    if (g_debug) 
      progress(const_cast<char *>("[%d] %d->%d\tpsi=%.3f\n"), i, source, target, psi_map[i]);
  }

  time_t end = time(NULL);

  if (g_debug)
    std::cerr << "compatibility graph: " << (end-beg) << " s" << std::endl;

  return true;
}

/*
template<class T>
bool create_compatibility_graph(directed_acyclic_graph_t& bundle_dag,
 SplicingGraph<T>& splicing_graph, std::vector<pair_t>& splicing_graph_edges) {

  //cerr << "Creating an compatibility graph..." << endl;                       
  splicing_graph_edges.clear();
  for (size_t i = 0; i < splicing_graph.node_set_.size(); ++i) {
    if (splicing_graph.node_set_[i].children.empty())
      continue;
    for (size_t j = 0; j < splicing_graph.node_set_[i].children.size(); ++j)
      splicing_graph_edges.push_back(pair_t(i,splicing_graph.node_set_[i].children[j]));
  }

  std::vector<double> coverage_of_edges(splicing_graph_edges.size());
  std::vector<double> total_coverage_from_sources(splicing_graph.get_size(), 0.001);
  std::vector<double> total_coverage_into_targets(splicing_graph.get_size(), 0.001);

  for (size_t i = 0; i < splicing_graph_edges.size(); ++i) {
 
    int source = splicing_graph_edges[i].first;
    int target = splicing_graph_edges[i].second;
    // compute the edge from source to target
    int start = splicing_graph.node_set_[source].sequence.size() > (g_kmer_length-1) ?
       (splicing_graph.node_set_[source].sequence.size()-g_kmer_length+1) : 0;
    std::string edge = splicing_graph.node_set_[source].sequence.substr(start); //left part
    int length = g_kmer_length-1;
    if (splicing_graph.node_set_[target].sequence.length() < length) { 
      length = splicing_graph.node_set_[target].sequence.length();
      if (splicing_graph.node_set_[target].children.size() == 1) {
	int child = splicing_graph.node_set_[target].children[0];
        int remain = g_kmer_length-1-length;
	edge = edge + splicing_graph.node_set_[target].sequence
	 + splicing_graph.node_set_[child].sequence.substr(0, remain);
      } else {
	edge += splicing_graph.node_set_[target].sequence.substr(0, length);
      }
    } else {
      edge += splicing_graph.node_set_[target].sequence.substr(0, length);
    }
   
    //cerr << "edge=" << edge << endl;
    int sum = 0; // sum of all kmer coverage 
    int size = static_cast<int>(edge.size())-g_kmer_length+1;
   
    if (size > 0) {
      for (int j = 0; j < size; ++j) {
        std::string kmer = edge.substr(j, g_kmer_length);
        kmer_int_type_t intval = kmer_to_intval(kmer);
	 
	//if (used_kmers.find(intval) == used_kmers.end()) {
        if (! splicing_graph.is_used(intval)) {
	  std::cerr << kmer << "("<< intval <<")\tedge=" << edge << std::endl;
          //assert (used_kmers.find(intval) != used_kmers.end());
	  if ( j >= 19) {
	    size = j;
	    break;
	  } else {
	    std::cerr << "Warning: this kmer does not exists" << std::endl;
	    // delete this edge
	    //splicing_graph.node_set_[source].delete_child(target);
	    //splicing_graph.node_set_[target].delete_parent(source);
	    //splicing_graph_edges.erase(splicing_graph_edges.begin()+i); // potential error(--i)
	    std::cerr <<"should delete edge(" << source << "," << target<< ")" << std::endl;
	    break;
	  }
	}  
        sum += splicing_graph.get_kmer_count(intval); 
      }
    } else {
      err("node %d and %d are too small!\n", source, target);
    }
    if (size > 0) {
      coverage_of_edges[i] = (double)sum*1.0/size;
    }
    if (g_debug) 
      progress(const_cast<char *>("[%d] %d->%d\tcov=%.3f\n"),i, source, target, coverage_of_edges[i]);
    
    total_coverage_from_sources[source] += coverage_of_edges[i];
    total_coverage_into_targets[target] += coverage_of_edges[i];
  }

  // build compatibility graph
  bundle_dag = directed_acyclic_graph_t(splicing_graph_edges.size());
  // add edges to dag
  for (size_t i = 0; i < splicing_graph_edges.size(); ++i) {
    for (size_t j = 0; j < splicing_graph_edges.size(); ++j) {
      //if ( (i == 13 && j == 11 )||( i ==1 && j ==15) || (i==8 &&j ==6)) continue; // test gene1
      if (splicing_graph_edges[i].second == splicing_graph_edges[j].first)
        boost::add_edge(i, j, bundle_dag);
    }
  }

  // compute percent-spliced-in metric for each node of dag
  PercentSplicedInMap
    psi_map = boost::get(vertex_percent_spliced_in_t(), bundle_dag);

  for (size_t i = 0; i < splicing_graph_edges.size(); ++i) {
    int source = splicing_graph_edges[i].first;
    int target = splicing_graph_edges[i].second;
    //assert(total_coverage_from_sources[source] > 0);
    psi_map[i] = coverage_of_edges[i]/total_coverage_from_sources[source];
    if ( splicing_graph.node_set_[source].children.size() == 1 
      && psi_map[i] > coverage_of_edges[i]/total_coverage_into_targets[target]) {
        psi_map[i] = coverage_of_edges[i]/total_coverage_into_targets[target];
    }

    if (g_debug) 
      progress(const_cast<char *>("[%d] %d->%d\tpsi=%.3f\n"), i, source, target, psi_map[i]);

  }
  
  return true;
}
*/

std::pair<directed_acyclic_graph_node_t, directed_acyclic_graph_node_t> 
 add_terminal_nodes(directed_acyclic_graph_t& bundle_dag);

#endif

