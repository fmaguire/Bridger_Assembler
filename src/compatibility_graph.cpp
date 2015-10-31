/*
 *  compatibility_graph.cpp
 *
 */


#include "compatibility_graph.h"
// header of tie()
#include <boost/tuple/tuple.hpp> 
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
//#include "transitive_reduction.h"

using namespace std;
using namespace boost;


// add terminal nodes to a dag
// copy from scaffold_graph.cpp in cufflinks software
// Copyright 2009 Cole Trapnell.
std::pair<directed_acyclic_graph_node_t, directed_acyclic_graph_node_t> 
 add_terminal_nodes(directed_acyclic_graph_t& bundle_dag) {

  std::vector<char> has_parent(num_vertices(bundle_dag)+2, false);
  std::vector<char> has_child (num_vertices(bundle_dag)+2, false);

  graph_traits < directed_acyclic_graph_t >::vertex_iterator u, uend;
  // tie() allows the assignment of two values of the pair to two separate variables
  // vertices() return a pair of type std::pair<vertex_iterator,vertex_iterator>.
  // The pair of iterators is assigned to the iterator u and uend.
  for (boost::tie(u, uend) = vertices(bundle_dag); u != uend; ++u) {
    graph_traits < directed_acyclic_graph_t >::adjacency_iterator v, vend;
    //adjacent_vertices() return a pair of type std::pair<adjacency_iterator,adjacency_iterator>
    for (tie(v,vend) = adjacent_vertices(*u, bundle_dag); v != vend; ++v) {
      directed_acyclic_graph_node_t U = *u;
      directed_acyclic_graph_node_t V = *v;
      has_parent[V] = true;
      has_child[U] = true;
    }
  }

  directed_acyclic_graph_node_t source = add_vertex(bundle_dag);
  directed_acyclic_graph_node_t sink = add_vertex(bundle_dag);

  //int num_attached_to_source = 0;
  //int num_attached_to_sink = 0;
  for (size_t i = 0; i < num_vertices(bundle_dag); ++i) {
    if (!has_parent[i] && i != sink && i != source) {
      //num_attached_to_source++;
      add_edge(source, i, bundle_dag);
    }
    if (!has_child[i] && i != source && i != sink) {
      //num_attached_to_sink++;
      add_edge(i, sink, bundle_dag);
    }
  }

  return make_pair(source, sink);
}




