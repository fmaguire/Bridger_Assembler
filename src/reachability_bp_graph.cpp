// This file was modified from file matching_merge.cpp in Cufflinks
// The original copyright info is listed below
// Copyright 2009 Cole Trapnell.
// Distributed under the  Distributed under Cufflinks LICENSE
// (See accompanying file LICENSE)

/*
 *  reachability_bp_graph.cpp 
 *
 */

#include "utility.h"
#include "reachability_bp_graph.h"


void create_reachability_bp_graph(directed_acyclic_graph_t& directed_acyclic_graph,
 reachability_bp_graph_t& reachability_bp_graph, const adjacency_list<>& transitive_closure) {

  std::vector<reachability_bp_graph_t::BNode> b_to_a;
  dag_to_bp_t dag_to_bp;
  //fprintf (stdout, "\tclosure edges:\t\t\t\%d\n", num_edges(transitive_closure));
  graph_traits < adjacency_list<> >::vertex_iterator v, vend;
  b_to_a.resize(num_vertices(transitive_closure));
  for (boost::tie(v, vend) = vertices(transitive_closure); v != vend; ++v) {
    dag_to_bp_t::iterator itr = dag_to_bp.find(*v);
    if (itr == dag_to_bp.end()) {
      reachability_bp_graph_t::ANode A = reachability_bp_graph.addANode();
      int a = reachability_bp_graph.aNodeId(A);
      reachability_bp_graph_t::BNode B = reachability_bp_graph.addBNode();
      int b = reachability_bp_graph.bNodeId(B);
      b_to_a[b] = A;
      dag_to_bp[*v] = make_pair(a, b);
    }
  }

  reachability_bp_graph.reserveEdge(num_edges(transitive_closure));
  reachability_bp_graph.reserveANode(num_vertices(transitive_closure));
  reachability_bp_graph.reserveBNode(num_vertices(transitive_closure));

  graph_traits < adjacency_list<> >::edge_iterator i, end;
  for (boost::tie(i, end) = edges(transitive_closure); i != end; ++i) {
    int a_id = -1;
    int b_id = -1;
    directed_acyclic_graph_node_t s = source(*i, transitive_closure);
    directed_acyclic_graph_node_t t = target(*i, transitive_closure);
    dag_to_bp_t::iterator itr = dag_to_bp.find(s);
    if (itr == dag_to_bp.end()) {
      assert (false);
    } else {
      a_id = itr->second.first;
    }
    itr = dag_to_bp.find(t);
    if (itr == dag_to_bp.end()) {
      assert(false);
    } else {
      b_id = itr->second.second;
    }

    if (in_degree(s, directed_acyclic_graph) == 0  /* virtual "source"? */
        || out_degree(t, directed_acyclic_graph) == 0 /* virtual "sink"?*/)
      continue;

    assert (a_id != -1);
    assert (b_id != -1);

    reachability_bp_graph_t::ANode a_node = reachability_bp_graph.nodeFromANodeId(a_id);
    reachability_bp_graph_t::BNode b_node = reachability_bp_graph.nodeFromBNodeId(b_id);
    reachability_bp_graph_t::ANode a_for_b = b_to_a[b_id];
    assert (a_for_b != a_node);

    reachability_bp_graph.addEdge(a_node, b_node);
  }
}


void add_weights_to_reachability_bp_graph(reachability_bp_graph_t& reachability_bp_graph,
 const PercentSplicedInMap& psi_map, reachability_bp_graph_t::UEdgeMap<long long>& weights) {

  for (reachability_bp_graph_t::UEdgeIt i(reachability_bp_graph); i!=lemon::INVALID; ++i) {

    reachability_bp_graph_t::ANode a = reachability_bp_graph.source(i);
    reachability_bp_graph_t::BNode b = reachability_bp_graph.target(i);
    
    directed_acyclic_graph_node_t a_dag = reachability_bp_graph.aNodeId(a);
    int a_id_for_b = reachability_bp_graph.aNodeId(b);
    reachability_bp_graph_t::ANode a_for_b = reachability_bp_graph.nodeFromANodeId(a_id_for_b);
    assert (a_for_b != lemon::INVALID);
    directed_acyclic_graph_node_t b_dag = a_id_for_b;
    double a_psi = psi_map[a_dag];
    double b_psi = psi_map[b_dag];

    double score = log (1.0 - abs(a_psi-b_psi));
    assert(score <= 0.0);
    if (score >= -1e-6)
      score = -1e-6;

    long long weight = (long long)(score * -1e6);
    if (weight < 0) {
    //if (weight < 0 || abs(1-a_psi) < 1e-6 || abs(1-b_psi) < 1e-6) {
      weight = 999999999;
    }

    weights[i] = weight;

    if (g_debug)
      progress(const_cast<char *>("[%d->%d] a_psi=%.3f, b_psi=%.3f, weight=%lld\n"), (int)a_dag, (int)b_dag, a_psi, b_psi, weight);
  }

}



