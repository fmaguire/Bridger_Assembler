// This file was modified from file assemble.cpp in Cufflinks
// The original copyright info is listed below
// Copyright 2009 Cole Trapnell.
// Distributed under the  Distributed under Cufflinks LICENSE
// (See accompanying file LICENSE)

#ifndef REACHABILITY_BP_GRAPH_H
#define REACHABILITY_BP_GRAPH_H

/*
 *  reachablity_bp_graph.h
 *
 */

#include <map>
#include <list>
#include <vector>
#include <iostream>

#define LEMON_ONLY_TEMPLATES

#include <lemon/topology.h>
#include <lemon/graph_utils.h>
#include <lemon/smart_graph.h>
#include <lemon/bipartite_matching.h>

#include <boost/ref.hpp>
// DON'T move this, or mystery compiler errors will result. Affects gcc >= 4.1
#include <boost/graph/vector_as_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/version.hpp>

#if (BOOST_VERSION < 103800)
#include <boost/vector_property_map.hpp>
#else
#include <boost/property_map/vector_property_map.hpp>
#endif

//#include <boost/math/distributions/normal.hpp> // for normal_distribution
//using boost::math::normal; // typedef provides default type is double.

#include "compatibility_graph.h"
#include "transitive_closure.h"

using namespace boost;
using namespace std;

typedef lemon::SmartBpUGraph reachability_bp_graph_t;
typedef map<directed_acyclic_graph_node_t, pair<int, int> > dag_to_bp_t;

void create_reachability_bp_graph(directed_acyclic_graph_t& directed_acyclic_graph,
  reachability_bp_graph_t& reachability_bp_graph, const adjacency_list<>& transitive_closure);

void add_weights_to_reachability_bp_graph(reachability_bp_graph_t& bp,
  const PercentSplicedInMap& psi_for_node, reachability_bp_graph_t::UEdgeMap<long long>& weights);

#endif

