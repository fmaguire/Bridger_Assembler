/*
 *  matching_merge.cpp
 *
 */

#include "matching_merge.h"

using namespace std;
using namespace boost;


// check if there exists a path from p to q using BFS
bool has_path(int p, int q, std::vector<pair_t >& splicing_graph_edges) {

  vector<int> descendant;
  descendant.clear();
  for (size_t i = 0; i< splicing_graph_edges.size(); ++i) {
    pair_t edge = splicing_graph_edges[i];
    if (edge.first == p) {
      if (edge.second == q) 
	return true;
      else
	descendant.push_back(edge.second);
    }
  }

  while (!descendant.empty()) {
    p = descendant.front();
    descendant.erase(descendant.begin());
    for (size_t i = 0; i< splicing_graph_edges.size(); ++i) {
      pair_t edge = splicing_graph_edges[i];
      if (edge.first == p) {
        if (edge.second == q )
          return true;
        else
          descendant.push_back(edge.second);
      }
    }
  }

  return false;
}

/*
template class<T>
void find_path(const directed_acyclic_graph_t& bundle_dag,
               const adjacency_list<>& TC,
               const directed_acyclic_graph_node_t& source,
               const directed_acyclic_graph_node_t& target,
               std::vector<directed_acyclic_graph_node_t>& path,
	       std::vector<pair_t >& splicing_graph_edges,
               SplicingGraph<T>& splicing_graph,
	       std::map<pair_t,double> & edge_cov_map) {

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
        // add by myself from here 
	if ( (i+1) != iend ) {
	  directed_acyclic_graph_node_t J = *(i+1);
	  pair<adjacency_list<>::edge_descriptor, bool> q;
	  q = edge(J, target, TC);
	  if (q.second) {  
            // target node of two edge
	    int tp = splicing_graph_edges[*i].second;
	    int tq = splicing_graph_edges[*(i+1)].second;
	    // choose a longer path
	    if (splicing_graph.has_path(tq, tp)) { 
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
	} 
        path.push_back(*i);
        curr = *i;
        break;
      }

      if (*i == target) {
        path.push_back(*i);
        done = true;
        break;
      }
    }
  }
}
*/


/*
template class<T>
void extend_chains_to_paths(const directed_acyclic_graph_t& bundle_dag,
                            std::vector<std::vector<directed_acyclic_graph_node_t> >& chains,
                            adjacency_list<>& transitive_closure,
                            directed_acyclic_graph_node_t source_node,
                            directed_acyclic_graph_node_t sink_node,
                            std::vector<std::vector<directed_acyclic_graph_node_t> >& paths,
			    std::vector<pair_t >& splicing_graph_edges,
                            SplicingGraph<T>& splicing_graph,
			    std::map<pair_t, double>& edge_cov_map) {

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
}
*/


// output sequence in multiple line, each line has the same length
std::string add_fasta_seq_line_breaks(const std::string& sequence, int interval) {

  std::stringstream fasta_seq;
  int counter = 0;
  for (std::string::const_iterator it = sequence.begin(); 
   it != sequence.end(); ++it) {
    counter++;
    fasta_seq << *it;
    if (counter % interval == 0 && (it+1) != sequence.end()) 
      fasta_seq << std::endl; 
  }
  return fasta_seq.str();
}


// describe the information of a path 
std::string describe_path(const std::vector<directed_acyclic_graph_node_t>& path) {

  std::stringstream str_path;
  for (size_t i = 0; i < path.size(); ++i) {
    str_path << path[i] ;
    if (i != path.size() -1 ) 
      str_path << "->";
  }
  return str_path.str();
}




