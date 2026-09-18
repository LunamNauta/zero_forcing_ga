#include "../include/zero_forcing.hpp"

#include <list>

std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &closure) {
  std::vector<std::size_t> uncolored_degree(graph.order(), 0);

  std::vector<Vertex> current_ready; 
  std::vector<Vertex> next_ready;
  current_ready.reserve(graph.order());
  next_ready.reserve(graph.order());

  for (Vertex u = 0; u < graph.order(); u++) {
    for (Vertex adjacent: graph.adjacent(u)) uncolored_degree[u] += !closure.test(adjacent);
    if (closure.test(u) && uncolored_degree[u] == 1) current_ready.push_back(u);
  }

  std::size_t pt = 0;
  while (!current_ready.empty()) {
    while (!current_ready.empty()) {
      Vertex u = current_ready.back();
      current_ready.pop_back();
      if (uncolored_degree[u] != 1) continue;

      Vertex v = graph.order(); 
      for (Vertex vi : graph.adjacent(u)) {
        if (closure.test(vi)) continue;
        v = vi;
        break;
      }
      if (v == graph.order()) continue; 

      closure.set(v);
      for (Vertex w : graph.adjacent(v)) {
        uncolored_degree[w]--;      
        if (closure.test(w) && uncolored_degree[w] == 1) next_ready.push_back(w);
      }
          
      if (uncolored_degree[v] == 1) next_ready.push_back(v);
    }

    pt++;
    
    current_ready = std::move(next_ready);
    next_ready.clear();
  }

  if (closure.count() != graph.order()) return std::numeric_limits<std::size_t>::max();
  return pt;
}


std::size_t zero_forcing_closure(const Graph &graph, VertexSet &filled) {
  VertexSet next;
  std::size_t pt;

  for (pt = 0; pt < graph.order(); pt++){
    next.clear();

    for (Vertex u : filled) {
      std::size_t white_count = 0;
      Vertex forced_vertex;

      for (Vertex v : graph.adjacent(u)) {
        if (filled.find(v) != filled.end()) continue;
        if (++white_count > 1) break;
        forced_vertex = v;
      }
      if (white_count == 1) next.insert(forced_vertex);
    }
    if (next.empty()) break;
    
    filled.insert(next.begin(), next.end());
  }

  if (filled.size() != graph.order()) return std::numeric_limits<std::size_t>::max();
  return pt;
}

std::size_t zero_forcing_wavefront(const Graph &graph, std::size_t upper_bound) {
  upper_bound = std::min(graph.order(), upper_bound);

  std::list<std::pair<VertexBitset, std::size_t>> cl_pairs;
  cl_pairs.emplace_back(VertexBitset(graph.order(), false), 0);

  for (std::size_t R = 1; R <= upper_bound - 1; R++) {
    std::list<std::pair<VertexBitset, std::size_t>> next_cl_pairs;

    // Use structured binding for cleaner access
    for (const auto& [S, r] : cl_pairs) {
      for (Vertex v = 0; v < graph.order(); v++) {
        std::size_t r_new = r;
        if (!S.test(v)) r_new++;

        const VertexSet& neighbors = graph.adjacent(v);
        std::size_t neighbors_outside_S = 0;
        
        for (Vertex neighbor : neighbors) {
          if (!S.test(neighbor)) neighbors_outside_S++;
        }
        
        if (neighbors_outside_S > 1) r_new += (neighbors_outside_S - 1);

        // Prune paths early
        if (r_new != R) continue;

        VertexBitset S_new = S; 
        S_new.set(v);
        for (Vertex u : neighbors) S_new.set(u);

        zero_forcing_closure(graph, S_new);

        if (S_new.count() == graph.order()) return r_new;

        bool already_present = false;
        
        // Check the main list for duplicates or better paths
        for (const auto& existing : cl_pairs) {
          if (existing.first == S_new && existing.second <= r_new) {
            already_present = true;
            break;
          }
        }

        // Check the pending generation list to prevent duplicate branching
        if (!already_present) {
          for (const auto& pending : next_cl_pairs) {
            if (pending.first == S_new && pending.second <= r_new) {
              already_present = true;
              break;
            }
          }
        }

        if (!already_present) {
          next_cl_pairs.emplace_back(std::move(S_new), r_new);
        }
      }
    }
    cl_pairs.splice(cl_pairs.end(), next_cl_pairs);
  }

  return upper_bound;
}
