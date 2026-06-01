#include "zero_forcing.hpp"

#include <limits>
#include <list>

#include "graph.hpp"

std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &filled) {
  VertexBitset next(graph.get_order());
  std::size_t pt;

  for (pt = 0; pt < graph.get_order(); pt++){
    next.reset();

    for (Vertex u = 0; u < graph.get_order(); u++) {
      if (!filled.test(u)) continue;
      std::size_t white_count = 0;
      Vertex forced_vertex;

      for (Vertex v : graph.get_adjacent(u)) {
        if (filled.test(v)) continue;
        if (++white_count > 1) break;
        forced_vertex = v;
      }
      if (white_count == 1) next.set(forced_vertex);
    }
    if (next.count() == 0) break;
    
    filled |= next;
  }

  if (filled.count() != graph.get_order()) return std::numeric_limits<std::size_t>::max();
  return pt;
}

std::size_t zero_forcing_closure(const Graph &graph, VertexSet &filled) {
  VertexSet next;
  std::size_t pt;

  for (pt = 0; pt < graph.get_order(); pt++){
    next.clear();

    for (Vertex u : filled) {
      std::size_t white_count = 0;
      Vertex forced_vertex;

      for (Vertex v : graph.get_adjacent(u)) {
        if (filled.find(v) != filled.end()) continue;
        if (++white_count > 1) break;
        forced_vertex = v;
      }
      if (white_count == 1) next.insert(forced_vertex);
    }
    if (next.empty()) break;
    
    filled.insert(next.begin(), next.end());
  }

  if (filled.size() != graph.get_order()) return std::numeric_limits<std::size_t>::max();
  return pt;
}

std::size_t zero_forcing_wavefront(const Graph &graph, std::size_t upper_bound) {
  if (graph.get_order() == 0) return 0;
  upper_bound = std::min(graph.get_order(), upper_bound);

  std::list<std::pair<VertexBitset, std::size_t>> cl_pairs;
  cl_pairs.emplace_back(VertexBitset(graph.get_order(), false), 0);

  for (std::size_t R = 1; R <= upper_bound - 1; R++) {
    std::list<std::pair<VertexBitset, std::size_t>> next_cl_pairs;

    // Use structured binding for cleaner access
    for (const auto& [S, r] : cl_pairs) {
      for (Vertex v = 0; v < graph.get_order(); v++) {
        std::size_t r_new = r;
        if (!S.test(v)) r_new++;

        const VertexSet& neighbors = graph.get_adjacent(v);
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

        if (S_new.count() == graph.get_order()) return r_new;

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

std::size_t zero_forcing_wavefront_old(const Graph &graph) {
  if (graph.get_order() == 0) return 0;

	std::list<std::pair<VertexSet, std::size_t>> cl_pairs;
	cl_pairs.push_back(std::pair<VertexSet, std::size_t>({}, 0));

	for (std::size_t R = 1; R <= graph.get_order(); R++) {
		for (auto j = cl_pairs.cbegin(); j != cl_pairs.cend(); j++) {
			const VertexSet &S = j->first;
			std::size_t r = j->second;

			for (Vertex v = 0; v < graph.get_order(); v++) {
				VertexSet S_new(S.cbegin(), S.cend());

				S_new.insert(v);

				VertexSet neighbors = graph.get_adjacent(v);
				S_new.insert(neighbors.cbegin(), neighbors.cend());

				zero_forcing_closure(graph,S_new);

				std::size_t r_new = r;
				if (S.find(v) == S.cend()) r_new++;

				std::size_t neighbors_outside_S = 0;
				for (VertexSet::const_iterator it = neighbors.cbegin(); it != neighbors.cend(); it++) {
					if (S.find(*it) != S.cend()) continue;

					neighbors_outside_S++;
				}

				if (neighbors_outside_S > 1) r_new += (neighbors_outside_S - 1);

				if (r_new <= R) {
					bool already_present = false;

					for (auto it = cl_pairs.cbegin(); it != cl_pairs.cend(); it++) {
						if (it->first != S_new || it->second > r_new) continue;

						already_present = true;
						break;
					}

					if (!already_present) {
						cl_pairs.push_back(std::pair<VertexSet, int>(S_new, r_new));

						if (S_new.size() == graph.get_order()) return r_new;
					}
				}
			}
		}
	}

	return graph.get_order();
}