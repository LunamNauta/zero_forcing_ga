#include "zero_forcing.hpp"
#include "graph.hpp"

#include <limits>
#include <list>

std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &filled) {
  std::vector<std::size_t> white_degree(graph.get_order());
  std::vector<Vertex> newly_filled;
  std::vector<Vertex> active;
  std::vector<Vertex> next;
  std::size_t pt = 0;

  for (Vertex u = 0; u < graph.get_order(); ++u) {
    for (Vertex v : graph.get_adjacent(u)) {
      if (!filled[v]) white_degree[u]++;
    }
    if (filled[u] && white_degree[u] == 1) active.push_back(u);
  }

  while (!active.empty()) {
    newly_filled.clear();
    next.clear();

    for (Vertex u : active) {
      if (white_degree[u] != 1) continue;

      for (Vertex v : graph.get_adjacent(u)) {
        if (filled[v]) continue;
        filled.set(v);
        newly_filled.push_back(v);
        break; 
      }
    }

    if (newly_filled.empty()) break;

    for (Vertex v : newly_filled) {
      for (Vertex neighbor : graph.get_adjacent(v)) {
        white_degree[neighbor]--;
        if (filled[neighbor] && white_degree[neighbor] == 1) next.push_back(neighbor);
      }
      if (white_degree[v] == 1) next.push_back(v);
    }

    active = std::move(next);
    pt++;
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