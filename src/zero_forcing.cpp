#include "zero_forcing.hpp"

#include <list>

bool is_valid_zero_forcing(const Graph &graph, const VertexSet &filled) {
  return zero_forcing_closure(graph, filled).size() == graph.get_order();
}

VertexSet zero_forcing_closure(const Graph &graph, VertexSet filled, std::size_t *pt) {
  std::size_t propagation_time;
	Vertex vert;

	for (propagation_time = 0; propagation_time < graph.get_order(); propagation_time++) {
		VertexSet active;

		for (Vertex u : filled) {
			std::size_t count = 0;
			for (Vertex v : graph.get_adjacent(u)) {
				if (filled.find(v) != filled.cend()) continue;
				count++;
				vert = v;
			}
			if (count == 1) active.insert(vert);
		}
		if (active.empty()) break;

		filled.insert(active.cbegin(), active.cend());
	}

  if (pt) {
	  if (filled.size() == graph.get_order()) *pt = propagation_time;
	  else *pt = std::numeric_limits<std::size_t>::max();
  }

  return filled;
}

std::size_t zero_forcing_wavefront(const Graph &graph) {
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

				S_new = std::move(zero_forcing_closure(graph,S_new));

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
;
}