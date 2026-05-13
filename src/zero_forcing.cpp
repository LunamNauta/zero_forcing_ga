#include "zero_forcing.hpp"

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

  typedef std::pair<VertexSet, std::size_t> ZFSolution;
  std::vector<ZFSolution> cl_pairs;
  cl_pairs.emplace_back(VertexSet(), 0);

  for (std::size_t rank_limit = 1; rank_limit <= graph.get_order(); rank_limit++) {
    for (std::vector<ZFSolution>::const_iterator it = cl_pairs.cbegin(); it != cl_pairs.cend(); it++) {
      auto [source, source_rank] = *it;

      for (Vertex u = 0; u < graph.get_order(); u++) {
        VertexSet candidate(source.cbegin(), source.cend());

        candidate.insert(u);
        VertexSet neighbors = graph.get_adjacent(u);
        candidate.insert(neighbors.cbegin(), neighbors.cend());

        zero_forcing_closure(graph, candidate);

        std::size_t candidate_rank = source_rank;
        if (source.find(u) == source.cend()) candidate_rank++;

        int neighbors_outside_source = 0;
        for (VertexSet::const_iterator it = neighbors.cbegin(); it != neighbors.cend(); it++) {
          if (source.find(*it) != source.cend()) continue;
          neighbors_outside_source++;
        }

        if (neighbors_outside_source > 1) candidate_rank += (neighbors_outside_source - 1);

        if (candidate_rank <= rank_limit) {
          bool already_present = false;

          for (std::vector<ZFSolution>::const_iterator it = cl_pairs.cbegin(); it != cl_pairs.cend(); it++) {
            if (it->first != candidate || it->second > candidate_rank) continue;
            already_present = true;
            break;
          }

          if (!already_present) {
            cl_pairs.emplace_back(candidate, candidate_rank);
            if (candidate.size() == graph.get_order()) return candidate_rank;
          }
        }
      }
    }
  }

  return graph.get_order();
}