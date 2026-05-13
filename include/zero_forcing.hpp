#ifndef ZERO_FORCING_HEADER
#define ZERO_FORCING_HEADER

#include "graph.hpp"

bool is_valid_zero_forcing(const Graph &graph, const VertexSet &filled);
VertexSet zero_forcing_closure(const Graph &graph, VertexSet filled, std::size_t *pt = nullptr);
std::size_t zero_forcing_wavefront(const Graph &graph);

#endif