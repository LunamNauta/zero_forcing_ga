#ifndef ZERO_FORCING_HEADER
#define ZERO_FORCING_HEADER

#include "graph.hpp"

std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &filled);
std::size_t zero_forcing_closure(const Graph &graph, VertexSet &filled);

std::size_t zero_forcing_wavefront(const Graph &graph);

#endif