#ifndef ZERO_FORCING_HEADER
#define ZERO_FORCING_HEADER

#include "graph.hpp"

std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &filled);
std::size_t zero_forcing_closure(const Graph &graph, VertexSet &filled);

std::size_t zero_forcing_wavefront(const Graph &graph, std::size_t upper_bound = std::numeric_limits<std::size_t>::max());
std::size_t zero_forcing_wavefront_old(const Graph &graph);

VertexBitset minimum_fort_ip_subgraph(const Graph &graph, const VertexBitset &filled);

#endif