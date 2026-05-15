#ifndef ZERO_FORCING_HEADER
#define ZERO_FORCING_HEADER

#include "graph.hpp"

/**
 * @brief Computes the zero forcing closure of a set of vertices.
 * 
 * This function iteratively applies the color-change rule: if a "colored" vertex 
 * has exactly one neighbor that is not colored, that neighbor becomes colored. 
 * This continues until no more forces are possible.
 * 
 * @param graph The graph to play the game on.
 * @param[in,out] filled The initial set of vertices; modified to include all vertices forced during the process.
 * @return The number of steps needed to propogate over the entire graph.
 */
std::size_t zero_forcing_closure(const Graph &graph, VertexSet &filled);

/** @brief Overload of zero_forcing_closure for VertexSet input. */
std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &filled);

/**
 * @brief Applies the wavefront algorithm to the graph
 * 
 * @param graph The graph to analyze.
 * @return The order of smallest found zero forcing set on @p graph
 */
std::size_t zero_forcing_wavefront(const Graph &graph);

#endif