#ifndef RANDOM_SAMPLER_HEADER
#define RANDOM_SAMPLER_HEADER

#include <random>
#include <vector>

#include "graph.hpp"

/** 
 * @brief Base weight added to all vertices to ensure a non-zero probability of selection. 
 * This prevents the sampler from completely ignoring vertices not yet found in a fort.
 */
inline constexpr double RS_EPSILON = 0.01;

/**
 * @class RandomSampler
 * @brief Manages weighted vertex selection based on fort frequency.
 * 
 * This class tracks how often each vertex appears in discovered forts and 
 * uses that data to bias the selection of vertices. 
 * Vertices that frequently appear in forts are weighted more heavily.
 */
class RandomSampler {
private:
  const Graph *graph;          //!< Pointer to the graph being sampled.
  std::vector<std::size_t> cv; //!< How many discovered forts contain each vertex.
  std::size_t nf;              //!< The total number of forts processed by the sampler.
  std::mt19937 gen;            //!< The Mersenne Twister engine for generating random samples.

  /**
   * @brief Calculates the selection weight for a specific vertex.
   * @param u The vertex to evaluate.
   * @return A double representing the relative probability of selecting @p u.
   */
  double get_weight(Vertex u) const;

  /**
   * @brief Computes the total sum of weights for all vertices not in the ignored set.
   * @param ignored A bitset of vertices to exclude from the sum.
   * @return The combined weight of all eligible vertices.
   */
  double sum_weights(const VertexSet &ignored) const;

  /** @brief Overload of sum_weights for VertexBitset input. */
  double sum_weights(const VertexBitset &ignored) const;

public:
  /**
   * @brief Constructs a sampler for a specific graph.
   * @param gi The graph to use for vertex indexing.
   */
  RandomSampler(const Graph *gi);

  /**
   * @brief Updates internal frequency counts (cv) using a new fort.
   * 
   * Incrementing these counts increases the likelihood that vertices 
   * within this fort will be sampled in future iterations.
   * @param fort The bitset representing a discovered fort.
   */
  void update_weights(const VertexSet &fort);

  /** @brief Overload of update_weights for VertexBitset input. */
  void update_weights(const VertexBitset &fort);

  /**
   * @brief Generates a random VertexSet using weighted sampling.
   * 
   * @param num_samples The number of vertices to attempt to include in the bitset.
   * @param ignored Vertices that should never be selected during this sampling.
   * @return A VertexSet populated with sampled vertices.
   */
  VertexSet sample_set(std::size_t num_samples, const VertexSet &ignored = {});

  /**
   * @brief Generates a random VertexBitset using weighted sampling.
   * 
   * @param num_samples The number of vertices to attempt to include in the bitset.
   * @param ignored Vertices that should never be selected during this sampling.
   * @return A VertexBitset populated with sampled vertices.
   */
  VertexBitset sample_bitset(std::size_t num_samples, VertexBitset ignored = {});
};

#endif