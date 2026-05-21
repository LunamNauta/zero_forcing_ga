#ifndef RANDOM_SAMPLER_HEADER
#define RANDOM_SAMPLER_HEADER

#include <random>
#include <vector>

#include "graph.hpp"

inline constexpr double RS_EPSILON = 0.01;

class RandomSampler {
private:
  const Graph *graph;
  std::vector<std::size_t> cv;
  std::size_t nf;
  std::mt19937 gen;

  double get_weight(Vertex u) const;

  double sum_weights(const VertexSet &ignored) const;
  double sum_weights(const VertexBitset &ignored) const;

public:
  RandomSampler(const Graph *gi);

  void update_weights(const VertexSet &fort);
  void update_weights(const VertexBitset &fort);

  VertexBitset sample_bitset(std::size_t num_samples, VertexBitset ignored = {}, bool invert = false);
  VertexSet sample_set(std::size_t num_samples, const VertexSet &ignored = {}, bool invert = false);
};

#endif