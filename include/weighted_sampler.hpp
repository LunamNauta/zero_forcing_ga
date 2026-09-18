#ifndef WEIGHTED_SAMPLER_HEADER
#define WEIGHTED_SAMPLER_HEADER

#include <random>
#include <vector>

#include "graph.hpp"

inline constexpr double RS_EPSILON = 0.01;
inline constexpr double RS_ALPHA = 0.9;

class WeightedSampler {
private:
  const Graph *graph;
  std::vector<double> cv;
  std::size_t nf;
  std::mt19937 generator;

  double sum_weights(const VertexBitset &ignored) const;
  double sum_weights(const VertexSet &ignored) const;

public:
  WeightedSampler(const Graph *graph);

  double get_weight(Vertex u) const;

  void reset_weights();

  void update_weights(const VertexBitset &fort);
  void update_weights(const VertexSet &fort);

  VertexBitset sample_bitset(std::size_t num_samples, const VertexBitset &ignored = {}, bool invert = false);
  VertexSet sample(std::size_t num_samples, const VertexSet &ignored = {}, bool invert = false);
};

#endif
