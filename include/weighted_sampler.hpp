#ifndef WEIGHTED_SAMPLER_HEADER
#define WEIGHTED_SAMPLER_HEADER

#include <random>
#include <vector>

#include "graph.hpp"

inline constexpr double WS_EPSILON = 0.9;
inline constexpr double WS_ALPHA = 0.01;

class WeightedSampler {
private:
  const Graph *graph;
  std::vector<double> wv;
  std::mt19937 generator;

  double sum_weights(const VertexBitset &ignored) const;

public:
  WeightedSampler(const Graph *graph);

  double get_weight(Vertex u) const;

  void update_weights(const std::unordered_set<VertexBitset> &forts);

  VertexBitset sample(std::size_t num_samples, VertexBitset ignored = {}, bool invert = false);
};

#endif
