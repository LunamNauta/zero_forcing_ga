#ifndef RANDOM_SAMPLER_HEADER
#define RANDOM_SAMPLER_HEADER

#include <random>

#include "graph.hpp"

class RandomSampler {
private:
  const Graph *graph;
  std::vector<double> weights;
  double sum_weights;
  std::mt19937 gen;

public:
  RandomSampler(const Graph *gi);

  void update_weight(VertexIndex vert, std::size_t delta);
  void update_weights(const std::vector<std::pair<VertexIndex, std::size_t>> &vert_deltas);

  double get_weight(VertexIndex vert) const;

  Graph operator()(std::size_t num_samples, std::size_t max_attempts);
};

#endif