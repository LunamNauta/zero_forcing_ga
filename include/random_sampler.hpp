#ifndef RANDOM_SAMPLER_HEADER
#define RANDOM_SAMPLER_HEADER

#include <random>

#include "graph.hpp"

inline constexpr double RS_EPSILON = 0.01;

class RandomSampler {
private:
  const Graph *graph;

  std::vector<std::size_t> num_fort;
  std::size_t fort_count;

  std::mt19937 gen;

public:
  RandomSampler(const Graph *gi);

  void update_weights(const VertexSet &fort);

  double get_weight(Vertex vert) const;
  double sum_weights() const;

  VertexSet operator()(std::size_t num_samples, std::size_t max_attempts, const VertexSet &ignore = {});
};

#endif