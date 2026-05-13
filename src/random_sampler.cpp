#include "random_sampler.hpp"

#include <algorithm>
#include <cassert>

#include "graph.hpp"

RandomSampler::RandomSampler(const Graph *gi) :
  graph(gi), weights(gi->get_order(), 1.0f), sum_weights(gi->get_order()), gen(std::random_device{}())
{}

void RandomSampler::update_weight(VertexIndex vert, std::size_t delta) {
  assert(vert < weights.size());
  weights[vert] += delta;
  sum_weights += delta; 
}

void RandomSampler::update_weights(const std::vector<std::pair<VertexIndex, std::size_t>> &vert_deltas) {
  for (const std::pair<VertexIndex, std::size_t> &vert_delta : vert_deltas) {
    update_weight(vert_delta.first, vert_delta.second);
  }
}

double RandomSampler::get_weight(VertexIndex vert) const {
  return (sum_weights > 0) ? (weights[vert] / sum_weights) : 0;
}

Graph RandomSampler::operator()(std::size_t num_samples, std::size_t max_attempts) {
  if (num_samples <= 0 || graph->get_order() == 0) return Graph();
  num_samples = std::min(num_samples, graph->get_order());

  std::discrete_distribution<int> dist(weights.begin(), weights.end());

  std::unordered_set<VertexIndex> sampled_indices;
  std::size_t attempts = 0;
  while (sampled_indices.size() < num_samples && attempts < max_attempts) {
    sampled_indices.insert(dist(gen));
    attempts++;
  }

  return graph->subgraph(sampled_indices);
}