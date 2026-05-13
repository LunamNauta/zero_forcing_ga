#include "random_sampler.hpp"

#include <algorithm>
#include <cassert>

#include "graph.hpp"

// epsilon + nv / nf per weights
// normalize weight
// take 1 / normalized
// take random (uniform value)
// raise it to the power of 1 / normalized

RandomSampler::RandomSampler(const Graph *gi) :
  graph(gi),
  num_fort(gi->get_order(), 0),
  fort_count(0),
  gen(std::random_device{}())
{}

void RandomSampler::update_weights(const VertexSet &fort) {
  for (Vertex vert : fort) {
    num_fort[vert]++;
  }
  fort_count++;
}

double RandomSampler::get_weight(Vertex vert) const {
  return RS_EPSILON + static_cast<double>(num_fort[vert]) / fort_count;
}

double RandomSampler::sum_weights() const {
  double sum = 0;
  for (std::size_t a = 0; a < graph->get_order(); a++) {
    sum += get_weight(a);
  }
  return sum;
}

VertexSet RandomSampler::operator()(std::size_t num_samples, std::size_t max_attempts, const VertexSet &ignored) {
  if (num_samples == 0 || graph->get_order() == 0) return {};
  num_samples = std::min(num_samples, graph->get_order());

  std::uniform_real_distribution<double> dist1(0.0, 1.0);
  double base = dist1(gen);
  double total_weight = sum_weights();

  struct WeightedVertex {
    double weight;
    double tie_breaker;
    std::size_t index;
  };

  std::vector<WeightedVertex> weighted_vertices(graph->get_order());
  for (std::size_t a = 0; a < graph->get_order(); a++) {
    if (ignored.find(a) == ignored.cend()) {
      double factor = 1.0 / (get_weight(a) / total_weight);
      weighted_vertices[a] = {std::pow(base, factor), dist1(gen), a};
    }
    else weighted_vertices[a] = {0, 0, a};
  }

  std::sort(weighted_vertices.begin(), weighted_vertices.end(), 
            [](const WeightedVertex& a, const WeightedVertex& b) {
              if (std::abs(a.weight - b.weight) > 1e-9) {
                return a.weight > b.weight;
              }
              return a.tie_breaker > b.tie_breaker;
            });

  VertexSet samples;
  std::size_t limit = std::min(num_samples, weighted_vertices.size());
  for (std::size_t i = 0; i < limit; i++) {
    samples.insert(static_cast<Vertex>(weighted_vertices[i].index));
  }

  return samples;
}