#include "random_sampler.hpp"

#include <algorithm>
#include <cassert>

#include "graph.hpp"

RandomSampler::RandomSampler(const Graph *gi) :
  graph(gi),
  num_fort(gi->get_order(), 0),
  fort_count(0),
  gen(std::random_device{}())
{}

void RandomSampler::update_weights(const VertexSet &fort) {
  for (Vertex u : fort) num_fort[u]++;
  fort_count++;
}

double RandomSampler::get_weight(Vertex vert) const {
  return RS_EPSILON + static_cast<double>(num_fort[vert]) / fort_count;
}

double RandomSampler::sum_weights(const VertexSet &ignored) const {
  double sum = 0;
  for (Vertex u = 0; u < graph->get_order(); u++){
    if (ignored.find(u) != ignored.cend()) continue;
    sum += get_weight(u);
  }
  return sum;
}

VertexSet RandomSampler::operator()(std::size_t num_samples, std::size_t max_attempts, const VertexSet &ignored) {
  if (num_samples == 0 || graph->get_order() == 0) return {};
  num_samples = std::min(num_samples, graph->get_order());

  // Find the sum of all the weights (for normalization)
  double total_weight = sum_weights(ignored);

  struct WeightedVertex {
    double weight;
    double tie_breaker;
    Vertex u;
  };

  std::uniform_real_distribution<double> dist1(0.0, 1.0);

  std::vector<WeightedVertex> weighted_vertices(graph->get_order());
  for (std::size_t u = 0; u < graph->get_order(); u++) {
    if (ignored.find(u) == ignored.cend()) { // If we are not ignoring this vertex
      double u_rand = dist1(gen); 
      double weight = get_weight(u);
      // Weight for each vertex is u^(cv/nf)
      weighted_vertices[u] = { std::pow(u_rand, 1.0 / weight), dist1(gen), u };
    }
    else weighted_vertices[u] = {-1.0, 0, u}; // Ignored vertices should never be sampled
  }

  // Sort based on weight
  std::sort(weighted_vertices.begin(), weighted_vertices.end(), [](const WeightedVertex& a, const WeightedVertex& b) {
    if (std::abs(a.weight - b.weight) > 1e-9) return a.weight > b.weight;
    return a.tie_breaker > b.tie_breaker;
  });

  // Select the top 'num_samples' weighted vertices
  VertexSet samples;
  for (std::size_t a = 0; a < num_samples; a++) {
    samples.insert(weighted_vertices[a].u);
  }

  return samples;
}