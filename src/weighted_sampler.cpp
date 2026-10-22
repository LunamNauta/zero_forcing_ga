#include "weighted_sampler.hpp"

WeightedSampler::WeightedSampler(const Graph *graph) :
  graph(graph),
  wv(graph->order(), WS_ALPHA),
  generator(std::random_device{}())
{}

double WeightedSampler::sum_weights(const VertexBitset &ignored = {}) const {
  double sum = 0.0;
  for (Vertex u = 0; u < graph->order(); u++){
    if (ignored[u]) continue;
    sum += get_weight(u);
  }
  return sum;
}

double WeightedSampler::get_weight(Vertex u) const {
  return wv[u];
}

void WeightedSampler::update_weights(const std::unordered_set<VertexBitset> &forts) {
  for (Vertex u = 0; u < graph->order(); u++) {
    wv[u] *= WS_EPSILON;
    for (const VertexBitset &fort : forts) {
      if (!fort.test(u)) continue;
      wv[u] += 1.0 / fort.count();
    }
  }
}

VertexBitset WeightedSampler::sample(std::size_t num_samples, VertexBitset ignored, bool invert) {
  if (num_samples == 0) return VertexBitset(graph->order());
  ignored.resize(graph->order());

  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  double total_weight = sum_weights(ignored);
  std::vector<std::pair<Vertex, double>> candidates;
  candidates.reserve(graph->order());

  for (Vertex u = 0; u < graph->order(); ++u) {
    if (ignored[u]) continue;

    double mu = distribution(generator);
    double weight = get_weight(u);

    double exponent;
    if (!invert) exponent = total_weight / weight; // 1 / a_v
    else exponent = weight / total_weight;

    candidates.emplace_back(u, std::pow(mu, exponent));
  }

  std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  const std::size_t limit = std::min(num_samples, candidates.size());
  VertexBitset sample(graph->order());
  for (std::size_t a = 0; a < limit; a++) sample.set(candidates[a].first);

  return sample;
}
