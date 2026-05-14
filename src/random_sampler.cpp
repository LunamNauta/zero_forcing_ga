#include "random_sampler.hpp"
#include "graph.hpp"

RandomSampler::RandomSampler(const Graph *gi) :
  graph(gi),
  cv(gi->get_order(), 0),
  nf(0),
  gen(std::random_device{}())
{}

double RandomSampler::get_weight(Vertex u) const {
  return RS_EPSILON + static_cast<double>(cv[u]) / nf;
}

double RandomSampler::sum_weights(const VertexBitset &ignored) const {
  double sum = 0;
  for (Vertex u = 0; u < graph->get_order(); u++){
    if (ignored[u]) continue;
    sum += get_weight(u);
  }
  return sum;
}

double RandomSampler::sum_weights(const VertexSet &ignored) const {
  double sum = 0;
  for (Vertex u = 0; u < graph->get_order(); u++){
    if (ignored.find(u) != ignored.cend()) continue;
    sum += get_weight(u);
  }
  return sum;
}

void RandomSampler::update_weights(const VertexBitset &fort) {
  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (!fort[u]) continue;
    cv[u]++;
  }
  nf++;
}

void RandomSampler::update_weights(const VertexSet &fort) {
  for (Vertex u : fort) cv[u]++;
  nf++;
}

VertexBitset RandomSampler::sample_bitset(std::size_t num_samples, std::size_t max_attempts, const VertexBitset &ignored) {
  if (num_samples == 0 || graph->get_order() == 0) return {};
  num_samples = std::min(num_samples, graph->get_order());

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<std::pair<Vertex, double>> candidates;  
  double total_weight = sum_weights(ignored);

  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (ignored[u]) continue;
    double base = dist(gen);
    double weight = get_weight(u);
    candidates.emplace_back(u, std::pow(base, 1.0 / (weight / total_weight)));
  }

  std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b){
    return a.second > b.second;
  });

  VertexBitset sample(graph->get_order(), false);
  for (std::size_t a = 0; a < std::min(candidates.size(), num_samples); a++) {
    sample[candidates[a].first] = true;
  }
  return sample;
}

VertexSet RandomSampler::sample_set(std::size_t num_samples, std::size_t max_attempts, const VertexSet &ignored) {
  if (num_samples == 0 || graph->get_order() == 0) return {};
  num_samples = std::min(num_samples, graph->get_order());

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<std::pair<Vertex, double>> candidates;  
  double total_weight = sum_weights(ignored);

  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (ignored.find(u) != ignored.cend()) continue;
    double base = dist(gen);
    double weight = get_weight(u);
    candidates.emplace_back(u, std::pow(base, 1.0 / (weight / total_weight)));
  }

  std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b){
    return a.second > b.second;
  });

  VertexSet sample;
  for (std::size_t a = 0; a < std::min(candidates.size(), num_samples); a++) {
    sample.insert(candidates[a].first);
  }
  return sample;
}