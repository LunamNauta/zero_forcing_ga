#include "../include/weighted_sampler.hpp"

WeightedSampler::WeightedSampler(const Graph *graph) :
  graph(graph),
  cv(graph->order(), 0),
  nf(0),
  generator(std::random_device{}())
{}

double WeightedSampler::sum_weights(const VertexBitset &ignored) const {
  double sum = 0.0;
  for (Vertex u = 0; u < graph->order(); u++){
    if (ignored[u]) continue;
    sum += get_weight(u);
  }
  return sum;
}

double WeightedSampler::sum_weights(const VertexSet &ignored) const {
  double sum = 0.0;
  for (Vertex u = 0; u < graph->order(); u++){
    if (ignored.find(u) != ignored.cend()) continue;
    sum += get_weight(u);
  }
  return sum;
}

double WeightedSampler::get_weight(Vertex u) const {
  return RS_EPSILON + static_cast<double>(cv[u]) / nf;
}

void WeightedSampler::reset_weights() {
  for (double &c : cv) {
    c = 0;
  }
  nf = 0;
}

void WeightedSampler::update_weights(const VertexBitset &fort) {
  for (Vertex u = 0; u < graph->order(); u++) {
    if (!fort[u]) continue;
    cv[u] = RS_ALPHA * cv[u] + 1.0 / fort.count();
  }
  nf++;
}

void WeightedSampler::update_weights(const VertexSet &fort) {
  for (Vertex u : fort) {
    cv[u] = RS_ALPHA * cv[u] + 1.0 / fort.size();
  }
  nf++;
}

VertexBitset WeightedSampler::sample_bitset(std::size_t num_samples, const VertexBitset &ignored, bool invert) {
  // Default conditions for empty samples
  if (num_samples == 0) return VertexBitset(graph->order());

  // Bound number of samples
  num_samples = std::min(num_samples, graph->order());
  VertexBitset ignored_tmp = ignored;
  ignored_tmp.resize(graph->order());

  // Distribution for base of selection weight
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  // Buffer for potential sample vertices (and their weight)
  std::vector<std::pair<Vertex, double>> candidates;  
  double total_weight = sum_weights(ignored_tmp);

  // Compute weights
  for (Vertex u = 0; u < graph->order(); u++) {
    if (ignored_tmp[u]) continue;
    // Vertex weight is rand^(1/(weight / total_weight))
    double base = distribution(generator);
    double weight = get_weight(u);
    double ratio = invert ? weight / total_weight : total_weight / weight;
    candidates.emplace_back(u, std::pow(base, ratio));
  }

  // Sort the potential samples based on their weight
  std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b){
    return a.second > b.second;
  });

  // Select the top @num_samples candidates for the sample
  VertexBitset sample(graph->order());
  for (std::size_t a = 0; a < std::min(candidates.size(), num_samples); a++) {
    sample[candidates[a].first] = true;
  }
  return sample;
}

VertexSet WeightedSampler::sample(std::size_t num_samples, const VertexSet &ignored, bool invert) {
  // Default conditions for empty samples or graph
  if (num_samples == 0) return {};

  // Bound number of samples
  num_samples = std::min(num_samples, graph->order());

  // Distribution for base of selection weight
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  // Buffer for potential sample vertices (and their weight)
  std::vector<std::pair<Vertex, double>> candidates;  
  double total_weight = sum_weights(ignored);

  for (Vertex u = 0; u < graph->order(); u++) {
    if (ignored.find(u) != ignored.end()) continue;
    // Vertex weight is rand^(1/(weight / total_weight))
    double base = distribution(generator);
    double weight = get_weight(u);
    double ratio = invert ? weight / total_weight : total_weight / weight;
    candidates.emplace_back(u, std::pow(base, ratio));
  }

  // Sort the potential samples based on their weight
  std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b){
    return a.second > b.second;
  });

  // Select the top @num_samples candidates for the sample
  VertexSet sample;
  for (std::size_t a = 0; a < std::min(candidates.size(), num_samples); a++) {
    sample.insert(candidates[a].first);
  }
  return sample;
}
