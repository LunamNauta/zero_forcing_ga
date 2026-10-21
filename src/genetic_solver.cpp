#include "genetic_solver.hpp"
#include "weighted_sampler.hpp"
#include "graph.hpp"

#include <functional>
#include <algorithm>
#include <numeric>
#include <utility>

// ------------------------------------------------------------------------------------------------

Individual::Individual(const Graph *gi, const WeightedSampler *si, const VertexBitset &ii) :
  graph(gi),
  sampler(si),
  _initial(gi->order()),
  closure(gi->order()),
  pt(0),
  dirty(true)
{
  initial(ii);
}

Individual::Individual(const Graph *gi, const WeightedSampler *si) :
  graph(gi),
  sampler(si),
  _initial(gi->order()),
  closure(gi->order()),
  pt(0),
  dirty(false)
{
  _initial.set();
  closure.set();
}

void Individual::update() {
  closure = _initial;
  pt = zero_forcing_closure(*graph, closure);
  dirty = false;
}

VertexBitset Individual::initial(const VertexBitset &ii) {
  _initial = ii;
  _initial.resize(graph->order(), false);
  dirty = true;
  return _initial;
}

VertexBitset Individual::initial() const {
  return _initial;
}

double Individual::score() {
  const std::size_t n = graph->order();
  const std::size_t z = _initial.count();
  const double size_term = static_cast<double>(n - (z > n ? n : z));
  constexpr double alpha = 1.0;
  constexpr double beta = 0.1;

  if (n == 0) return 0.0;

  double total_weight = 0.0;
  double selected_weight = 0.0;
  for (Vertex u = 0; u < n; ++u) {
    const double w = sampler->get_weight(u);
    if (!std::isfinite(w) || w <= 0.0) {
      // Treat an invalid weight vector as the paper's uniform initial state.
      return alpha * size_term + beta * static_cast<double>(z) / static_cast<double>(n);
    }
    total_weight += w;
    if (_initial.test(u)) selected_weight += w;
  }
  if (!(total_weight > 0.0) || !std::isfinite(total_weight) ||
      !std::isfinite(selected_weight)) {
    return alpha * size_term + beta * static_cast<double>(z) / static_cast<double>(n);
  }
  return alpha * size_term + beta * (selected_weight / total_weight);
}

std::size_t Individual::z() {
  return _initial.count();
}

bool Individual::forces() {
  if (dirty) update();
  return closure.count() == graph->order();
}

// ------------------------------------------------------------------------------------------------

GeneticSolverInfo::GeneticSolverInfo(const Graph *gi, const WeightedSampler *si) :
  _best_individual(gi, si),
  _best_variance(1.0),
  _best_score(0),
  _best_z(gi->order()),
  _variance(1.0),
  _since_better_variance(0),
  _since_better_score(0),
  _since_better_z(0),
  _upper_bound(gi->order()),
  _lower_bound(1)
{}

Individual GeneticSolverInfo::best_individual() const {
  return _best_individual;
}

double GeneticSolverInfo::best_variance() const {
  return _best_variance;
}

double GeneticSolverInfo::best_score() const {
  return _best_score;
}

std::size_t GeneticSolverInfo::best_z() const {
  return _best_z;
}

double GeneticSolverInfo::variance() const {
  return _variance;
}

std::size_t GeneticSolverInfo::since_better_variance() const {
  return _since_better_variance;
}

std::size_t GeneticSolverInfo::since_better_score() const {
  return _since_better_score;
}

std::size_t GeneticSolverInfo::since_better_z() const {
  return _since_better_z;
}

std::pair<std::size_t, std::size_t> GeneticSolverInfo::bound() const {
  return {_lower_bound, _upper_bound};
}

std::size_t GeneticSolverInfo::upper_bound() const {
  return _upper_bound;
}

std::size_t GeneticSolverInfo::lower_bound() const {
  return _lower_bound;
}

// ------------------------------------------------------------------------------------------------

VertexBitset minimum_hitting_set(const Graph &graph, const std::unordered_set<VertexBitset> &known_forts) {
  const std::size_t n = graph.order();
  const std::size_t m = known_forts.size();

  VertexBitset best(n);
  if (m == 0) {
    return best;
  }

  // Store only the vertices belonging to each fort.
  std::vector<std::vector<Vertex>> members;
  members.reserve(m);

  // Inverted index: forts_by_vertex[v] lists the forts containing v.
  std::vector<std::vector<std::size_t>> forts_by_vertex(n);

  for (const VertexBitset &fort : known_forts) {
    const std::size_t id = members.size();
    members.emplace_back();

    for (Vertex v = 0; v < n; ++v) {
      if (fort.test(v)) {
        members.back().push_back(v);
        forts_by_vertex[v].push_back(id);
      }
    }

    // No hitting set can hit an empty fort.
    if (members.back().empty()) {
      return best;
    }
  }

  // Number of selected vertices in each fort.
  std::vector<std::size_t> hit_count(m, 0);
  std::vector<bool> selected(n, false);

  // Build a greedy feasible solution to provide an upper bound.
  VertexBitset greedy(n);
  std::size_t greedy_count = 0;

  while (true) {
    std::size_t uncovered = 0;
    std::vector<std::size_t> coverage(n, 0);

    for (std::size_t f = 0; f < m; ++f) {
      if (hit_count[f] != 0) {
        continue;
      }

      ++uncovered;
      for (Vertex v : members[f]) {
        if (!selected[v]) {
          ++coverage[v];
        }
      }
    }

    if (uncovered == 0) {
      best = greedy;
      break;
    }

    Vertex best_vertex = n;
    std::size_t best_coverage = 0;

    for (Vertex v = 0; v < n; ++v) {
      if (!selected[v] && coverage[v] > best_coverage) {
        best_coverage = coverage[v];
        best_vertex = v;
      }
    }

    // Defensive guard: an uncovered fort must have an unselected vertex.
    if (best_vertex == n) {
      return VertexBitset(n);
    }

    selected[best_vertex] = true;
    greedy.set(best_vertex);
    ++greedy_count;

    for (std::size_t f : forts_by_vertex[best_vertex]) {
      ++hit_count[f];
    }
  }

  std::size_t best_count = greedy_count;

  // Start exact search from the empty set.
  std::fill(hit_count.begin(), hit_count.end(), 0);
  std::fill(selected.begin(), selected.end(), false);

  VertexBitset chosen(n);

  std::function<void(std::size_t)> search =
      [&](std::size_t chosen_count) {
    if (chosen_count >= best_count) {
      return;
    }

    // Choose the smallest currently unhit fort.
    std::size_t target = m;
    std::size_t target_size = n + 1;
    std::size_t uncovered = 0;

    for (std::size_t f = 0; f < m; ++f) {
      if (hit_count[f] != 0) {
        continue;
      }

      ++uncovered;

      if (members[f].size() < target_size) {
        target = f;
        target_size = members[f].size();
      }
    }

    if (uncovered == 0) {
      best = chosen;
      best_count = chosen_count;
      return;
    }

    if (chosen_count + 1 >= best_count) {
      // Every remaining solution needs at least one more vertex.
      return;
    }

    // Admissible lower bound:
    // every remaining unhit fort must be hit by a future vertex.
    // One vertex can hit at most max_coverage of these forts.
    std::vector<std::size_t> coverage(n, 0);

    for (std::size_t f = 0; f < m; ++f) {
      if (hit_count[f] != 0) {
        continue;
      }

      for (Vertex v : members[f]) {
        if (!selected[v]) {
          ++coverage[v];
        }
      }
    }

    std::size_t max_coverage = 0;
    for (Vertex v = 0; v < n; ++v) {
      if (!selected[v] && coverage[v] > max_coverage) {
        max_coverage = coverage[v];
      }
    }

    if (max_coverage == 0) {
      return;
    }

    const std::size_t lower_bound =
        (uncovered + max_coverage - 1) / max_coverage;

    if (chosen_count + lower_bound >= best_count) {
      return;
    }

    // Every hitting set must choose a vertex from the target fort.
    for (Vertex v : members[target]) {
      selected[v] = true;
      chosen.set(v);

      for (std::size_t f : forts_by_vertex[v]) {
        ++hit_count[f];
      }

      search(chosen_count + 1);

      for (std::size_t f : forts_by_vertex[v]) {
        --hit_count[f];
      }

      chosen.reset(v);
      selected[v] = false;
    }
  };

  search(0);
  return best;
}

bool population_sorter(Individual &a, Individual &b) {
  const bool a_forces = a.forces();
  const bool b_forces = b.forces();
  if (a_forces != b_forces) return a_forces;
  return a.score() > b.score();
}

GeneticSolver::GeneticSolver(const Graph *gi, std::size_t psi) :
  graph(gi),
  sampler(gi),
  _info(gi, &sampler),
  population(std::max<std::size_t>(3, psi), Individual(gi, &sampler)),
  gen(std::random_device{}()),
  known_forts({}),
  dead(true)
{
  initialize_population();
  repair_all_individuals();
  trim_all_individuals();
  trim_all_forts();
  acknowledge_all_forts();
  
  _info._best_individual = population.front();
  _info._best_variance = variance();
  _info._best_score = population.front().score();
  _info._best_z = population.front().z();
  force_upper_bound(_info._best_z);
  _info._variance = variance();
}

// MATCHES SOURCE (2.4)
void GeneticSolver::initialize_population() {
  double step = static_cast<double>(_info._upper_bound - _info._lower_bound) / (population.size() - 1);
  for (std::size_t a = 0; a < population.size(); a++) {
    std::size_t k = _info._lower_bound + a * step;
    VertexBitset subset = sampler.sample_bitset(k, VertexBitset(graph->order()));
    population[a].initial(subset);
    repair_individual(population[a]);
    trim_individual(population[a]);
  }
  std::sort(population.begin(), population.end(), population_sorter);
}

// MATCHES SOURCE (3.1)
std::pair<Individual&, Individual&> GeneticSolver::select_parents() {
  const std::size_t TOURNAMENT_SIZE = 3;
  const double GAMMA = 1.0;
  const double DELTA = 0.1;

  std::size_t tournament_size = std::min(TOURNAMENT_SIZE, graph->order());
  std::uniform_int_distribution<std::size_t> dist(0, population.size() - 1);
  std::vector<std::size_t> indices;
  indices.reserve(tournament_size);

  while (indices.size() < tournament_size) {
    std::size_t candidate = dist(gen);
    if (std::find(indices.begin(), indices.end(), candidate) != indices.end()) continue;
    indices.push_back(candidate);
  }

  Individual *curr_parent1;
  Individual *curr_parent2;
  double curr_score = 0;

  for (std::size_t a = 0; a < indices.size(); a++) {
    for (std::size_t b = a + 1; b < indices.size(); b++) {
      Individual &ind1 = population[indices[a]];
      Individual &ind2 = population[indices[b]];
      double fitness = GAMMA / _info.best_score() * (ind1.score() + ind2.score());
      double andr = (ind1.initial() & ind2.initial()).count();
      double orr = (ind1.initial() | ind2.initial()).count();
      double jaccard = DELTA * (1 - andr / orr);
      double score = fitness + jaccard;

      if (score <= curr_score) continue;
      curr_parent1 = &ind1;
      curr_parent2 = &ind2;
      curr_score = score;
    }
  }

  return {*curr_parent1, *curr_parent2};
}

// MATCHES SOURCE (3.2)
Individual GeneticSolver::crossover_individual(Individual &ind1, Individual &ind2) {
  const std::size_t R = 3;

  Individual best_child(graph, &sampler);
  Individual child(graph, &sampler);

  VertexBitset sample_space0 = ind1.initial() | ind2.initial();
  std::uniform_int_distribution<std::size_t> dist(1, graph->order() - sample_space0.count());
  VertexBitset epsilon = sampler.sample_bitset(dist(gen), sample_space0);
  VertexBitset sample_space1 = sample_space0 | epsilon;

  double step = static_cast<double>(_info._upper_bound - _info._lower_bound) / R;
  for (std::size_t a = 0; a < R; a++) {
    std::size_t k = _info._lower_bound + a * step;
    child.initial(sampler.sample_bitset(k, ~sample_space1));
    repair_individual(child);
    trim_individual(child);
    if (child.score() > best_child.score()) best_child = child;
  }

  return best_child;
}

// MATCHES SOURCE ~ (3.3)
// Uses a brute force hitting set problem,
// rather than an IP model for finding the
// final child. Currently, it does not even
// do that. My rationale is that it is too slow,
// and shows no improvement. Could be that my hitting
// set calculator is just not working
void GeneticSolver::crossover_all_individuals() {
  std::vector<Individual> next_population;
  next_population.reserve(population.size());

  std::size_t elites = std::max(1.0, population.size() * 0.1);
  next_population.insert(next_population.end(), population.begin(), population.begin() + elites);

  std::size_t offspring_count = population.size() - elites - 1;
  for (std::size_t a = 0; a < offspring_count; a++) {
    auto [parent1, parent2] = select_parents();
    Individual child = crossover_individual(parent1, parent2);
    next_population.emplace_back(child);
  }

  trim_all_forts();

  VertexBitset cover = minimum_hitting_set(*graph, known_forts);
  force_lower_bound(cover.count());

  Individual cover_ind(graph, &sampler, cover);
  repair_individual(cover_ind);
  next_population.emplace_back(cover_ind);

  acknowledge_all_forts();

  population = std::move(next_population);
  std::sort(population.begin(), population.end(), population_sorter);
  for (Individual &individual : population) {
    if (individual.z() >= _info._upper_bound) continue;
    _info._upper_bound = individual.z();
  }
}

// MATCHES SOURCE ~ (2.6)
// Does not trim the child at the end
// These are distinct operations
void GeneticSolver::repair_individual(Individual &ind) {
  if (ind.forces()) return;

  VertexBitset current(ind.initial());
  VertexBitset closure(current);
  zero_forcing_closure(*graph, closure);
  if (closure.count() == graph->order()) return;
  VertexBitset fort(~closure);

  std::vector<Vertex> indices(graph->order());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), gen);

  for (Vertex u : indices) {
    if (closure.test(u)) continue;

    current.set(u);
    closure = current;
    zero_forcing_closure(*graph, closure);
    if (closure.count() != graph->order()) fort = ~closure;
    else break;
  }

  known_forts.insert(fort);
  ind.initial(current);
}

void GeneticSolver::repair_all_individuals() {
  for (Individual &ind : population) repair_individual(ind);
}

// MATCHES SOURCE (2.7)
void GeneticSolver::trim_individual(Individual &ind) {
  if (!ind.forces()) return;

  VertexBitset current(ind.initial());
  VertexBitset closure;

  std::vector<Vertex> indices(graph->order());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), gen);
  
  for (Vertex u : indices) {
    if (!current.test(u)) continue;

    current.reset(u);
    closure = current;
    zero_forcing_closure(*graph, closure);
    if (closure.count() == graph->order()) continue;
    else current.set(u);
  }
  
  ind.initial(current);
}

void GeneticSolver::trim_all_individuals() {
  for (Individual &ind : population) trim_individual(ind);
}

void GeneticSolver::trim_fort(VertexBitset &fort) {
  VertexBitset current(~fort);
  VertexBitset closure;

  std::vector<Vertex> vertices(graph->order());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::shuffle(vertices.begin(), vertices.end(), gen);

  for (Vertex u : vertices) {
    if (current.test(u)) continue;

    current.set(u);
    closure = current; 
    zero_forcing_closure(*graph, closure);
    if (closure.count() == graph->order()) current.reset(u);
    else current = closure;
  }

  fort = ~current;
}

void GeneticSolver::trim_all_forts() {
  std::unordered_set<VertexBitset> next_forts;

  for (const VertexBitset &fort : known_forts) {
    VertexBitset minimal_fort = fort;
    trim_fort(minimal_fort);
    if (!minimal_fort.any()) continue;
    next_forts.emplace(minimal_fort);
  }

  known_forts = std::move(next_forts);
}

void GeneticSolver::acknowledge_fort(const VertexBitset &fort) {
  sampler.update_weights(fort);
}

void GeneticSolver::acknowledge_all_forts() {
  for (const VertexBitset &fort : known_forts) acknowledge_fort(fort);
  known_forts.clear();
}

void GeneticSolver::update_info() {
  Individual &current_best = population.front();
  if (current_best.score() > _info._best_score) {
    _info._best_score = current_best.score();
    _info._since_better_score = 0;
  } else _info._since_better_score++;

  if (current_best.forces() && current_best.z() < _info._best_z) {
    _info._best_z = current_best.z();
    _info._best_individual = current_best;
    _info._since_better_z = 0;
    } else _info._since_better_z++;

  for (Individual &ind : population) {
    if (ind.forces()) force_upper_bound(ind.z());
  }

  const double current_variance = variance();
  if (current_variance < _info._best_variance) {
    _info._best_variance = current_variance;
    _info._since_better_variance = 0;
  } else _info._since_better_variance++;
}

double GeneticSolver::variance() const {
  if (population.empty() || graph->order() == 0) return 0.0;

  double total_variance = 0.0;
  const std::size_t pop_size = population.size();
  const std::size_t n = graph->order();

  for (Vertex u = 0; u < n; ++u) {
    std::size_t bit_count = 0;
    for (const Individual &ind : population) {
      if (ind.initial().test(u)) ++bit_count;
    }
    const double p = static_cast<double>(bit_count) / static_cast<double>(pop_size);
    total_variance += p * (1.0 - p);
  }
  return total_variance / static_cast<double>(n);
}

void GeneticSolver::force_upper_bound(std::size_t ub) {
  ub = std::max(_info._lower_bound, std::min(ub, _info._upper_bound));
  _info._upper_bound = std::min(ub, _info._upper_bound);
}

void GeneticSolver::force_lower_bound(std::size_t lb) {
  lb = std::min(std::max(lb, _info._lower_bound), _info._upper_bound);
  _info._lower_bound = std::max(lb, _info._lower_bound);
}

GeneticSolverInfo GeneticSolver::info() const {
  return _info;
}

void GeneticSolver::run(std::size_t generations) {
  dead = false;
  for (std::size_t generation = 0; !dead && generation < generations; ++generation) {
    crossover_all_individuals();
    update_info();
  }
  dead = true;
}

void GeneticSolver::kill() {
  dead = true;
}
