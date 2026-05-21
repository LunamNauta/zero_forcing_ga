#include "genetic_algorithm.hpp"

#include <random>

#include "zero_forcing.hpp"
#include "graph.hpp"

// Initialization ----------------------------------------------------------------
GeneticSolver::GeneticSolver(const Graph *gi, std::size_t psi) :
  graph(gi),
  population(psi),
  min_mutation(0.01),
  max_mutation(0.8),
  mutation_pct(max_mutation),
  min_elite(0.01),
  max_elite(0.2),
  elite_pct(max_elite),
  sampler(gi),
  gen(std::random_device{}()),
  _lower_bound(std::min(1UL, gi->get_order())),
  _upper_bound(gi->get_order()),
  _since_better_score(0),
  _since_better_z(0),
  _best_score(0),
  _best_z(0)
{
  initialize_population();
  fix_population();
  reduce_population();
  score_population();

  _best_score = population[0].score;
  if (population[0].forces) {
    _best_z = population[0].gene.count();
    _best_ind = population[0];
  }
  else {
    _best_ind.gene.resize(graph->get_order());
    _best_ind.gene.set();
    _best_ind.forces = true;
    _best_ind.dirty = false;
    _best_z = _best_ind.gene.count();
    score_individual(_best_ind);
  }
}

void GeneticSolver::initialize_population() {
  double step = static_cast<double>(graph->get_order()) / population.size();
  double order = 1;
  VertexBitset initial;
 
  for (std::size_t a = 0; a < population.size(); a++, order += step) {
    initial = std::move(sampler.sample_bitset(order));
    population[a].gene = initial;
    population[a].dirty = true;
  }
}

// Population Functions ----------------------------------------------------------------
void GeneticSolver::crossover_population() {
  std::vector<Individual> next_population;
  next_population.reserve(population.size());

  std::size_t elites = std::max(1.0, population.size() * elite_pct);
  next_population.insert(next_population.end(), population.begin(), population.begin() + elites);

  while (next_population.size() < population.size()) {
    auto [parent1, parent2] = select_parents();
    Individual child = crossover_individual(parent1, parent2);
    next_population.emplace_back(child);
  }

  population = std::move(next_population);
}

void GeneticSolver::score_population() {
  for (Individual &ind : population) {
    score_individual(ind);
  }
  std::sort(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
    return std::tie(a.forces, a.score) > std::tie(b.forces, b.score);
  });
}

void GeneticSolver::fix_population() {
  for (Individual &ind : population) {
    fix_individual(ind);
  }
}

void GeneticSolver::reduce_population() {
  for (Individual &ind : population) {
    reduce_individual(ind);
  }
}

// Individual Functions ----------------------------------------------------------------
Individual GeneticSolver::crossover_individual(const Individual &ind1, const Individual &ind2) {
  VertexBitset sample;
  Individual best_child;
  Individual child;
  std::size_t min_order = std::min(ind1.gene.count(), ind2.gene.count());

  VertexBitset sample_space0 = ind1.gene | ind2.gene;

  VertexBitset epsilon = sampler.sample_bitset((graph->get_order() - sample_space0.count())*mutation_pct, sample_space0);
  /*
  // VertexBitset epsilon = sampler.sample_bitset((graph->get_order() - sample_space0.count())*mutation_pct, sample_space0);
  std::bernoulli_distribution dist(mutation_pct);
  VertexBitset epsilon(graph->get_order(), false);
  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (!dist(gen)) continue;
    epsilon.set(u);
  }
  */

  VertexBitset sample_space1 = sample_space0 | epsilon;

  best_child.gene = sample_space1;
  best_child.forces = false;
  best_child.dirty = true;

  for (std::size_t a = _lower_bound; a < min_order; a++) {
    sample = std::move(sampler.sample_bitset(a, ~sample_space1));

    child.gene = sample;
    child.dirty = true;
    score_individual(child);

    if (child.forces && child.score > best_child.score) best_child = child;
  }
 
  reduce_individual(best_child);
  score_individual(best_child);

  return best_child;
}

void GeneticSolver::score_individual(Individual &ind) {
  const double penalty_weight = 0.75;
  const double base_weight = 0.25;

  VertexBitset forced;
  std::size_t ptime;
  
  if (ind.dirty || !ind.forces) {
    forced = ind.gene;
    ind.ptime = zero_forcing_closure(*graph, forced);
    acknowledge_fort(~forced);
    ind.forces = forced.count() == graph->get_order();
    ind.dirty = false;
  }
  else forced.resize(graph->get_order(), true);

  double base_score = graph->get_order() - ind.gene.count(); 
  double penalty_score = graph->get_order() - forced.count();
  // ind.score = (base_score * base_weight) - (penalty_score * penalty_weight);
  ind.score = (base_score*base_weight) / (1.0 + penalty_score*penalty_weight);
}

void GeneticSolver::fix_individual(Individual &ind) {
  if (!ind.dirty && ind.forces) return;

  VertexBitset forced(ind.gene);
  VertexBitset sample;

  while (forced.count() < graph->get_order()) {
    zero_forcing_closure(*graph, forced);
    if (forced.count() == graph->get_order()) break;
    acknowledge_fort(~forced);

    sample = std::move(sampler.sample_bitset(1, forced));
    ind.gene |= sample;
    forced |= sample;
  }

  ind.forces = true;
  ind.dirty = false;
}

void GeneticSolver::reduce_individual(Individual &ind) {
  if (ind.dirty || !ind.forces) return;

  std::vector<Vertex> vertices(graph->get_order());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::shuffle(vertices.begin(), vertices.end(), gen);

  VertexBitset forced;
  for (Vertex u : vertices) {
    if (!ind.gene.test(u)) continue;
    ind.gene.reset(u);
    
    forced = ind.gene;
    zero_forcing_closure(*graph, forced);
    if (forced.count() == graph->get_order()) continue;
    acknowledge_fort(~forced);

    ind.gene.set(u);
  }
}

// Selection Functions ----------------------------------------------------------------
std::pair<const Individual&, const Individual&> GeneticSolver::select_parents() {
  const std::size_t tournament_size = std::sqrt(population.size());
  std::uniform_int_distribution<std::size_t> dist(0, population.size() - 1);

  Individual *parent1 = &population[dist(gen)];
  Individual *parent2 = &population[dist(gen)];
  Individual *candidate;

  for (std::size_t a = 0; a < tournament_size; a++) {
    candidate = &population[dist(gen)];
    if (candidate->score <= parent1->score) continue;
    parent1 = candidate;
  }
  for (std::size_t a = 0; a < tournament_size; a++) {
    candidate = &population[dist(gen)];
    if (candidate->score <= parent2->score) continue;
    if (candidate == parent1) continue;
    parent2 = candidate;
  }

  return {*parent1, *parent2};
}

// Fort Acknowledgement ----------------------------------------------------------------
void GeneticSolver::acknowledge_fort(const VertexBitset &fort, bool reduce) {
  VertexBitset smaller_fort = reduce ? reduced_fort(fort) : fort;
  if (!smaller_fort.any()) return;
  sampler.update_weights(smaller_fort);
}

VertexBitset GeneticSolver::reduced_fort(const VertexBitset &fort) {
  std::vector<Vertex> vertices(graph->get_order());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::shuffle(vertices.begin(), vertices.end(), gen);
  VertexBitset closure = ~fort;
  VertexBitset forced;

  for (Vertex u : vertices) {
    if (closure.test(u)) continue;
    closure.set(u);

    forced = closure; 
    zero_forcing_closure(*graph, forced);
    if (forced.count() == graph->get_order()) {
      closure.reset(u);
      continue;
    }
    acknowledge_fort(~forced, false);

    closure = forced;
  }

  return ~closure;
}

// Running ----------------------------------------------------------------
void GeneticSolver::run(std::size_t generations) {
  for (std::size_t a = 0; a < generations; a++) {
    crossover_population();
    score_population();

    if (_since_better_score > std::sqrt(population.size())) {
      mutation_pct *= 1.15;
      elite_pct *= 0.85;
    }
    else {
      mutation_pct *= 0.85;
      elite_pct *= 1.15;
    }
    mutation_pct = std::max(min_mutation, std::min(mutation_pct, max_mutation));
    elite_pct = std::max(min_elite, std::min(elite_pct, max_elite));

    _since_better_score++;
    for (Individual &ind : population) {
      if (ind.score <= _best_score) continue;
      _best_score = ind.score;
      _since_better_score = 0;
    }

    _since_better_z++;
    for (Individual &ind : population) {
      if (!ind.forces || ind.gene.count() >= _best_z) continue;
      _best_z = ind.gene.count();
      _upper_bound = _best_z;
      _since_better_z = 0;
      _best_ind = ind;
    }
  }
}

// Introspection ----------------------------------------------------------------
std::size_t GeneticSolver::since_better_score() const {
  return _since_better_score;
}

std::size_t GeneticSolver::since_better_z() const {
  return _since_better_z;
}

double GeneticSolver::best_score() const {
  return _best_score;
}

std::size_t GeneticSolver::best_z() const {
  return _best_z;
}

Individual GeneticSolver::best_individual() const {
  return _best_ind;
}

std::pair<std::size_t, std::size_t> GeneticSolver::bound() const {
  return {_lower_bound, _upper_bound};
}