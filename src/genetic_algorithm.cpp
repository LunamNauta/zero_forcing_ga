#include "genetic_algorithm.hpp"

#include <random>
#include <algorithm>
#include <tuple>

#include "zero_forcing.hpp"
#include "graph.hpp"

// Initialization ----------------------------------------------------------------
GeneticSolver::GeneticSolver(const Graph *gi, std::size_t psi) :
  graph(gi),
  population(psi, Individual(gi)),
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
  _best_z(0),
  _best_ind(gi)
{
  initialize_population();
  // fix_population();
  reduce_population();

  std::sort(population.begin(), population.end(), [](Individual &a, Individual &b) {
    bool a_forces = a.forces();
    bool b_forces = b.forces();
    double a_score = a.get_score();
    double b_score = b.get_score();
    return std::tie(a_forces, a_score) > std::tie(b_forces, b_score);
  });

  if (population[0].forces()) _best_ind = population[0];
  _best_score = _best_ind.get_score();
  _best_z = _best_ind.get_z();
}

void GeneticSolver::initialize_population() {
  double step = static_cast<double>(graph->get_order()) / population.size();
  double order = 1;
  VertexBitset initial;
 
  for (std::size_t a = 0; a < population.size(); a++, order += step) {
    initial = std::move(sampler.sample_bitset(order));
    population[a].set_initial(initial);
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
Individual GeneticSolver::crossover_individual(Individual &ind1, Individual &ind2) {
  VertexBitset sample;
  Individual best_child(graph);
  Individual child(graph);
  std::size_t min_order = std::min(ind1.get_z(), ind2.get_z());

  VertexBitset sample_space0 = ind1.get_initial() | ind2.get_initial();
  VertexBitset epsilon = sampler.sample_bitset((graph->get_order() - sample_space0.count())*mutation_pct, sample_space0);
  VertexBitset sample_space1 = sample_space0 | epsilon;

  best_child.set_initial(sample_space1);
  for (std::size_t a = _lower_bound; a < min_order; a++) {
    sample = std::move(sampler.sample_bitset(a, ~sample_space1));
    child.set_initial(sample);
    if (child.forces() && child.get_score() > best_child.get_score()) best_child = child;
  }
  reduce_individual(best_child);

  return best_child;
}

void GeneticSolver::fix_individual(Individual &ind) {
  if (ind.forces()) return;

  VertexBitset current(ind.get_initial());
  VertexBitset forced = current;
  VertexBitset sample;

  while (forced.count() < graph->get_order()) {
    zero_forcing_closure(*graph, forced);
    if (forced.count() == graph->get_order()) break;
    acknowledge_fort(~forced);

    sample = std::move(sampler.sample_bitset(1, forced));
    current |= sample;
    forced |= sample;
  }

  ind.set_initial(current);
}

void GeneticSolver::reduce_individual(Individual &ind) {
  if (!ind.forces()) return;

  std::vector<Vertex> vertices(graph->get_order());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::shuffle(vertices.begin(), vertices.end(), gen);

  VertexBitset current(ind.get_initial());
  VertexBitset forced;
  
  for (Vertex u : vertices) {
    if (!current.test(u)) continue;
    current.reset(u);
    
    forced = current;
    zero_forcing_closure(*graph, forced);
    if (forced.count() == graph->get_order()) continue;
    known_forts.insert(~forced);

    current.set(u);
  }
  
  ind.set_initial(current);
}

// Selection Functions ----------------------------------------------------------------
std::pair<Individual&, Individual&> GeneticSolver::select_parents() {
  const std::size_t tournament_size = 3; //std::sqrt(population.size());
  std::uniform_int_distribution<std::size_t> dist(0, population.size() / 2);

  Individual *parent1 = &population[dist(gen)];
  Individual *parent2 = &population[dist(gen)];
  Individual *candidate;

  for (std::size_t a = 0; a < tournament_size; a++) {
    candidate = &population[dist(gen)];
    if (candidate->get_score() <= parent1->get_score()) continue;
    parent1 = candidate;
  }
  for (std::size_t a = 0; a < tournament_size; a++) {
    candidate = &population[dist(gen)];
    if (candidate->get_score() <= parent2->get_score()) continue;
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
  std::uniform_int_distribution<std::size_t> dist_v;

  for (std::size_t a = 0; a < generations; a++) {
    known_forts.clear();
    crossover_population();
    for (const VertexBitset &fort : known_forts) {
      acknowledge_fort(fort);
    }
    
    std::sort(population.begin(), population.end(), [](Individual &a, Individual &b) {
      bool a_forces = a.forces();
      bool b_forces = b.forces();
      double a_score = a.get_score();
      double b_score = b.get_score();
      return std::tie(a_forces, a_score) > std::tie(b_forces, b_score);
    });

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

    Individual &current_best = population.front();

    if (current_best.get_score() > _best_score) {
      _best_score = current_best.get_score();
      _since_better_score = 0;
    } 
    else _since_better_score++;

    if (current_best.forces() && current_best.get_z() < _best_z) {
      _best_z = current_best.get_z();
      _upper_bound = _best_z;
      _since_better_z = 0;
      _best_ind = current_best;
    } 
    else _since_better_z++;
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

double GeneticSolver::variance() const {
  if (population.empty() || graph->get_order() == 0) return 0.0;

  double total_variance = 0.0;
  std::size_t pop_size = population.size();
  std::size_t num_vertices = graph->get_order();

  for (Vertex u = 0; u < num_vertices; u++) {
    std::size_t bit_count = 0;
    
    for (const Individual &ind : population) {
      if (ind.initial.test(u)) bit_count++;
    }

    double p = static_cast<double>(bit_count) / pop_size;
    
    total_variance += p * (1.0 - p);
  }

  return total_variance / num_vertices;
}