#include "genetic_algorithm.hpp"

#include <random>

#include "random_sampler.hpp"
#include "zero_forcing.hpp"
#include "graph.hpp"

GeneticSolver::GeneticSolver(const Graph *gi, std::size_t population_size) :
  graph(gi),
  population(population_size),
  sampler(gi),
  gen(std::random_device{}())
{
  initialize_population();
}

void GeneticSolver::initialize_population() {
  std::bernoulli_distribution dist(0.85);

  for (Individual &ind : population) {
    ind.genes.resize(graph->get_order(), false);
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (!dist(gen)) continue;
      ind.genes[u] = true;
    }
  }
}

void GeneticSolver::score_population() {
  double penalty_weight = 1.0;

  for (Individual &ind : population) {
    VertexBitset forced(ind.genes);
    zero_forcing_closure(*graph, forced);

    double base_score = graph->get_order() - ind.genes.count(); 
    double unforced_count = graph->get_order() - forced.count();
    double penalty = penalty_weight * unforced_count;
    
    ind.score = 1000.0 + base_score - penalty;
  }
}

void GeneticSolver::crossover_population() {
  std::sort(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
    return a.score > b.score;
  });

  std::vector<Individual> next_generation;
  next_generation.reserve(population.size());

  std::size_t elites = population.size() * GA_ELITE_PCT;
  for (std::size_t a = 0; a < elites; a++) {
    next_generation.push_back(population[a]);
  }

  while (next_generation.size() < population.size()) {
    Individual &parent1 = select_parent();
    Individual &parent2 = select_parent();
    next_generation.push_back(crossover_individual(parent1, parent2));
  }

  population = std::move(next_generation);
}

void GeneticSolver::mutate_population() {
  std::size_t elites = population.size() * GA_ELITE_PCT;
  for (std::size_t a = elites; a < population.size(); a++) {
    mutate_individual(population[a]);
  }
}

void GeneticSolver::fix_population() {
  std::size_t elites = population.size() * GA_ELITE_PCT;
  for (std::size_t a = elites; a < population.size(); a++) {
    fix_individual(population[a]);
  }
}

Individual& GeneticSolver::select_parent() {
  std::uniform_int_distribution<std::size_t> dist(0, population.size() - 1);
  const std::size_t tournament_size = 3;

  Individual &best = population[dist(gen)];
  for (std::size_t a = 1; a < tournament_size; a++) {
    Individual &candidate = population[dist(gen)];
    if (candidate.score <= best.score) continue;
    best = candidate;
  }

  return best;
}

Individual GeneticSolver::crossover_individual(const Individual &parent1, const Individual &parent2) {
  Individual offspring;
  offspring.genes = parent1.genes & parent2.genes;
  return offspring;
}

void GeneticSolver::mutate_individual(Individual &ind) {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  double prob = 1.0 / graph->get_order();

  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (dist(gen) >= prob) continue;
    ind.genes[u] = !ind.genes[u];
  }
}

void GeneticSolver::fix_individual(Individual &ind) {
  VertexBitset ignored(ind.genes);
  VertexBitset forced(ind.genes);
  zero_forcing_closure(*graph, forced);
  
  while (forced.count() < graph->get_order()) {
    VertexBitset fort = ~forced;
    sampler.update_weights(fort);

    VertexBitset verts = sampler.sample_bitset((graph->get_order() - forced.count()) * 0.1 + 1, 100, ignored);
    ind.genes |= verts;
    ignored |= verts;
    forced |= verts;
    zero_forcing_closure(*graph, forced);
  }
}

VertexBitset GeneticSolver::chromosome_to_bitset(const Chromosome &genes) {
  return genes;
}

VertexSet GeneticSolver::chromosome_to_set(const Chromosome &genes) {
  VertexSet output;
  for (Vertex u = 0; u < genes.size(); u++) {
    if (!genes[u]) continue;
    output.insert(u);
  }
  return output;
}

Individual GeneticSolver::run(std::size_t num_generations) {
  for (std::size_t a = 0; a < num_generations; a++) {
    score_population();
    crossover_population();
    mutate_population();
    fix_population();
  }

  score_population();
  auto best = std::max_element(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
    return a.score < b.score;
  });

  return *best;
}

VertexBitset GeneticSolver::run_bitset(std::size_t num_generations) {
  return run(num_generations).genes;
}

VertexSet GeneticSolver::run_set(std::size_t num_generations) {
  return chromosome_to_set(run(num_generations).genes);
}