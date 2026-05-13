#include "genetic_algorithm.hpp"

#include <algorithm>
#include <random>

#include "random_sampler.hpp"
#include "zero_forcing.hpp"
#include "graph.hpp"

ZFGeneticSolver::ZFGeneticSolver(const Graph *gi, std::size_t population_size) :
  graph(gi),
  population(population_size),
  sampler(gi)
{
  initialize_population();
}

void ZFGeneticSolver::initialize_population() {
  std::mt19937 gen(std::random_device{}());
  std::bernoulli_distribution dist(0.15);

  for (auto &ind : population) {
    ind.genes.assign(graph->get_order(), false);
    for (std::size_t a = 0; a < graph->get_order(); a++) {
      ind.genes[a] = dist(gen);
    }
  }
}

void ZFGeneticSolver::score_population() {
  for (auto &ind : population) {
    VertexSet start_nodes = chromosome_to_set(ind.genes);
    ind.size = start_nodes.size();

    VertexSet forced_nodes = zero_forcing_closure(*graph, start_nodes);
    std::size_t forced_count = forced_nodes.size();

    if (forced_count == graph->get_order()) ind.fitness = 1000.0 + (double)(graph->get_order() - ind.size);
    else ind.fitness = (double)forced_count / graph->get_order();
  }
}

Individual ZFGeneticSolver::crossover_individual(const Individual &ind1, const Individual &ind2) {
  Individual offspring;
  offspring.genes.assign(graph->get_order(), false);
    
  std::vector<Vertex> parent_union;

  for (std::size_t a = 0; a < graph->get_order(); a++) {
    if (ind1.genes[a] && ind2.genes[a]) offspring.genes[a] = true;
    if (ind1.genes[a] || ind2.genes[a]) parent_union.push_back(a);
  }

  VertexSet current_start = chromosome_to_set(offspring.genes);

  VertexSet forced = zero_forcing_closure(*graph, current_start);
    
  while (forced.size() < graph->get_order()) {
    VertexSet verts = sampler(1, 100, forced);
    current_start.insert(*verts.begin());
    forced.insert(*verts.begin());
    forced = std::move(zero_forcing_closure(*graph, forced));
  }

  offspring.size = current_start.size();
  return offspring;
}

void ZFGeneticSolver::mutate_individual(Individual &ind) {
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
    
  double prob = 1.0 / graph->get_order();

  for (std::size_t a = 0; a < ind.genes.size(); a++) {
    if (dist(gen) >= prob) continue;
    ind.genes[a] = !ind.genes[a];
  }
}

void ZFGeneticSolver::crossover_population() {
  std::sort(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
    return a.fitness > b.fitness;
  });

  std::size_t pop_size = population.size();
  std::size_t half = pop_size / 2;
  std::vector<Individual> next_gen;
  next_gen.reserve(pop_size);

  std::size_t elites = pop_size / 10;
  for (std::size_t i = 0; i < elites; ++i) {
    next_gen.push_back(population[i]);
  }

  std::mt19937 gen(std::random_device{}());
  std::uniform_int_distribution<std::size_t> parent_dist(0, half);

  while (next_gen.size() < pop_size) {
    next_gen.push_back(crossover_individual(population[parent_dist(gen)], population[parent_dist(gen)]));
  }

  population = std::move(next_gen);
}

void ZFGeneticSolver::mutate_population() {
  std::size_t elites = population.size() / 10;
  for (std::size_t a = elites; a < population.size(); a++) {
    mutate_individual(population[a]);
  }
}

VertexSet ZFGeneticSolver::chromosome_to_set(const Chromosome& genes) {
  VertexSet v_set;
  for (std::size_t a = 0; a < genes.size(); a++) {
    if (genes[a]) v_set.insert(a);
  }
  return v_set;
}

VertexSet ZFGeneticSolver::run(std::size_t generations) {
  for (std::size_t a = 0; a < generations; a++) {
    score_population();
    crossover_population();
    mutate_population();
  }

  score_population();
  auto best = std::max_element(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
    return a.fitness < b.fitness;
  });

  return chromosome_to_set(best->genes);
}