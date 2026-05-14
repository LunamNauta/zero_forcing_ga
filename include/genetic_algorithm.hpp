#ifndef GENETIC_ALGORITHM_HEADER
#define GENETIC_ALGORITHM_HEADER

#include <random>

#include "random_sampler.hpp"
#include "graph.hpp"

typedef VertexBitset Chromosome;

inline constexpr double GA_ELITE_PCT = 0.15;

struct Individual {
  Chromosome genes;
  double score;
};

class GeneticSolver {
private:
  const Graph *graph;
  std::vector<Individual> population;
  RandomSampler sampler;
  std::mt19937 gen;

  Individual run(std::size_t num_generations);

  void initialize_population();
  void score_population();
  void crossover_population();
  void mutate_population();
  void fix_population();

  Individual& select_parent();
  Individual crossover_individual(const Individual &parent1, const Individual &parent2);
  void mutate_individual(Individual &ind);
  void fix_individual(Individual &ind);

  static VertexBitset chromosome_to_bitset(const Chromosome &genes);
  static VertexSet chromosome_to_set(const Chromosome &genes);

public:
  GeneticSolver(const Graph *gi, std::size_t population_size);

  VertexBitset run_bitset(std::size_t num_generations);
  VertexSet run_set(std::size_t num_generations);
};

#endif