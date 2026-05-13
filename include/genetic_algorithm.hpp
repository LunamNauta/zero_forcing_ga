#ifndef GENETIC_ALGORITHM_HEADER
#define GENETIC_ALGORITHM_HEADER

#include <vector>

#include "graph.hpp"
#include "random_sampler.hpp"

typedef std::vector<bool> Chromosome;

struct Individual {
    Chromosome genes;
    double fitness;
    std::size_t size;
};

class ZFGeneticSolver {
private:
  const Graph *graph;
  std::vector<Individual> population;
  RandomSampler sampler;

  void initialize_population();
  void score_population();
  void crossover_population();
  void mutate_population();

  Individual crossover_individual(const Individual &ind1, const Individual &ind2);
  void mutate_individual(Individual &ind);

  static VertexSet chromosome_to_set(const Chromosome& genes);

public:
  ZFGeneticSolver(const Graph *gi, std::size_t population_size);

  VertexSet run(std::size_t generations);
};

#endif