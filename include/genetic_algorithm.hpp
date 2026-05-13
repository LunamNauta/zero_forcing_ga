#ifndef GENETIC_ALGORITHM_HEADER
#define GENETIC_ALGORITHM_HEADER

#include <unordered_set>
#include <vector>

#include "graph.hpp"

typedef std::vector<bool> Chromosome;

struct Individual {
    Chromosome genes;
    double fitness;
    std::size_t size;
};

class ZFGeneticSolver {
public:
  ZFGeneticSolver(const Graph *gi, std::size_t population);

  std::unordered_set<VertexIndex> run(std::size_t generations);

private:
  const Graph *graph;
  std::vector<Individual> population;

  void initialize_population();
  void score_population();
  void crossover_population();
  void mutate_population();

  static Individual crossover_individual(const Individual &ind1, const Individual &ind2);
  static void mutate_individual(Individual &ind);

  static std::unordered_set<VertexIndex> chromosome_to_set(const Chromosome& genes);
};

#endif