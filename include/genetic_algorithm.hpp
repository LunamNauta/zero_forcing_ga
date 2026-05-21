#ifndef GENETIC_ALGORITHM_HEADER
#define GENETIC_ALGORITHM_HEADER

#include "random_sampler.hpp"
#include "graph.hpp"
#include <random>

typedef VertexBitset Gene;

struct Individual {
  Gene gene;
  double score = 0;
  bool forces = false;
  bool dirty = true;
  std::size_t ptime = std::numeric_limits<std::size_t>::max();
};

class GeneticSolver {
private:
  const Graph *graph;
  std::vector<Individual> population;
  double min_mutation;
  double max_mutation;
  double mutation_pct;
  double min_elite;
  double max_elite;
  double elite_pct;
  RandomSampler sampler;
  std::mt19937 gen;

  std::size_t _lower_bound;
  std::size_t _upper_bound;

  std::size_t _since_better_score;
  std::size_t _since_better_z;
  double _best_score;
  std::size_t _best_z;
  Individual _best_ind;

  // Initialization ----------------------------------------------------------------
  void initialize_population();

  // Population Functions ----------------------------------------------------------------
  void crossover_population();
  void score_population();
  void fix_population();
  void reduce_population();

  // Individual Functions ----------------------------------------------------------------
  Individual crossover_individual(const Individual &ind1, const Individual &ind2);
  void score_individual(Individual &ind);
  void fix_individual(Individual &ind);
  void reduce_individual(Individual &ind);

  // Selection Functions ----------------------------------------------------------------
  std::pair<const Individual&, const Individual&> select_parents();

  // Fort Acknowledgement ----------------------------------------------------------------
  void acknowledge_fort(const VertexBitset &fort, bool reduce = true);
  VertexBitset reduced_fort(const VertexBitset &fort);

public:
  // Initialization ----------------------------------------------------------------
  GeneticSolver(const Graph *gi, std::size_t psi);

  // Running ----------------------------------------------------------------
  void run(std::size_t generations);

  // Introspection ----------------------------------------------------------------
  std::size_t since_better_score() const;
  std::size_t since_better_z() const;
  double best_score() const;
  std::size_t best_z() const;
  Individual best_individual() const;
  std::pair<std::size_t, std::size_t> bound() const;
};

#endif