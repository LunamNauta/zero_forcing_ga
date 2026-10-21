#ifndef GENETIC_SOLVER_HEADER
#define GENETIC_SOLVER_HEADER

#include "../include/zero_forcing.hpp"
#include "../include/graph.hpp"
#include "weighted_sampler.hpp"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

using Gene = VertexBitset;

struct Individual {
  const Graph *graph;
  const WeightedSampler *sampler;
  VertexBitset _initial;
  VertexBitset closure;
  std::size_t pt;
  bool dirty;

  Individual(const Graph *gi, const WeightedSampler *si, const VertexBitset &ii);
  Individual(const Graph *gi, const WeightedSampler *si);

  void update();

  VertexBitset initial(const VertexBitset &ii);
  VertexBitset initial() const;

  double score();
  std::size_t z();
  bool forces();
};

class GeneticSolverInfo {
private:
  Individual _best_individual;
  double _best_variance;
  double _best_score;
  std::size_t _best_z;
  double _variance;

  std::size_t _since_better_variance;
  std::size_t _since_better_score;
  std::size_t _since_better_z;

  std::size_t _upper_bound;
  std::size_t _lower_bound;

public:
  GeneticSolverInfo(const Graph *gi, const WeightedSampler *si);

  Individual best_individual() const;
  double best_variance() const;
  double best_score() const;
  std::size_t best_z() const;
  double variance() const;

  std::size_t since_better_variance() const;
  std::size_t since_better_score() const;
  std::size_t since_better_z() const;

  std::pair<std::size_t, std::size_t> bound() const;
  std::size_t upper_bound() const;
  std::size_t lower_bound() const;

  friend class GeneticSolver;
};

class GeneticSolver {
private:
  const Graph *graph;
  WeightedSampler sampler;
  GeneticSolverInfo _info;
  std::vector<Individual> population;
  std::mt19937 gen;
  std::unordered_set<VertexBitset> known_forts;
  bool dead;

  void initialize_population();

  std::pair<Individual&, Individual&> select_parents();
  Individual crossover_individual(Individual &ind1, Individual &ind2);
  void crossover_all_individuals();

  void repair_individual(Individual &ind);
  void repair_all_individuals();

  void trim_individual(Individual &ind);
  void trim_all_individuals();

  void trim_fort(VertexBitset &fort);
  void trim_all_forts();

  void acknowledge_fort(const VertexBitset &fort);
  void acknowledge_all_forts();

  void update_info();
  double variance() const;

public:
  GeneticSolver(const Graph *gi, std::size_t psi);

  void force_upper_bound(std::size_t ub);
  void force_lower_bound(std::size_t lb);

  GeneticSolverInfo info() const;
  void run(std::size_t generations);
  void kill();
};

#endif
