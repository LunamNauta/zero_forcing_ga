#ifndef GENETIC_SOLVER_HEADER
#define GENETIC_SOLVER_HEADER

#include "weighted_sampler.hpp"
#include "graph.hpp"

#include <cstddef>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

struct GeneticSolverConfig {
  std::size_t population_size = 10;

  // GeneticSolver::select_parents
  std::size_t tournament_size = 3;
  double gamma = 0.8;
  double delta = 0.1;

  // Individual::score
  double alpha = 0.9;
  double beta = 0.1;

  // GeneticSolver::crossover_individual
  std::size_t r = 3;

  // GeneticSolver::crossover_all_individuals
  // Usually keep disabled. Allows calculation of lower bounds
  // Only gets ~10% of the way to Z(G), making it largely a
  // waste of CPU cycles
  bool perform_hitting_set = false;

  bool repair_in_init = true;
  bool repair_in_crossover = true;

  bool trim_in_init = true;
  bool trim_in_crossover = true;

  GeneticSolverConfig(std::size_t psi);
};

struct Individual {
  const Graph *graph;
  const WeightedSampler *sampler;
  const GeneticSolverConfig *_config;
  VertexBitset _initial;
  VertexBitset closure;
  std::size_t pt;
  bool dirty;

  Individual(const Graph *gi, const WeightedSampler *si, const GeneticSolverConfig *ci, const VertexBitset &ii);
  Individual(const Graph *gi, const WeightedSampler *si, const GeneticSolverConfig *ci);

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
  GeneticSolverInfo(const Graph *gi, const WeightedSampler *si, const GeneticSolverConfig *ci);

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
  GeneticSolverConfig _config;
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

  void acknowledge_all_forts();

  void update_info();
  double variance() const;

public:
  GeneticSolver(const Graph *gi, GeneticSolverConfig config);

  void force_upper_bound(std::size_t ub);
  void force_lower_bound(std::size_t lb);

  GeneticSolverConfig config() const;
  GeneticSolverInfo info() const;

  void run(std::size_t generations);
  void kill();
};

#endif
