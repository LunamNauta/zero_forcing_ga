#ifndef GENETIC_ALGORITHM_HEADER
#define GENETIC_ALGORITHM_HEADER

#include <random>

#include "random_sampler.hpp"
#include "graph.hpp"

/** @brief Alias for the genetic representation of a solution. */
typedef VertexBitset Chromosome;

/** @brief The (maximum) percentage chance of a vertex being added/removed during mutation. */
inline constexpr double GA_MUTATION_PCT = 0.5;

/** @brief The (maximum) percentage of the population preserved as "elites" during crossover and mutation. */
inline constexpr double GA_ELITE_PCT = 0.33;

/** @brief Represents a single solution in the genetic population. */
struct Individual {
  Chromosome genes; //!< The bitset representing the potential solution
  double score;     //!< Fitness score calculated by score_population()
  bool forces;      //!< Cache flag: true if this set forces the entire graph
};

/**
 * @class GeneticSolver
 * @brief A heuristic solver using a genetic algorithm to compute Z(G).
 * 
 * This class implements a genetic algorithm to optimize individuals for minimal
 * Z(G) values.
 */
class GeneticSolver {
private:
  const Graph *graph;                             //!< The graph being analyzed for Z(G) and forts.
  std::vector<Individual> population;             //!< The current generation of candidate solutions.
  RandomSampler sampler;                          //!< Utility for weighted random selection of vertices.
  std::mt19937 gen;                               //!< The Mersenne Twister engine used for all stochastic operations.
  double mutation_rate;                           //!< The probability (0.0 to 1.0) of a mutation occurring during crossover.
  double elite_pct;                               //!< The percent of the population that are elites (excluded from operations such as mutation).
  std::unordered_set<VertexBitset> known_forts;   //!< A collection of unique, small forts discovered across all generations.
  std::size_t since_better_score;                 //!< The number of generations since best_score changed.
  std::size_t since_better_z;                     //!< The number of generations since best_z changed.
  double best_score;                              //!< The highest score found so far
  std::size_t best_z;                             //!< The lowest Z(G) found so far

  /**
  * @brief Executes the genetic evolution loop.
  * @param num_generations The number of generations to execute.
  * @return The best individual found during the execution.
  */
  Individual run(std::size_t num_generations);

  /** @brief Populates the initial generation with candidate solutions. */
  void initialize_population();

  /**
   * @brief Evaluates the fitness of every individual in the current population.
   * 
   * Higher scores are assigned to individuals with smaller zero forcing sets,
   * and failed zero forcing sets are heavily penalized
   */
  void score_population();

  /**
   * @brief Updates the population with crossovers (ignoring elites).
   * 
   * Retains the top performing individuals and replaces the lower performing 
   * half with offspring generated from 'crossover_individual()'.
   */
  void crossover_population();

  /**
   * @brief Applies random mutations across the population.
   * 
   * This introduces genetic diversity to prevent the solver from stalling 
   * in local optima.
   */
  void mutate_population();

  /**
   * @brief Repairs invalid candidate solutions.
   * 
   * Ensures all individuals represent valid zero forcing sets
   */
  void fix_population();

  /**
  * @brief Selects a high-fitness individual to act as a parent.
  * @return A reference to the chosen Individual.
  */
  Individual& select_parent();

  void score_individual(Individual &ind);

  /**
   * @brief Combines two parents to produce a new offspring.
   * 
   * Crossover focuses on merging structural traits from successful solutions.
   *
   * @param parent1 The first genetic donor.
   * @param parent2 The second genetic donor.
   * @return A new Individual representing the child.
   */
  Individual crossover_individual(const Individual &parent1, const Individual &parent2);

 /**
  * @brief Randomly toggles vertices within an individual's chromosome.
  * @param ind The individual to modify.
  */
  void mutate_individual(Individual &ind);

 /**
  * @brief Adds vertices to a set until it satisfies the forcing property.
  * @param ind The individual to be repaired.
  */
  void fix_individual(Individual &ind);

  /** @brief Converts internal chromosome format to a VertexBitset. */
  static VertexBitset chromosome_to_bitset(const Chromosome &genes);

  /** @brief Converts internal chromosome format to a VertexSet. */
  static VertexSet chromosome_to_set(const Chromosome &genes);

 /**
  * @brief Processes and stores a new fort.
  * 
  * The fort is reduced to a smaller form via 'extract_smaller_fort'.
  * before being stored in 'known_forts'.
  * @param fort The newly discovered fort.
  */
  void acknowledge_fort(const VertexBitset &fort);

  /**
   * @brief Greedily minimizes a fort while preserving its structural integrity.
   * 
   * Attempts to remove vertices one by one, ensuring that the remaining set
   * still satisfies the "no vertex not in the fort has exactly 1 neighbor in the fort" property.
   * @param fort The fort to minimize.
   * @return A subset of the original fort that is smaller.
   */
  VertexBitset extract_smaller_fort(VertexBitset fort);

public:
  /**
   * @brief Constructs the solver for a specific graph.
   * @param gi Pointer to the graph to be analyzed.
   * @param population_size The number of individuals to maintain in each generation.
   */
  GeneticSolver(const Graph *gi, std::size_t population_size);

  /**
   * @brief Executes the solver and returns the result as a VertexBitset.
   * @param num_generations The number of generations to execute.
   * @return The best VertexBitset found.
   */
  VertexBitset run_bitset(std::size_t num_generations);

  /**
   * @brief Executes the solver and returns the result as a VertexSet.
   * @param num_generations The number of generations to execute.
   * @return The best VertexSet found.
   */
  VertexSet run_set(std::size_t num_generations);

  /**
   * @brief Tracks stagnation in the evolution process.
   * @return The number of generations elapsed since the last improvement in score.
   */
  std::size_t time_since_better_score() const;

  /**
   * @brief Tracks stagnation in the evolution process.
   * @return The number of generations elapsed since the last improvement in Z(G).
   */
  std::size_t time_since_better_z() const;

  /**
   * @brief Tracks the current best score
   * @return The highest score seen so far
   */
  double current_score() const;

  /**
   * @brief Tracks the current best Z(G)
   * @return The lowest Z(G) seen so far
   */
  std::size_t current_z() const;
};

#endif