#ifndef GENETIC_ALGORITHM_HEADER
#define GENETIC_ALGORITHM_HEADER

#include "random_sampler.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"
#include <random>

typedef VertexBitset Gene;

struct Individual {
  const Graph *graph;
  VertexBitset initial;
  VertexBitset closure;
  std::size_t pt;
  bool dirty;

  Individual(const Graph *gi) :
    graph(gi),
    initial(gi->get_order()),
    closure(gi->get_order()),
    pt(0),
    dirty(false)
  {
    initial.set();
    closure.set();
  }
  Individual(const Graph *gi, const VertexBitset &ii) : graph(gi) {
    set_initial(ii);
  }

  void update() {
    closure = initial;
    pt = zero_forcing_closure(*graph, closure);
    dirty = false;
  }

  void set_initial(const VertexBitset &ii) {
    initial = ii;
    initial.resize(graph->get_order(), false);
    dirty = true;
  }
  VertexBitset get_initial() {
    return initial;
  }

  double get_score() {
    if (dirty) update();

    double size_bonus = graph->get_order() - initial.count();
    // double force_penalty = graph->get_order() - closure.count();
    // double pt_bonus = 1.0 - (static_cast<double>(pt) / graph->get_order());
    return size_bonus; // - force_penalty + pt_bonus;
  }
  std::size_t get_z() {
    return initial.count();
  }
  bool forces() {
    if (dirty) update();
    return closure.count() == graph->get_order();
  }

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
  void fix_population();
  void reduce_population();

  // Individual Functions ----------------------------------------------------------------
  Individual crossover_individual(Individual &ind1, Individual &ind2);
  void fix_individual(Individual &ind);
  void reduce_individual(Individual &ind);

  // Selection Functions ----------------------------------------------------------------
  std::pair<Individual&, Individual&> select_parents();

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
  double variance() const;
};

#endif