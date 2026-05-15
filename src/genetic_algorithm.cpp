#include "genetic_algorithm.hpp"

#include <limits>
#include <random>

#include "random_sampler.hpp"
#include "zero_forcing.hpp"
#include "graph.hpp"

GeneticSolver::GeneticSolver(const Graph *gi, std::size_t population_size) :
  graph(gi),
  population(population_size),
  sampler(gi),
  gen(std::random_device{}()),
  best_score(std::numeric_limits<std::size_t>::max()),
  since_best(0),
  mutation_rate(GA_MUTATION_PCT),
  elite_pct(GA_ELITE_PCT)
{
  initialize_population();
}

void GeneticSolver::initialize_population() {
  double step = static_cast<double>(graph->get_order()) / population.size();
  double order = 1;

  for (std::size_t a = 0; a < population.size(); a++, order += step) {
    VertexBitset initial = sampler.sample_bitset(order);
    population[a].genes = initial;
  }
}

void GeneticSolver::score_individual(Individual &ind) {
  // How much to penalize failed zero forcing sets
  const double penalty_weight = 1;
  // How much to reward larger closures
  const double bonus_weight = 1.5;

  // Find the closure of the initial zero forcing set
  VertexBitset forced(ind.genes);
  if (!ind.forces) {
    VertexBitset forced(ind.genes);
    zero_forcing_closure(*graph, forced);
    acknowledge_fort(~forced);
  }
  else forced.set();

  // Update whether the individual forces the graph
  ind.forces = forced.count() == graph->get_order();

  // Compute scores
  // Positive points for graphs that cover less initially
  // Negative points based on how much of the graph wasn't forced
  double bonus_score = graph->get_order() - ind.genes.count(); 
  double penalty_score = graph->get_order() - forced.count();
    
  // Compute total score
  ind.score = (bonus_score*bonus_weight) / (1.0 + penalty_score*penalty_weight);
}

void GeneticSolver::score_population() {
  for (Individual &ind : population) score_individual(ind);

  std::sort(population.begin(), population.end(), [](const Individual &ind1, const Individual &ind2){
    return ind1.score > ind2.score;
  });
}

void GeneticSolver::crossover_population() {
  // Sort the population based on their score. Better scores come first
  std::sort(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
    return a.score > b.score;
  });

  // Declare the next generation
  std::vector<Individual> next;
  next.reserve(population.size());

  // Preserve the elite members of the population
  std::size_t elites = population.size() * GA_ELITE_PCT;
  for (std::size_t a = 0; a < elites; a++) {
    next.push_back(population[a]);
  }

  while (next.size() < population.size()) {
    // Select two parents (that are not the same individual)
    Individual &parent1 = select_parent();
    Individual &parent2 = select_parent();
    if (&parent1 == &parent2) continue;

    // Add the child to the next generation
    next.emplace_back(crossover_individual(parent1, parent2));
  }

  // Replace the old population with the new one
  population = std::move(next);
}

void GeneticSolver::mutate_population() {
  // Exclude the elite members of the population from mutation
  std::size_t elites = population.size() * GA_ELITE_PCT;

  // Mutate the rest of the population
  for (std::size_t a = elites; a < population.size(); a++) {
    mutate_individual(population[a]);
  }
}

void GeneticSolver::fix_population() {
  // Exclude the elite members of the population from being fixed
  std::size_t elites = population.size() * GA_ELITE_PCT;

  // Fix the rest of the population
  for (std::size_t a = elites; a < population.size(); a++) {
    fix_individual(population[a]);
  }
}

Individual& GeneticSolver::select_parent() {
  // The number of candidates to choose from
  const std::size_t tournament_size = std::sqrt(population.size());

  // Distribution to select the next parent
  std::uniform_int_distribution<std::size_t> dist(0, population.size() - 1);

  // Choose the highest scoring individual out of all candidates
  Individual &best = population[dist(gen)];
  for (std::size_t a = 1; a < tournament_size; a++) {
    Individual &candidate = population[dist(gen)];
    if (candidate.score <= best.score) continue;
    best = candidate;
  }

  return best;
}

Individual GeneticSolver::crossover_individual(const Individual &parent1, const Individual &parent2) {
  Individual best_candidate = parent1;
  
  std::size_t min_order = std::min(parent1.genes.count(), parent2.genes.count());
  std::bernoulli_distribution dist(0.5);

  for (std::size_t a = 1; a < min_order; a++) {
    VertexBitset subset = parent1.genes | parent2.genes;
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (!dist(gen)) continue;
      subset.set(u);
    }
    subset = sampler.sample_bitset(a, ~subset);

    Individual tmp_ind;
    tmp_ind.genes = std::move(subset);
    score_individual(tmp_ind);

    if (tmp_ind.score > best_candidate.score) best_candidate = tmp_ind;
  }

  return best_candidate;
}

void GeneticSolver::mutate_individual(Individual &ind) {
  std::bernoulli_distribution dist(mutation_rate);

  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (!dist(gen)) continue;
    ind.genes.flip(u);
  }
  ind.forces = false;
}

// TODO: Make this work for disconnected graphs
void GeneticSolver::fix_individual(Individual &ind) {
  // If we already know that this individual forces the graph, we can ignore it
  if (ind.forces) return;

  // Force the initial set of vertices
  VertexBitset forced(ind.genes);

  // Modify the individual until it forces the entire graph
  while (forced.count() < graph->get_order()) {
    // Find the closure of the remaining graph
    zero_forcing_closure(*graph, forced);
    if (forced.count() == graph->get_order()) break;
    acknowledge_fort(~forced);

    // Find vertices with exactly two uncolored neighbors
    // Coloring one of the uncolored neighbors then colors the second automatically
    bool found_clog = false;
    for (Vertex u = 0; u < graph->get_order(); ++u) {
      if (!forced[u]) continue;
      std::size_t unforced = 0;
      bool on1 = true;
      Vertex uncolored1 = 0;
      Vertex uncolored2 = 0;
      for (Vertex v : graph->get_adjacent(u)) {
        if (forced[v]) continue;
        if (on1) uncolored1 = v;
        else uncolored2 = v;
        on1 = !on1;
        unforced++;
      }
      if (unforced == 2) {
        ind.genes.set(uncolored1);
        forced.set(uncolored1);
        forced.set(uncolored2);
        found_clog = true;
        break;
      }
    }
    if (found_clog) continue;

    // Find all vertices that are directly preventing zero forcing
    VertexBitset edge(graph->get_order(), false);
    if (forced.count() != 0) {
      for (Vertex u = 0; u < graph->get_order(); u++){
        if (!forced[u]) continue;
        for (Vertex v : graph->get_adjacent(u)) {
          if (forced[v]) continue;
          edge.set(v);
        }
      }
    }
    else edge.set();

    // Add 1 vertex from each of the known forts that wasn't covered
    // Prioritize vertices that are both in the fort and preventing zero forcing
    for (std::unordered_set<VertexBitset>::const_iterator it_fort = known_forts.cbegin(); it_fort != known_forts.cend(); it_fort++) {
      if ((*it_fort & forced).any()) continue; // Ignore this fort if we've already covered it
      VertexBitset verts; // = sampler.sample_bitset(1, ~*it_fort);
      if ((*it_fort & edge).any()) verts = std::move(sampler.sample_bitset(1, ~(*it_fort & edge)));
      else verts = std::move(sampler.sample_bitset(1, ~*it_fort));
      ind.genes |= verts;
      forced |= verts;
    }
  }

  // Attempt to remove vertices
  for (Vertex v = 0; v < graph->get_order(); ++v) {
    if (ind.genes[v]) {
      ind.genes.reset(v);
      VertexBitset forced(ind.genes);
      zero_forcing_closure(*graph, forced);
      if (forced.count() < graph->get_order()) ind.genes.set(v);
    }
  }

  // We know for a fact that this individual forces the graph
  ind.forces = true;
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

void GeneticSolver::acknowledge_fort(const VertexBitset &fort) {
  // Add the fort to the list of known forts
  // Extract a smaller fort from it first
  VertexBitset smaller_fort = extract_smaller_fort(fort);
  if (smaller_fort.any()) {
    sampler.update_weights(smaller_fort);
    known_forts.insert(smaller_fort);
  } 
}

VertexBitset GeneticSolver::extract_smaller_fort(VertexBitset fort) {
  VertexBitset closure = ~fort;

  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (closure.test(u)) continue;

    VertexBitset forced(closure); 
    forced.set(u);
    
    zero_forcing_closure(*graph, forced);
    
    if (forced.count() != graph->get_order()) closure = forced;
  }

  return ~closure;
}

Individual GeneticSolver::run(std::size_t num_generations) {
  // Score the current population (either from initialization or from a previous run)
  score_population();
  std::size_t since_better = 0;
  double best_score_raw = 0;

  for (std::size_t a = 0; a < num_generations; a++) {
    // Crossover + mutate + re-score
    crossover_population();
    mutate_population();
    score_population();
    fix_population();
    score_population();
    
    // Update the time since an improvement was made
    if (population[0].genes.count() < best_score) {
      best_score = population[0].genes.count();
      since_best = 0;
    }
    else since_best++;

    if (population[0].score > best_score_raw) {
      best_score_raw = population[0].score;
      since_better = 0;
    }
    else since_better++;

    // If we haven't seen improvement, increase the mutation rate, and decrease the number of elites
    if (since_better > population.size()) {
      mutation_rate *= 1.01;
      elite_pct *= 0.99;
    }
    else {
      mutation_rate *= 0.99;
      elite_pct *= 1.01;
    }

    mutation_rate = std::max(0.0, std::min(mutation_rate, GA_MUTATION_PCT));
    elite_pct = std::max(0.0, std::min(elite_pct, GA_ELITE_PCT));
  }

  // Return the most optimal specimen
  return population[0];
}

VertexBitset GeneticSolver::run_bitset(std::size_t num_generations) {
  return run(num_generations).genes;
}

VertexSet GeneticSolver::run_set(std::size_t num_generations) {
  return chromosome_to_set(run(num_generations).genes);
}

std::size_t GeneticSolver::time_since_best() const {
  return since_best;
}