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
  mutation_rate(1)
{
  initialize_population();
}

void GeneticSolver::initialize_population() {
  // Distribution for the approximate number of vertices initially colored in an individual
  std::uniform_real_distribution<double> dist1(0.3, 0.9);

  for (Individual &ind : population) {
    ind.genes.resize(graph->get_order(), false);

    // Choose a random percent of vertices to cover
    double density = dist1(gen);
    std::bernoulli_distribution dist2(density);
    
    // Set approximately @density percent of vertices to colored
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (dist2(gen)) ind.genes.set(u);
    }
  }
}

void GeneticSolver::score_population() {
  // How much to penalize failed zero forcing sets
  const double penalty_weight = 100.0;

  for (Individual &ind : population) {
    // Find the closure of the initial zero forcing set
    VertexBitset forced(ind.genes);
    zero_forcing_closure(*graph, forced);
    acknowledge_fort(~forced);

    // Update whether the individual forces the graph
    ind.forces = forced.count() == graph->get_order();

    // Compute scores
    // Positive points for graphs that cover less initially
    // Negative points based on how much of the graph wasn't forced
    double base_score = graph->get_order() - ind.genes.count(); 
    double unforced_count = graph->get_order() - forced.count();
    double penalty = penalty_weight * unforced_count;
    
    // Compute total score
    ind.score = 1000.0 + base_score - penalty;
  }
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
  const std::size_t tournament_size = 3;

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
  Individual offspring;
  // 70% better than random (0.3%)
  offspring.genes = parent1.genes & parent2.genes;
  return offspring;

  // Random (1%)
  /*
  offspring.genes.resize(graph->get_order());
  std::bernoulli_distribution dist(0.5);
  for (Vertex u = 0; u < graph->get_order(); u++) {
    offspring.genes[u] = dist(gen);
  }
  */

  // 80% better than random (0.2%)
  /*
  offspring = parent2;
    
  std::uniform_int_distribution<Vertex> v_dist(0, graph->get_order() - 1);
  Vertex seed = v_dist(gen);

  std::queue<Vertex> queue;
  queue.push(seed);
  std::vector<bool> visited(graph->get_order(), false);
  visited[seed] = true;
    
  std::size_t count = 0;
  std::size_t limit = graph->get_order() * 0.7;

  while (!queue.empty() && count < limit) {
    Vertex curr = queue.front();
    queue.pop();
        
    offspring.genes[curr] = parent1.genes[curr];
    count++;

    for (Vertex neighbor : graph->get_adjacent(curr)) {
      if (visited[neighbor]) continue;
      visited[neighbor] = true;
    }
  }
  */

  // 75% better than random (0.25%)
  /*
  offspring.genes = parent1.genes | parent2.genes;
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (parent1.genes[u] != parent2.genes[u]) {
      if (dist(gen) >= 0.5) continue;
      offspring.genes[u] = false;
    }
  }
  */

  // 70% better than random (0.3%)
  /*
  offspring.genes.resize(graph->get_order(), false);
  std::bernoulli_distribution dist(0.5);
  
  for (Vertex u = 0; u < graph->get_order(); u++) {
    offspring.genes[u] = dist(gen) ? parent1.genes[u] : parent2.genes[u];
  }
  */
}

void GeneticSolver::mutate_individual(Individual &ind) {
  // Distribution for whether a vertex is flipped
  std::bernoulli_distribution dist(mutation_rate);

  // Invert approximately @mutation_rate percent of the individuals vertices
  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (!dist(gen)) continue;
    ind.genes.flip();
  }
}

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

    // Add 1 vertex from each of the known forts that wasn't covered
    for (std::unordered_set<VertexBitset>::const_iterator it_fort = minimal_forts.cbegin(); it_fort != minimal_forts.cend(); it_fort++) {
      if ((*it_fort & forced).any()) continue; // Ignore this fort if we've already covered it
      VertexBitset verts = sampler.sample_bitset(1, ~*it_fort);
      ind.genes |= verts;
      forced |= verts;
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
  VertexBitset minimal_fort = extract_smaller_fort(fort);
  if (minimal_fort.any()) {
    sampler.update_weights(minimal_fort);
    minimal_forts.insert(minimal_fort);
  } 
}

VertexBitset GeneticSolver::extract_smaller_fort(VertexBitset fort) {
  // Buffer to store the vertices that *could* be in a smaller fort
  std::vector<Vertex> candidates;
  candidates.reserve(fort.count()); 

  // Buffer to store the "fort degree" of each vertex
  // This is the number of vertices in this fort that each vertex is adjacent to
  std::vector<std::size_t> fort_degree(graph->get_order(), 0);
  for (Vertex u = 0; u < graph->get_order(); u++) {
    if (!fort[u]) continue;
    candidates.push_back(u);
    for (Vertex v : graph->get_adjacent(u)) fort_degree[v]++;
  }
  
  // Shuffle so we discover a random minimal sub-fort each time
  std::shuffle(candidates.begin(), candidates.end(), gen);
  
  // Attempt to remove each vertex from the fort
  for (Vertex u : candidates) {
    // If the fort degree of @u is 1, removing it would violate the fort property
    if (fort_degree[u] == 1) continue;
    
    // For each vertex not in the fort, will removing @u violate the fort property
    bool can_remove = true;
    for (Vertex v : graph->get_adjacent(u)) {
      // If we remove @u, @v would violate the fort property (@v could force the fort)
      if (fort[v] || fort_degree[v] != 2) continue;
      can_remove = false;
      break;
    }
    if (!can_remove) continue;
    
    // Remove @u from the fort, and update the fort degrees
    fort.reset(u);
    for (Vertex v : graph->get_adjacent(u)) fort_degree[v]--;
  }
  
  // @fort is now (hopefully) a smaller fort
  return fort;
}

Individual GeneticSolver::run(std::size_t num_generations) {
  // Score the current population (either from initialization or from a previous run)
  score_population();

  for (std::size_t a = 0; a < num_generations; a++) {
    // Crossover + mutate + re-score
    crossover_population();
    mutate_population();
    score_population();
    
    // Update the time since an improvement was made
    if (population[0].genes.count() < best_score) {
      best_score = population[0].genes.count();
      since_best = 0;
    }
    else since_best++;

    // If we haven't seen improvement, increase the mutation rate (up to 20%)
    if (since_best > 5) mutation_rate += 0.005;
    else mutation_rate -= 0.005;
    mutation_rate = std::max(0.0, std::min(mutation_rate, 0.2));
  }

  // Fix the population to ensure all individuals are zero forcing sets
  fix_population();

  // Re-score the population. The last step may have found better solutions
  score_population();

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