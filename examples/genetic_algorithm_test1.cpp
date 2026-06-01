#include <chrono>
#include <iostream>
#include <cassert>

#include <nauty/naututil.h>
#include <nauty/gtools.h>
#include <nauty/nauty.h>

#include "genetic_algorithm.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"

/*
// GeneticSolver::run -> 
//  O(G*P*(V*E))
//  If G ~ P ~ V ~ E -> O(x^4) on average
*/

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <number_of_graphs> <population_size> <generations_before_quit> <graph_order>" << "\n";
    return 1; 
  }

  std::size_t num_graphs;
  std::size_t population_size;
  std::size_t generations_before_quit;
  std::size_t order;

  num_graphs = std::stoull(argv[1]);
  population_size = std::stoull(argv[2]);
  generations_before_quit = std::stoull(argv[3]);
  order = std::stoull(argv[4]);

  std::vector<Graph> graphs = Graph::generate_random(order, num_graphs, 0.6);

  double total_error = 0;
  double num_failed = 0;
  double total_generations = 0;
  std::size_t generation = 0;
  std::size_t best_generation = 0;
  std::size_t generation_found = 0;
  std::size_t max_generation = 0;

  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  std::chrono::time_point<std::chrono::high_resolution_clock> end;

  for (const Graph &graph : graphs) {
    if (!graph.is_connected()) continue;

    start = std::chrono::high_resolution_clock::now();

    GeneticSolver solver(&graph, population_size);
    
    double best_variance = solver.variance();
 
    std::cout << "Gen 0: Better Z(G) -> " << solver.best_z() << "\n";
    std::cout << "Gen 0: Initial Variance -> " << best_variance << "\n";

    generation = 0;
    generation_found = 0;
    while (true) {
      solver.run(1);
      generation++;
      
      double current_variance = solver.variance();
      if (current_variance < best_variance) {
        std::cout << "Gen " << generation << ": Better Variance -> " << current_variance << "\n";
        best_variance = current_variance;
      }

      if (solver.since_better_z() == 0) {
        std::cout << "Gen " << generation << ": Better Z(G) -> " << solver.best_z() << "\n";
        generation_found = generation;
        max_generation = std::max(max_generation, generation);
      }
      else if (solver.since_better_score() == 0) std::cout << "Gen " << generation << ": Better Score -> " << solver.best_score() << "\n";
      // else std::cout << "Gen " << generation << "\n";
      
      if (current_variance < 0.1 || solver.since_better_score() > generations_before_quit) {
        std::cout << "Gen " << generation << ": Terminating early. Variance: " << current_variance 
                  << ", Gens since better score: " << solver.since_better_score() << "\n";
        break;
      }
    }

    end = std::chrono::high_resolution_clock::now();
    double duration = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000 / 1000;

    std::cout << "Observed Z(G): " << solver.best_z() << "\n";
    std::cout << "Compute Time: " << duration << " seconds" << "\n";

    std::size_t expected = zero_forcing_wavefront(graph);
    std::int64_t delta = solver.best_z() - expected;

    std::cout << "Real Z(G): " << expected << "\n";
    std::cout << "Delta: " << delta << "\n";
    if (delta < 0) {
      std::cerr << "Crital Error: Zero forcing closure logic failed" << "\n";
      return 1;
    }

    if (solver.best_z() != expected) std::cout << "Attempting to reach expected answer" << "\n";
    while (solver.best_z() != expected) {
      solver.run(1);
      generation++;
      
      double current_variance = solver.variance();
      if (current_variance < best_variance) {
        std::cout << "Gen " << generation << ": Better Variance -> " << current_variance << "\n";
        best_variance = current_variance;
      }

      if (solver.since_better_z() == 0) {
        std::cout << "Gen " << generation << ": Better Z(G) -> " << solver.best_z() << "\n";
        generation_found = generation;
        max_generation = std::max(max_generation, generation);
      }
      else if (solver.since_better_score() == 0) std::cout << "Gen " << generation << ": Better Score -> " << solver.since_better_score() << "\n";
      else std::cout << "Gen " << generation << "\n";
    }
    std::cout << "Expected answer took " << generation_found << " generation(s)" << "\n";

    std::cout << "\n";

    total_error += static_cast<double>(delta) / order;
    if (delta != 0) num_failed++;
    total_generations += generation_found;
  }

  double average_accuracy = 1.0 - total_error / num_graphs;
  double average_error = total_error / num_graphs;
  double fail_rate = num_failed / num_graphs;
  double average_generations = total_generations / num_graphs;

  std::cout << "Average Accuracy: " << average_accuracy*100 << "%" << "\n";
  std::cout << "Average Error: " << average_error*100 << "%" << "\n";
  std::cout << "Fail Rate: " << fail_rate*100 << "%" << "\n";
  std::cout << "Average #Generations: " << average_generations << "\n";
  std::cout << "Maximum #Generations: " << max_generation << "\n";
}

/* Population size 10: (Order -> #Generations)
20 -> 10
32 -> 20
64 -> 250
*/