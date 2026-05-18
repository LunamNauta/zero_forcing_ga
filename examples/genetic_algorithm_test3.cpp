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

std::size_t test_graphs(std::size_t order, std::size_t num_graphs, double edge_prob) {
  std::vector<Graph> graphs = Graph::generate_random(order, num_graphs, edge_prob);

  double total_error = 0;
  double num_failed = 0;
  double total_generations = 0;
  std::size_t generation = 0;
  std::size_t best_generation = 0;

  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  std::chrono::time_point<std::chrono::high_resolution_clock> end;

  for (const Graph &graph : graphs) {
    if (!graph.is_connected()) continue;
    // std::cout << "Graph: \n" << graph << "\n";

    start = std::chrono::high_resolution_clock::now();

    GeneticSolver solver(&graph, order);
    VertexSet result;
 
    generation = 0;
    while (true) {
      result = std::move(solver.run_set(1)); //std::sqrt(graph.get_order())));
      generation++;
      //if (solver.time_since_better_z() == 0) std::cout << "Gen " << generation << ": Better Z(G) -> " << solver.current_z() << "\n";
      //else if (solver.time_since_better_score() == 0) std::cout << "Gen " << generation << ": Better Score -> " << solver.current_score() << "\n";
      //else std::cout << "Gen " << generation << "\n";
      if (solver.time_since_better_score() > 8) break;
    }

    end = std::chrono::high_resolution_clock::now();
    double duration = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000 / 1000;

    //std::cout << "Observed Z(G): " << result.size() << "\n";
    //std::cout << "Compute Time: " << duration << " seconds" << "\n";

    std::size_t expected = zero_forcing_wavefront(graph);
    std::int64_t delta = result.size() - expected;

    //std::cout << "Real Z(G): " << expected << "\n";
    //std::cout << "Delta: " << delta << "\n";
    if (delta < 0) {
      //std::cerr << "Crital Error: Zero forcing closure logic failed" << "\n";
      return 1;
    }

    if (result.size() != expected) std::cout << "Attempting to reach expected answer" << "\n";
    while (result.size() != expected) {
      result = std::move(solver.run_set(1));
      generation++;
      //if (solver.time_since_better_z() == 0) std::cout << "Gen " << generation << ": Better Z(G) -> " << solver.current_z() << "\n";
      //else if (solver.time_since_better_score() == 0) std::cout << "Gen " << generation << ": Better Score -> " << solver.current_score() << "\n";
      //else std::cout << "Gen " << generation << "\n";
    }
    //std::cout << "Expected answer took " << generation << " generation(s)" << "\n";

    // std::cout << "\n";

    total_error += static_cast<double>(delta) / order;
    if (delta != 0) num_failed++;
    total_generations += best_generation;
  }

  double average_accuracy = 1.0 - total_error / num_graphs;
  double average_error = total_error / num_graphs;
  double fail_rate = num_failed / num_graphs;
  double average_generations = total_generations / num_graphs;

  std::cout << "Order: " << order << "\n";
  std::cout << "Number of Graphs: " << num_graphs << "\n";
  std::cout << "Edge Probability: " << edge_prob << "\n"; 
  std::cout << "Average Accuracy: " << average_accuracy*100 << "%" << "\n";
  std::cout << "Average Error: " << average_error*100 << "%" << "\n";
  std::cout << "Fail Rate: " << fail_rate*100 << "%" << "\n";
  std::cout << "Average #Generations: " << average_generations << "\n";
  std::cout << "\n";

  return num_failed;
}

int main(int argc, char* argv[]) {
  for (std::size_t a = 3; a < 32; a++) {
    for (double prob = 0.3; prob < 0.9; prob += 0.1) {
      if (test_graphs(a, a, prob)) return 1;
    }
  }
}