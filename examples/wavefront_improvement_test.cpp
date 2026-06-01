#include <iostream>
#include <chrono>
#include <nauty/nauty.h>

#include "zero_forcing.hpp"
#include "graph.hpp"

double total_old = 0.0;
double total_new = 0.0;
std::size_t num_tested;

double average_old() {
  return 0.001 * (total_old / num_tested);
}

double average_new() {
  return 0.001 * (total_new / num_tested);
}

bool compare_ver(const Graph &graph) {
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  std::chrono::time_point<std::chrono::high_resolution_clock> end;

  start = std::chrono::high_resolution_clock::now();
  std::size_t old_ver = zero_forcing_wavefront_old(graph);
  end = std::chrono::high_resolution_clock::now();
  total_old += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  std::size_t upper_bound = old_ver; //graph.get_order(); //old_ver; // std::min(graph.get_order(), static_cast<std::size_t>(old_ver * 1.5));

  start = std::chrono::high_resolution_clock::now();
  std::size_t new_ver = zero_forcing_wavefront(graph, upper_bound);
  end = std::chrono::high_resolution_clock::now();
  total_new += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  num_tested++;

  std::cout << old_ver << ", " << new_ver << "\n";
  return old_ver == new_ver;
}

int main() {
  std::size_t path_count = 0; //1000;
  std::size_t cycle_count = 0; //256;
  std::size_t complete_count = 0; //64;
  std::size_t cubic_count = 0; //45;

  std::size_t random_count = 64;
  std::size_t random_subcount = 1; //100;
  double random_min_pct = 0.5;
  double random_max_pct = 0.5;
  double random_delta_pct = 0.1;

  for (std::size_t a = 1; a < path_count; a++) {
    Graph graph = Graph::generate_path(a);
    std::cout << "Path(" << a << "): ";
    if (!compare_ver(graph)) {
      std::cout << "\nFail\n";
      return 1;
    }
  }
  std::cout << "\nAfter Paths (Old, New) = " << average_old() << "ms, " << average_new() << "ms\n\n"; 

  for (std::size_t a = 3; a < cycle_count; a++) {
    Graph graph = Graph::generate_cycle(a);
    std::cout << "Cycle(" << a << "): ";
    if (!compare_ver(graph)) {
      std::cout << "\nFail\n";
      return 1;
    }
  }
  std::cout << "\nAfter Cycles (Old, New) = " << average_old() << "ms, " << average_new() << "ms\n\n"; 

  for (std::size_t a = 3; a < complete_count; a++) {
    Graph graph = Graph::generate_complete(a);
    std::cout << "Complete(" << a << "): ";
    if (!compare_ver(graph)) {
      std::cout << "\nFail\n";
      return 1;
    }
  }
  std::cout << "\nAfter Completes (Old, New) = " << average_old() << "ms, " << average_new() << "ms\n\n"; 

  for (std::size_t a = 4; a < cubic_count; a += 2) {
    Graph graph = Graph::generate_cubic(a);
    std::cout << "Random Cubic(" << a << "): ";
    if (!compare_ver(graph)) {
      std::cout << "\nFail\n";
      return 1;
    }
  }
  std::cout << "\nAfter Cubics (Old, New) = " << average_old() << "ms, " << average_new() << "ms\n\n";
  
  for (std::size_t a = 64; a <= random_count; a++) {
    for (double prob = random_min_pct; prob <= random_max_pct; prob += random_delta_pct) {
      std::vector<Graph> graphs = Graph::generate_random(a, random_subcount, prob);    

      for (const Graph &graph : graphs) {
        std::cout << "Random(" << a << ", " << prob << "): ";
        if (!compare_ver(graph)) {
          std::cout << "\nFail\n";
          return 1;
        }
      }
    } 
  }
  std::cout << "\nAfter Randoms (Old, New) = " << average_old() << "ms, " << average_new() << "ms\n\n"; 

  std::cout << "\nSuccess\n";
}