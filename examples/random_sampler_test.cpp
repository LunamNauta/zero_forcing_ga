#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include "random_sampler.hpp"
#include "graph.hpp"

int main() {
  std::random_device rd;
  std::mt19937 local_gen(rd());
    
  std::uniform_int_distribution<std::size_t> order_dist(10, 60);
  std::size_t random_order = order_dist(local_gen);
    
  Graph graph;
  std::string g6_header = "";
  g6_header += (char)(random_order + 63);
  
  std::size_t num_bits = (random_order*(random_order - 1))/2;
  std::size_t num_chars = (num_bits + 5) / 6;
  g6_header.append(num_chars, '?');

  graph.from_graph6(g6_header);

  std::cout << "Testing Random Graph Scale: " << graph.get_order() << " vertices\n";
  RandomSampler sampler(&graph);

  std::uniform_int_distribution<std::size_t> weight_dist(1, 50);
  std::vector<std::pair<VertexIndex, std::size_t>> random_updates;
    
  for (std::size_t a = 0; a < graph.get_order(); a++) {
    random_updates.emplace_back(a, weight_dist(local_gen));
  }
  sampler.update_weights(random_updates);

  std::vector<std::size_t> counts(graph.get_order(), 0);
  std::size_t iterations = 30000;
  for (std::size_t a = 0; a < iterations; a++) {
    Graph sub = sampler(1, 10); 
    if (sub.get_order() == 0) continue;
    counts[sub.get_label(0)]++;
  }

  int failures = 0;
    
  std::cout << "\n" << std::setw(6) << "Vertex" 
            << " | " << std::setw(12) << "Target Prob" 
            << " | " << std::setw(13) << "Observed Prob" 
            << " | " << "Status" << "\n";
  std::cout << std::string(50, '-') << "\n";

  for (std::size_t a = 0; a < graph.get_order(); a++) {
    double target = sampler.get_weight(a);
    double observed = static_cast<double>(counts[a]) / iterations;
    double deviation = std::abs(target - observed);
        
    bool passed = (deviation <= 0.015); // 1.5% tolerance
    if (!passed) failures++;

    std::cout << std::setw(6) << a << " | " 
              << std::setw(12) << std::fixed << std::setprecision(4) << target << " | " 
              << std::setw(13) << observed << " | " 
              << (passed ? "OK" : "FAIL") << "\n";
  }

  std::cout << std::string(50, '-') << "\n";
  if (failures == 0) std::cout << "OVERALL RESULT: SUCCESS (All vertices within tolerance)\n";
  else std::cout << "OVERALL RESULT: FAIL (" << failures << " vertices outside tolerance)\n";
}