#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

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

  std::uniform_int_distribution<std::size_t> freq_dist(1, 50);
  std::vector<std::size_t> fort_appearances(graph.get_order(), 0);
    
  for (std::size_t a = 0; a < graph.get_order(); a++) {
    std::size_t frequency = freq_dist(local_gen);
    fort_appearances[a] = frequency;
    for (std::size_t f = 0; f < frequency; f++) {
      VertexSet fort = { static_cast<Vertex>(a) };
      sampler.update_weights(fort);
    }
  }

  std::vector<std::size_t> counts(graph.get_order(), 0);
  std::size_t iterations = 30000;
  for (std::size_t a = 0; a < iterations; a++) {
    Graph sub = graph.subgraph(sampler.sample_set(1, 10)); 
    if (sub.get_order() == 0) continue;
    counts[sub.get_label(0)]++;
  }

  std::cout << "\n" << std::setw(6) << "Vertex" 
            << " | " << std::setw(12) << "Fort Count" 
            << " | " << std::setw(15) << "Observed %" << "\n";
  std::cout << std::string(40, '-') << "\n";

  for (std::size_t a = 0; a < graph.get_order(); a++) {
    double observed_pct = (static_cast<double>(counts[a]) / iterations) * 100.0;

    std::cout << std::setw(6) << a << " | " 
              << std::setw(12) << fort_appearances[a] << " | " 
              << std::setw(14) << std::fixed << std::setprecision(2) << observed_pct << "%" << "\n";
  }

  std::cout << std::string(40, '-') << "\n";
  std::cout << "Test completed over " << iterations << " iterations.\n";
}