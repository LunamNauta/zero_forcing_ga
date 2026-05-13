#include <iostream>

#include "zero_forcing.hpp"
#include "fort_cover.hpp"

void generate_cubic_prism(Graph& g, int n) {
  g = Graph(2*n);
  // Add 2n vertices and connect them in two rings with cross-edges
  for (std::size_t a = 0; a < n; a++) {
    g.insert_edge(a, (a + 1) % n);         // Inner ring
    g.insert_edge(a + n, (a + 1) % n + n); // Outer ring
    g.insert_edge(a, a + n);               // Spokes connecting rings
  }
}

int main() {
  Graph graph;
  fort_cover_data fc_data;
  std::size_t z;
  
  double submodel_type1_time = 0;
  double submodel_type2_time = 0;
  double submodel_type3_time = 0;
  
  const int num_tests = 5;

  for (int a = 1; a <= num_tests; a++) {
    // Instead of loading from file, generate a cubic graph
    int vertices_per_ring = 4 + a; 
    generate_cubic_prism(graph, vertices_per_ring);
    
    std::cout << "\n--- Testing Cubic Prism Graph (n=" << graph.get_order() << ") ---" << "\n";

    // Version 1
    fort_cover_ip(graph, fc_data, 1);
    submodel_type1_time += fc_data.sub_model_time;
    z = zero_forcing_wavefront(graph);
    if (z != fc_data.val) {
      std::cout << "error: wavefront (" << z << ") != fort cover (" << fc_data.val << ")" << "\n";
      return -1;
    }

    // Version 2
    fort_cover_ip(graph, fc_data, 2);
    submodel_type2_time += fc_data.sub_model_time;
    z = zero_forcing_wavefront(graph);
    if (z != fc_data.val) {
      std::cout << "error: wavefront (" << z << ") != fort cover (" << fc_data.val << ")" << "\n";
      return -1;
    }

    // Version 3
    fort_cover_ip(graph, fc_data, 3);
    submodel_type3_time += fc_data.sub_model_time;
    z = zero_forcing_wavefront(graph);
    if (z != fc_data.val) {
      std::cout << "error: wavefront (" << z << ") != fort cover (" << fc_data.val << ")" << "\n";
      return -1;
    }
  }

  std::cout << "\n========================================" << "\n";
  std::cout << "Average Submodel Times (" << num_tests << " graphs):" << "\n";
  std::cout << "Type 1: " << submodel_type1_time / num_tests << "s" << "\n";
  std::cout << "Type 2: " << submodel_type2_time / num_tests << "s" << "\n";
  std::cout << "Type 3: " << submodel_type3_time / num_tests << "s" << "\n";
  std::cout << "========================================" << "\n";
}