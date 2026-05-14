#include <iostream>

#include <nauty/naututil.h>
#include <nauty/gtools.h>
#include <nauty/nauty.h>

#include "zero_forcing.hpp"
#include "fort_cover.hpp"
#include "graph.hpp"

std::vector<std::string> generate_random_graph6(int num_graphs, int n) {
  std::vector<std::string> results;
  int m = SETWORDSNEEDED(n);
    
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  char* buffer = nullptr;
  size_t size = 0;

  for (std::size_t a = 0; a < num_graphs; a++) {
    DYNALLSTAT(graph, g, g_sz);
    DYNALLOC2(graph, g, g_sz, m, n, "Allocation failed");
    EMPTYGRAPH(g, m, n);

    for (int v = 0; v < n; ++v) {
      for (int w = v + 1; w < n; ++w) {
        if ((double)rand() / RAND_MAX <= 0.5) continue;
        ADDONEEDGE(g, v, w, m);
      }
    }

    FILE* mem_stream = open_memstream(&buffer, &size);
    if (mem_stream) {
      writeg6(mem_stream, g, m, n);
      fclose(mem_stream);
            
      std::string g6_str(buffer);
      if (!g6_str.empty() && g6_str.back() == '\n') {
        g6_str.pop_back();
      }
      results.push_back(g6_str);
            
      free(buffer);
      buffer = nullptr;
    }

    DYNFREE(g, g_sz);
  }

  return results;
}

int main() {
  Graph graph;
  FortCoverData fc_data;
  std::size_t z;
  
  double submodel_type1_time = 0;
  double submodel_type2_time = 0;
  double submodel_type3_time = 0;
  
  const int num_tests = 1;

  for (int a = 1; a <= num_tests; a++) {
    // Instead of loading from file, generate a ramndom, graph
    std::string graph6 = generate_random_graph6(1, 50)[0];
    graph.from_graph6(graph6);
    
    std::cout << "\nTesting Graph: " << graph6 << "\n";

    // Version 1
    fort_cover_ip(graph, fc_data, 1);
    submodel_type1_time += fc_data.sub_model_time;
    z = zero_forcing_wavefront(graph);
    if (z != fc_data.val) {
      std::cout << "Error: wavefront (" << z << ") != fort cover (" << fc_data.val << ")" << "\n";
      return -1;
    }

    // Version 2
    fort_cover_ip(graph, fc_data, 2);
    submodel_type2_time += fc_data.sub_model_time;
    z = zero_forcing_wavefront(graph);
    if (z != fc_data.val) {
      std::cout << "Error: wavefront (" << z << ") != fort cover (" << fc_data.val << ")" << "\n";
      return -1;
    }

    // Version 3
    fort_cover_ip(graph, fc_data, 3);
    submodel_type3_time += fc_data.sub_model_time;
    z = zero_forcing_wavefront(graph);
    if (z != fc_data.val) {
      std::cout << "Error: wavefront (" << z << ") != fort cover (" << fc_data.val << ")" << "\n";
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