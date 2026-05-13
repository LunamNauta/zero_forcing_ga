#include <iostream>
#include <cassert>

#include "genetic_algorithm.hpp"
#include "graph.hpp"

#include <nauty/naututil.h>
#include <nauty/gtools.h>
#include <nauty/nauty.h>

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

// 10 -> 0.25
// 20 -> 0.5
// 40 -> 2.5
// 80 -> 23

int main() {
  std::vector<std::string> graph_strs = generate_random_graph6(100, 40);
  Graph graph;

  for (const std::string &graph6 : graph_strs) {
    graph.from_graph6(graph6);

    ZFGeneticSolver solver(&graph, 100);
    VertexSet result = solver.run(100);

    std::cout << "Graph: " << graph6 << "\n";
    std::cout << "Z(G): " << result.size() << "\n";
  }

}