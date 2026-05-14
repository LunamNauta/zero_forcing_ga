#include <chrono>
#include <iostream>
#include <cassert>

#include <nauty/naututil.h>
#include <nauty/gtools.h>
#include <nauty/nauty.h>

#include "genetic_algorithm.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"

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

/*
4, 4.6e-05
8, 0.000224
16, 0.003195
32, 0.045251
64, 0.893052
128, 14.2069
256, 241.401

// GeneticSolver::run -> 
//  O(G*P*(V*E))
//  If G ~ P ~ V ~ E -> O(x^4) on average
*/

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <number_of_graphs> <graph_order>" << "\n";
    return 1; 
  }

  std::size_t num_graphs;
  std::size_t order;

  try { num_graphs = std::stoull(argv[1]); } 
  catch (const std::invalid_argument& e) {
    std::cerr << "Error: First argument is not an integer" << "\n";
    return 1;
  }
  catch (const std::out_of_range& e) {
    std::cerr << "Error: First argument is too large" << "\n";
    return 1;
  }

  try { order = std::stoull(argv[2]); } 
  catch (const std::invalid_argument& e) {
    std::cerr << "Error: Second argument is not an integer" << "\n";
    return 1;
  }
  catch (const std::out_of_range& e) {
    std::cerr << "Error: Second argument is too large" << "\n";
    return 1;
  }

  std::vector<std::string> graph_strs = generate_random_graph6(num_graphs, order);
  Graph graph;

  double total_error = 0;

  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  std::chrono::time_point<std::chrono::high_resolution_clock> end;

  for (const std::string &graph6 : graph_strs) {
    graph.from_graph6(graph6);

    start = std::chrono::high_resolution_clock::now();

    GeneticSolver solver(&graph, graph.get_order());
    VertexSet result = solver.run_set(graph.get_order());

    end = std::chrono::high_resolution_clock::now();
    double duration = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000 / 1000;

    std::cout << "Graph: " << graph6 << "\n";
    std::cout << "Observed Z(G): " << result.size() << "\n";
    std::cout << "Compute Time: " << duration << " seconds" << "\n";

    std::size_t expected = zero_forcing_wavefront(graph);
    std::int64_t delta = result.size() - expected;

    std::cout << "Real Z(G): " << expected << "\n";
    std::cout << "Delta: " << delta << "\n";
    std::cout << "\n";

    if (delta < 0) {
      std::cerr << "Crital Error: Zero forcing closure logic failed" << "\n";
      return 1;
    }

    total_error += static_cast<double>(delta) / order;
  }

  double average_accuracy = 1.0 - total_error / num_graphs;
  double average_error = total_error / num_graphs;

  std::cout << "Average Accuracy: " << average_accuracy*100 << "%" << "\n";
  std::cout << "Average Error: " << average_error*100 << "%" << "\n";
}