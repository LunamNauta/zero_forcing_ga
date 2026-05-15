#include <chrono>
#include <iostream>
#include <cassert>

#include <limits>
#include <nauty/naututil.h>
#include <nauty/gtools.h>
#include <nauty/nauty.h>

#include "genetic_algorithm.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"

std::vector<std::string> generate_random_graph6(int num_graphs, int n) {
  srand(time(0));
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

  std::vector<std::string> graph_strs = generate_random_graph6(num_graphs, order);
  Graph graph;

  double total_error = 0;
  double num_failed = 0;
  std::size_t generation = 0;

  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  std::chrono::time_point<std::chrono::high_resolution_clock> end;

  for (const std::string &graph6 : graph_strs) {
    graph.from_graph6(graph6);
    std::cout << "Graph: " << graph6 << "\n";

    start = std::chrono::high_resolution_clock::now();

    GeneticSolver solver(&graph, population_size);
    VertexSet result;
    std::size_t best = std::numeric_limits<std::size_t>::max();
 
    generation = 0;
    while (true) {
      result = std::move(solver.run_set(1));
      if (result.size() < best) {
        best = result.size();
        std::cout << "Gen " << generation++ << ": New Best -> " << best << "\n";
      }
      if (solver.time_since_best() > generations_before_quit) break;
    }

    end = std::chrono::high_resolution_clock::now();
    double duration = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000 / 1000;

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
    if (delta != 0) num_failed++;
  }

  double average_accuracy = 1.0 - total_error / num_graphs;
  double average_error = total_error / num_graphs;
  double fail_rate = num_failed / num_graphs;

  std::cout << "Average Accuracy: " << average_accuracy*100 << "%" << "\n";
  std::cout << "Average Error: " << average_error*100 << "%" << "\n";
  std::cout << "Fail Rate: " << fail_rate*100 << "%" << "\n";
}