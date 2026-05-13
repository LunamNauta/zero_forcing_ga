#include <nauty/naututil.h>
#include <nauty/gtools.h>
#include <nauty/nauty.h>

#include <iostream>
#include <vector>
#include <string>

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
  srand(time(nullptr));

  int order = 50;
  int count = 100;

  std::vector<std::string> graphs = generate_random_graph6(count, order);

  std::cout << "Generated " << graphs.size() << " graphs." << std::endl;
  std::cout << "First 5 examples:" << std::endl;
  for (std::size_t a = 0; a < 5; a++) {
    std::cout << graphs[a] << std::endl;
  }

  return 0;
}