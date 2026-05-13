#include <iostream>
#include <cassert>

#include "graph.hpp"

void test_empty() {
  Graph graph;
  // "" -> empty
  graph.from_graph6("");
  assert(graph.get_order() == 0);
  assert(graph.get_edge_count() == 0);
  std::cout << "Test Empty: Passed" << "\n";
}

void test_single() {
  Graph graph;
  // "@" -> single
  graph.from_graph6("@"); 
  assert(graph.get_order() == 1);
  assert(graph.get_edge_count() == 0);
  std::cout << "Test Single: Passed" << "\n";
}

void test_k3() {
  Graph graph;
  // "Bw" -> K3
  graph.from_graph6("Bw"); 
  assert(graph.get_order() == 3);
  assert(graph.get_edge_count() == 3);
  assert(graph.has_edge(0, 1));
  assert(graph.has_edge(0, 2));
  assert(graph.has_edge(1, 2));
  std::cout << "Test K3: Passed" << "\n";
}

void test_p4() {
  Graph graph;
  // "Ch" -> P4
  graph.from_graph6("Ch");
  assert(graph.get_order() == 4);
  assert(graph.get_edge_count() == 3);
  assert(graph.has_edge(0, 1));
  assert(graph.has_edge(1, 2));
  assert(graph.has_edge(2, 3));
  assert(!graph.has_edge(0, 3));
  std::cout << "Test P4: Passed" << "\n";
}

void test_star5() {
  Graph graph;
  // "Ds_" -> star on 5
  graph.from_graph6("Ds_");
  assert(graph.get_order() == 5);
  assert(graph.get_edge_count() == 4);
  for (std::size_t a = 1; a <= 4; a++) {
    assert(graph.has_edge(0, a));
  }
  std::cout << "Test Star5: Passed" << "\n";
}

void test_bull5() {
  Graph graph;
  // "DyG" -> bull on 5
  graph.from_graph6("DyG");
  assert(graph.get_order() == 5);
  assert(graph.get_edge_count() == 5);
  // Bull graph contains a triangle and two pendant "horns"
  assert(graph.get_pendants().size() == 2);
  std::cout << "Test Bull: Passed" << "\n";
}

int main() {
  test_empty();
  test_single();
  test_k3();
  test_p4();
  test_star5();
  test_bull5();
        
  std::cout << "\nAll Graph6 Parsing Tests Passed!" << "\n";
}