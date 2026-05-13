#include "graph.hpp"

#include <iostream>
#include <cassert>

void test_empty() {
  Graph graph;
  // "" -> empty
  graph.from_edge_list("");
  assert(graph.get_order() == 0);
  assert(graph.get_edge_count() == 0);
  std::cout << "Test Empty: Passed" << "\n";
}

void test_single() {
  Graph graph;
  // "1 0" -> single
  graph.from_edge_list("1 0");
  assert(graph.get_order() == 1);
  assert(graph.get_edge_count() == 0);
  std::cout << "Test Single: Passed" << "\n";
}

void test_k2() {
  Graph graph;
  // "2 1\n0 1" -> k2
  graph.from_edge_list("2 1\n0 1");
  assert(graph.get_order() == 2);
  assert(graph.get_edge_count() == 1);
  assert(graph.has_edge(0, 1));
  std::cout << "Test K2: Passed" << "\n";
}

void test_p3() {
  Graph graph;
  // "3 2\n0 1\n1 2" -> p3
  graph.from_edge_list("3 2\n0 1\n1 2");
  assert(graph.get_order() == 3);
  assert(graph.get_edge_count() == 2);
  assert(graph.has_edge(0, 1));
  assert(graph.has_edge(1, 2));
  std::cout << "Test P3: Passed" << "\n";
}

void test_star5() {
  Graph graph;
  // "5 4\n0 1\n0 2\n0 3\n0 4" -> star5
  graph.from_edge_list("5 4\n0 1\n0 2\n0 3\n0 4");
  assert(graph.get_order() == 5);
  assert(graph.get_edge_count() == 4);
  for (std::size_t a = 1; a <= 4; a++) {
    assert(graph.has_edge(0, a));
  }
  std::cout << "Test Star5: Passed" << "\n";
}

void test_c4() {
  Graph graph;
  // "4 4\n0 1\n1 2\n2 3\n3 0" -> c4
  graph.from_edge_list("4 4\n0 1\n1 2\n2 3\n3 0");
  assert(graph.get_order() == 4);
  assert(graph.get_edge_count() == 4);
  assert(graph.has_edge(0, 1));
  assert(graph.has_edge(1, 2));
  assert(graph.has_edge(2, 3));
  assert(graph.has_edge(3, 0));
  std::cout << "Test Cycle 4: Passed" << "\n";
}

int main() {
  test_empty();
  test_single();
  test_k2();
  test_p3();
  test_star5();
  test_c4();

  std::cout << "\nAll Edge List Parsing Tests Passed!" << "\n";
  return 0;
}