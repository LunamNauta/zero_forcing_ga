#include <iostream>
#include <cassert>

#include "graph.hpp"

void test_empty() {
  Graph graph;
  // : -> empty
  graph.from_sparse6(":");
  assert(graph.get_order() == 0);
  std::cout << "Test Empty: Passed" << "\n";
}

void test_single() {
  Graph graph;
  // ":@" -> single
  graph.from_sparse6(":@");   
  assert(graph.get_order() == 1);
  assert(graph.get_edge_count() == 0);
  std::cout << "Test Single: Passed" << "\n";
}

void test_k2() {
  Graph graph;
  // ":A_" -> k2
  graph.from_sparse6(":A_");
  assert(graph.get_order() == 2);
  assert(graph.get_edge_count() == 1);
  assert(graph.has_edge(0, 1));
  std::cout << "Test K2: Passed" << "\n";
}

void test_p3() {
  Graph graph;
  // ":B_O" -> p3
  graph.from_sparse6(":Bd");
  assert(graph.get_order() == 3);
  assert(graph.get_edge_count() == 2);
  assert(graph.has_edge(0, 1));
  assert(graph.has_edge(1, 2));
  std::cout << "Test P3: Passed" << "\n";
}

void test_star5() {
  Graph graph;
  // ":DaG_ -> star5"
  graph.from_sparse6(":DaG_");
  assert(graph.get_order() == 5);
  assert(graph.get_edge_count() == 4);
  for (std::size_t a = 1; a <= 4; a++) {
    assert(graph.has_edge(0, a));
  }
  std::cout << "Test Star5: Passed" << "\n";
}

void test_c4() {
  Graph graph;
  // ":Cdo" -> C4
  graph.from_sparse6(":Cdo");
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

  std::cout << "\nAll Sparse6 Parsing Tests Passed!" << "\n";
}