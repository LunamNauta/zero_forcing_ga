#include <iostream>
#include <cassert>

#include "genetic_algorithm.hpp"
#include "graph.hpp"

void test_empty() {
  Graph graph;
  // : -> empty
  graph.from_sparse6(":");
  
  ZFGeneticSolver solver(&graph, 100);
  VertexSet result = solver.run(100);

  assert(graph.is_valid_zf(result));
  assert(result.size() == 0);

  std::cout << "Test Empty: Passed" << "\n";
}

void test_single() {
  Graph graph;
  // ":@" -> single
  graph.from_sparse6(":@");   

  ZFGeneticSolver solver(&graph, 100);
  VertexSet result = solver.run(100);

  assert(graph.is_valid_zf(result));
  assert(result.size() == 1);

  std::cout << "Test Single: Passed" << "\n";
}

void test_k2() {
  Graph graph;
  // ":A_" -> k2
  graph.from_sparse6(":A_");

  ZFGeneticSolver solver(&graph, 100);
  VertexSet result = solver.run(100);

  assert(graph.is_valid_zf(result));
  assert(result.size() == 1);

  std::cout << "Test K2: Passed" << "\n";
}

void test_p3() {
  Graph graph;
  // ":B_O" -> p3
  graph.from_sparse6(":Bd");

  ZFGeneticSolver solver(&graph, 100);
  VertexSet result = solver.run(100);

  assert(graph.is_valid_zf(result));
  assert(result.size() == 1);

  std::cout << "Test P3: Passed" << "\n";
}

void test_star5() {
  Graph graph;
  // ":DaG_ -> star5"
  graph.from_sparse6(":DaG_");

  ZFGeneticSolver solver(&graph, 100);
  VertexSet result = solver.run(100);

  assert(graph.is_valid_zf(result));
  assert(result.size() == 3);

  std::cout << "Test Star5: Passed" << "\n";
}

void test_c4() {
  Graph graph;
  // ":Cdo" -> C4
  graph.from_sparse6(":Cdo");

  ZFGeneticSolver solver(&graph, 100);
  VertexSet result = solver.run(100);

  assert(graph.is_valid_zf(result));
  assert(result.size() == 2);

  std::cout << "Test Cycle 4: Passed" << "\n";
}

int main() {
  test_empty();
  test_single();
  test_k2();
  test_p3();
  test_star5();
  test_c4();

  std::cout << "\nAll Genetic Algorithm Tests Completed!" << "\n";
}