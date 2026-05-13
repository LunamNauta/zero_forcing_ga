#include <iostream>
#include <cassert>

#include "zero_forcing.hpp"
#include "graph.hpp"

void test_p4() {
    Graph graph;
    graph.from_graph6("Ch"); // P4: 0-1-2-3
    // An endpoint should force the whole path
    VertexSet initial = {0};
    VertexSet forced = zero_forcing_closure(graph, initial);
    assert(forced.size() == 4);
    // A middle vertex (1) has two uncolored neighbors (0, 2), so it forces nothing
    initial = {1};
    forced = zero_forcing_closure(graph, initial);
    assert(forced.size() == 1); 
    std::cout << "Test P4: Passed" << std::endl;
}

void test_c4() {
    Graph graph;
    graph.from_sparse6(":Cdo"); // C4: 0-1, 1-2, 2-3, 3-0
    // One vertex in a cycle has two uncolored neighbors, forces nothing
    VertexSet initial = {0};
    assert(zero_forcing_closure(graph, initial).size() == 1);
    // Two adjacent vertices force the whole cycle
    initial = {0, 1};
    assert(zero_forcing_closure(graph, initial).size() == 4);
    std::cout << "Test C4: Passed" << std::endl;
}

void test_star5() {
    Graph graph;
    // Star on 5: Center 0, Leaves 1,2,3,4
    graph.from_edge_list("5 4\n0 1\n0 2\n0 3\n0 4"); 
    // Center alone forces nothing (4 uncolored neighbors)
    VertexSet initial = {0};
    assert(zero_forcing_closure(graph, initial).size() == 1);
    // Center plus all but one leaf forces the last leaf
    initial = {0, 1, 2, 3};
    assert(zero_forcing_closure(graph, initial).size() == 5);
    std::cout << "Test Star 5: Passed" << std::endl;
}

void test_complete_k4() {
    Graph graph;
    graph.from_graph6("C~"); // K4
    // In Kn, you need n-1 vertices to force the last one
    VertexSet initial = {0, 1};
    assert(zero_forcing_closure(graph, initial).size() == 2);
    initial = {0, 1, 2};
    assert(zero_forcing_closure(graph, initial).size() == 4);
    std::cout << "Test K4: Passed" << std::endl;
}

void test_disjoint() {
    // 4 isolated vertices
    Graph graph(4); 
    VertexSet initial = {0, 1};
    // Isolated vertices have no neighbors to force
    assert(zero_forcing_closure(graph, initial).size() == 2);
    std::cout << "Test Disjoint: Passed" << std::endl;
}

int main() {
    test_p4();
    test_c4();
    test_star5();
    test_complete_k4();
    test_disjoint();

    std::cout << "\nAll zf_closure Tests Passed!" << std::endl;
}