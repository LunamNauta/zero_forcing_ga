
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <thread>

#include "genetic_algorithm.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"

struct SuiteResults {
    double accuracy;
    double error;
    double fail_rate;
    double avg_generations;
    std::size_t max_generations;
};

void test_graph_suite(const std::string& suite_name, const std::vector<Graph>& graphs, std::size_t population_size, std::size_t generations_before_quit) {
    if (graphs.empty()) return;

    std::cout << "========================================================\n";
    std::cout << "Starting Suite: " << suite_name << "\n";
    std::cout << "========================================================\n";

    double total_error = 0;
    double num_failed = 0;
    double total_generations = 0;
    std::size_t max_generation = 0;
    std::size_t num_tested = 0;

    for (const Graph &graph : graphs) {
        if (!graph.is_connected()) continue;

        std::cout << "--- Testing " << suite_name << " (Order: " << graph.get_order() << ") ---\n";
        auto start = std::chrono::high_resolution_clock::now();
        
        FortCoverData fc_data;
        FortCoverSolver fc(&graph, 1, 3600);
        GeneticSolver ga(&graph, 10);
        fc.set_ga_solver(&ga);
        ga.set_fc_solver(&fc);

        std::thread thread_fc(&FortCoverSolver::run, &fc, std::ref(fc_data));
        std::thread thread_ga(&GeneticSolver::run, &ga, std::numeric_limits<std::size_t>::max());

        thread_fc.join();
        ga.kill();
        thread_ga.join();
    }
}

int main(int argc, char* argv[]) {
    std::vector<Graph> paths, cycles, completes, cubics;
    for (std::size_t a = 2; a <= 32; a++) {
        paths.push_back(Graph::generate_path(a));
        if (a >= 3) cycles.push_back(Graph::generate_cycle(a));
    }
    for (std::size_t a = 3; a < 48; a++) {
        if (a >= 3) completes.push_back(Graph::generate_complete(a));
        // if (a >= 4 && a % 2 == 0) cubics.push_back(Graph::generate_cubic(a));
    }

    test_graph_suite("Paths", paths,10, 10);
    test_graph_suite("Cycles", cycles, 10, 10);
    test_graph_suite("Completes", completes, 10, 10);
    test_graph_suite("Cubics", cubics, 10, 10);

    return 0;
}