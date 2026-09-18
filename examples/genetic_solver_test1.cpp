#include <chrono>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>
#include <cassert>

#include "genetic_solver.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"

struct PerGraphResults {
    std::size_t generations;
    std::size_t distance;
};

struct SuiteResults {
    std::vector<PerGraphResults> results;
    double accuracy;
    double error;
    double fail_rate;
    double avg_generations;
    std::size_t max_generations;
};

SuiteResults test_graph_suite(const std::string& suite_name, const std::vector<Graph>& graphs, std::size_t population_size, std::size_t generations_before_quit) {
    if (graphs.empty()) return {{}, 0, 0, 0, 0, 0};

    std::cout << "========================================================\n";
    std::cout << "Starting Suite: " << suite_name << "\n";
    std::cout << "========================================================\n";

    std::vector<PerGraphResults> results;
    double total_error = 0;
    double num_failed = 0;
    double total_generations = 0;
    std::size_t max_generation = 0;
    std::size_t num_tested = 0;

    for (const Graph &graph : graphs) {
        if (!GraphMiscInfo::is_connected(graph)) continue;

        std::cout << "--- Testing " << suite_name << " (Order: " << graph.order() << ") ---\n";
        auto start = std::chrono::high_resolution_clock::now();
        GeneticSolver solver(&graph, population_size);
        
        std::size_t generation = 0;
        std::size_t generation_found = 0;
        double best_variance = solver.variance();

        while (true) {
            solver.run(1);
            generation++;
            if (solver.since_better_z() == 0) generation_found = generation;
            
            double current_variance = solver.variance();
            if (current_variance < 0.5 / population_size || solver.since_better_score() > generations_before_quit) break;
        }

        std::cout << "Generation Found: " << generation_found << "\n";
        std::size_t expected = zero_forcing_wavefront(graph);
        std::int64_t delta = solver.best_z() - expected;
        std::cout << "Distance: " << delta << "\n\n";
        
        if (delta < 0) { std::cerr << "Critical Error: Logic failed\n"; exit(1); }
        if (delta != 0) num_failed++;

        results.push_back({generation_found, static_cast<std::size_t>(delta)});
        total_error += static_cast<double>(delta) / graph.order();
        total_generations += generation_found;
        max_generation = std::max(max_generation, generation_found);
        num_tested++;
    }

    double acc = 1.0 - (total_error / num_tested);
    double err = total_error / num_tested;
    double fail = num_failed / num_tested;
    double avg_gen = total_generations / num_tested;

    std::cout << "--- " << suite_name << " Results ---\n";
    std::cout << "Accuracy: " << acc * 100 << "%, Fail Rate: " << fail * 100 << "%, Avg Gens: " << avg_gen << "\n\n";

    return {results, acc, err, fail, avg_gen, max_generation};
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pop_size> <gens_quit> <max_order>\n";
        return 1; 
    }

    std::size_t population_size = std::stoull(argv[1]);
    std::size_t generations_before_quit = std::stoull(argv[2]);
    std::size_t max_order = std::stoull(argv[3]);

    std::vector<Graph> randoms, paths, cycles, completes, cubics;
    for (std::size_t a = 2; a <= max_order; a++) {
        for (double p = 0.5; p < 0.9; p++) {
          std::vector<Graph> graphs = GraphGenerator::random(a, 16, p);
          randoms.insert(randoms.end(), graphs.begin(), graphs.end());
        }
        paths.push_back(GraphGenerator::path(a));
        if (a >= 3) cycles.push_back(GraphGenerator::cycle(a));
        if (a >= 3) completes.push_back(GraphGenerator::complete(a));
        // if (a >= 4 && a % 2 == 0) cubics.push_back(GraphGenerator::cubic(a));
    }

    std::vector<SuiteResults> results = {
        test_graph_suite("Randoms", randoms, population_size, generations_before_quit),
        test_graph_suite("Paths", paths, population_size, generations_before_quit),
        test_graph_suite("Cycles", cycles, population_size, generations_before_quit),
        test_graph_suite("Completes", completes, population_size, generations_before_quit),
        // test_graph_suite("Cubics", cubics, population_size, generations_before_quit)
    };

    double g_acc = 0;
    double g_err = 0;
    double g_fail = 0;
    double total = 0;
    for (const auto& r : results) {
      std::cout << r.accuracy << "\n";
        g_acc += (float)r.accuracy;
        g_err += (float)r.error;
        g_fail += (float)r.fail_rate;
        total += r.results.size();
    }
    g_acc /= results.size();
    g_err /= results.size();
    g_fail /= results.size();

    std::cout << "========================================================\n";
    std::cout << "RANDOMS PERFORMANCE SUMMARY\n";
    std::cout << "========================================================\n";
    std::cout << "Average Accuracy: " << results[0].accuracy * 100 << "%\n";
    std::cout << "Average Error: " << results[0].error * 100 << "%\n";
    std::cout << "Global Fail Rate: " << results[0].fail_rate * 100 << "%\n";
    std::cout << "\n";

    std::cout << "========================================================\n";
    std::cout << "PATHS PERFORMANCE SUMMARY\n";
    std::cout << "========================================================\n";
    std::cout << "Average Accuracy: " << results[1].accuracy * 100 << "%\n";
    std::cout << "Average Error: " << results[1].error * 100 << "%\n";
    std::cout << "Global Fail Rate: " << results[1].fail_rate * 100 << "%\n";
    std::cout << "\n";

    std::cout << "========================================================\n";
    std::cout << "CYCLES PERFORMANCE SUMMARY\n";
    std::cout << "========================================================\n";
    std::cout << "Average Accuracy: " << results[2].accuracy * 100 << "%\n";
    std::cout << "Average Error: " << results[2].error * 100 << "%\n";
    std::cout << "Global Fail Rate: " << results[2].fail_rate * 100 << "%\n";
    std::cout << "\n";

    std::cout << "========================================================\n";
    std::cout << "COMPLETES PERFORMANCE SUMMARY\n";
    std::cout << "========================================================\n";
    std::cout << "Average Accuracy: " << results[3].accuracy * 100 << "%\n";
    std::cout << "Average Error: " << results[3].error * 100 << "%\n";
    std::cout << "Global Fail Rate: " << results[3].fail_rate * 100 << "%\n";
    std::cout << "\n";

    std::cout << "========================================================\n";
    std::cout << "GLOBAL PERFORMANCE SUMMARY\n";
    std::cout << "========================================================\n";
    std::cout << "Average Accuracy: " << g_acc * 100 << "%\n";
    std::cout << "Average Error: " << g_err * 100 << "%\n";
    std::cout << "Global Fail Rate: " << g_fail * 100 << "%\n";
    std::cout << "\n";

    return 0;
}
