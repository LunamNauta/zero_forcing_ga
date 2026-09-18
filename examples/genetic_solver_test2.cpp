#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <iomanip>

#include "../include/graph.hpp"
#include "../include/genetic_solver.hpp"

struct BenchmarkResult {
    std::string exp_type;
    std::string graph_type;
    std::size_t order;
    std::size_t pop_size;
    std::size_t generations;
    double runtime_ms;
    std::size_t best_z;
};

// Benchmark runner with detailed real-time logging
BenchmarkResult run_trial(std::size_t current_trial,
                          std::size_t total_trials,
                          const std::string& exp_type,
                          const std::string& type_name,
                          const Graph& graph,
                          std::size_t pop_size,
                          std::size_t generations) 
{
    // Log trial setup before starting solver
    std::cout << "[" << std::setw(2) << current_trial << "/" << total_trials << "] "
              << "Exp: " << std::left << std::setw(12) << exp_type
              << " | Type: " << std::left << std::setw(12) << type_name
              << " | N=" << std::setw(3) << graph.order()
              << " | P=" << std::setw(3) << pop_size
              << " | G=" << std::setw(3) << generations
              << " ... " << std::flush; // std::flush ensures log outputs immediately

    auto start = std::chrono::high_resolution_clock::now();
    
    GeneticSolver solver(&graph, pop_size);
    solver.run(generations);
    
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double, std::milli>(end - start).count();

    // Log trial results immediately upon completion
    std::cout << "DONE in " << std::right << std::fixed << std::setprecision(2) 
              << std::setw(7) << duration << " ms"
              << " | Best Z: " << std::setw(3) << solver.best_z() 
              << " | Var: " << std::setprecision(4) << solver.variance() << "\n";

    return {exp_type, type_name, graph.order(), pop_size, generations, duration, solver.best_z()};
}

int main() {
    std::ofstream csv("benchmark_results.csv");
    csv << "exp_type,graph_type,order,pop_size,generations,runtime_ms,best_z\n";

    auto write_result = [&](const BenchmarkResult& r) {
        csv << r.exp_type << "," << r.graph_type << "," << r.order << ","
            << r.pop_size << "," << r.generations << "," << r.runtime_ms << ","
            << r.best_z << "\n";
    };

    // Configuration sets
    std::vector<std::size_t> orders = {10, 20, 30, 40, 50, 60, 80, 100};
    std::vector<std::size_t> pop_sizes = {10, 20, 50, 100, 200, 400};
    std::vector<std::size_t> generations_list = {10, 25, 50, 100, 200, 500};
    std::size_t default_pop = 10;
    std::size_t default_gens = 10;

    // Calculate total trials for progress tracking
    std::size_t exp1_trials = 0;
    for (std::size_t N : orders) {
        exp1_trials += 4; // Path, Cycle, Complete, Random
        // if (N >= 4 && N % 2 == 0) exp1_trials++; // Cubic
    }
    std::size_t total_trials = exp1_trials + pop_sizes.size() + generations_list.size();
    std::size_t trial_counter = 1;

    std::cout << "====================================================================================\n";
    std::cout << "                  GENETIC SOLVER BENCHMARK SUITE - VERBOSE LOG                       \n";
    std::cout << "====================================================================================\n";
    std::cout << "Total Scheduled Trials: " << total_trials << "\n\n";

    // -------------------------------------------------------------
    // Experiment 1: Graph Order
    // -------------------------------------------------------------
    std::cout << "--- Stage 1: Varying Graph Order (N) [Default P=" << default_pop << ", G=" << default_gens << "] ---\n";
    for (std::size_t N : orders) {
        Graph path = GraphGenerator::path(N);
        write_result(run_trial(trial_counter++, total_trials, "Order", "Path", path, default_pop, default_gens));

        Graph cycle = GraphGenerator::cycle(N);
        write_result(run_trial(trial_counter++, total_trials, "Order", "Cycle", cycle, default_pop, default_gens));

        Graph complete = GraphGenerator::complete(N);
        write_result(run_trial(trial_counter++, total_trials, "Order", "Complete", complete, default_pop, default_gens));

        Graph random_g = GraphGenerator::random(N, 1, 0.8)[0];
        write_result(run_trial(trial_counter++, total_trials, "Order", "Random_p0.8", random_g, default_pop, default_gens));

        // if (N >= 4 && N % 2 == 0) {
        //     Graph cubic_g = GraphGenerator::cubic(N, 1)[0];
        //     write_result(run_trial(trial_counter++, total_trials, "Order", "Cubic", cubic_g, default_pop, default_gens));
        // }
    }

    // -------------------------------------------------------------
    // Experiment 2: Population Size
    // -------------------------------------------------------------
    std::cout << "\n--- Stage 2: Varying Population Size (P) [Fixed Random Graph N=40, G=" << default_gens << "] ---\n";
    Graph fixed_graph_pop = GraphGenerator::random(40, 1, 0.3)[0];
    for (std::size_t P : pop_sizes) {
        write_result(run_trial(trial_counter++, total_trials, "Population", "Random_N40", fixed_graph_pop, P, default_gens));
    }

    // -------------------------------------------------------------
    // Experiment 3: Generation Count
    // -------------------------------------------------------------
    std::cout << "\n--- Stage 3: Varying Generations (G) [Fixed Random Graph N=40, P=" << default_pop << "] ---\n";
    Graph fixed_graph_gen = GraphGenerator::random(40, 1, 0.3)[0];
    for (std::size_t G : generations_list) {
        write_result(run_trial(trial_counter++, total_trials, "Generations", "Random_N40", fixed_graph_gen, default_pop, G));
    }

    csv.close();
    std::cout << "\n====================================================================================\n";
    std::cout << "All " << total_trials << " trials successfully executed. Output saved to 'benchmark_results.csv'.\n";
    std::cout << "====================================================================================\n";

    return 0;
}
