#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>

#include "../include/graph.hpp"
#include "../include/genetic_solver.hpp"

struct TrialResult {
    std::string param_name;
    double param_value;
    std::string graph_type;
    std::size_t order;
    std::size_t pop_size;
    double runtime_ms;
    std::size_t gens_to_solution;
    double time_per_gen_ms;
    std::size_t best_z;
};

TrialResult run_adaptive_trial(const std::string& param_name,
                               double param_val,
                               const std::string& graph_type,
                               const Graph& graph,
                               std::size_t pop_size,
                               std::size_t max_gens = 1000000000) 
{
    std::cout << "[Trial] Param: " << std::left << std::setw(15) << param_name << " = " << std::setw(5) << param_val
              << " | Graph: " << std::left << std::setw(10) << graph_type << " (N=" << graph.order() << ")"
              << " | Pop: " << std::setw(3) << pop_size << " ... " << std::flush;

    auto start = std::chrono::high_resolution_clock::now();

    GeneticSolver solver(&graph, pop_size);
    
    std::size_t best_z_val = solver.best_z();
    std::size_t gens_to_sol = 1;
    std::size_t stable_gens = 0;
    std::size_t total_gens = 0;
    std::size_t patience = 16;

    for (std::size_t g = 1; g <= max_gens; ++g) {
        solver.run(1);
        total_gens = g;
        std::size_t current_z = solver.best_z();
        
        if (current_z < best_z_val) {
            best_z_val = current_z;
            gens_to_sol = g;
            stable_gens = 0;
        } else {
            stable_gens++;
        }

        if (stable_gens > patience) break;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double, std::milli>(end - start).count();
    double time_per_gen = duration / static_cast<double>(total_gens > 0 ? total_gens : 1);

    std::cout << "DONE in " << std::right << std::fixed << std::setprecision(2) << std::setw(7) << duration << " ms"
              << " | Gens to Sol: " << std::setw(3) << gens_to_sol
              << " | Time/Gen: " << std::setprecision(4) << time_per_gen << " ms"
              << " | Best Z: " << std::setw(3) << best_z_val << "\n";

    return {param_name, param_val, graph_type, graph.order(), pop_size, duration, gens_to_sol, time_per_gen, best_z_val};
}

void generate_python_plotter() {
    std::ofstream py("plot_results.py");
    py << "import pandas as pd\n";
    py << "import matplotlib.pyplot as plt\n\n";
    py << "df = pd.read_csv('benchmark_parameter_results.csv')\n\n";
    
    py << "fig, axs = plt.subplots(2, 2, figsize=(16, 12))\n";
    py << "fig.suptitle('Genetic Solver Performance Analysis', fontsize=16, fontweight='bold')\n\n";
    
    py << "# 1. Generations vs Order\n";
    py << "for g_type, group in df[df['param_name'] == 'Graph_Order'].groupby('graph_type'):\n";
    py << "    axs[0, 0].plot(group['param_value'], group['gens_to_solution'], marker='o', linewidth=2, label=g_type)\n";
    py << "axs[0, 0].set_title('Graph Order (N) vs Generations to Solution (Fixed Pop=10)')\n";
    py << "axs[0, 0].set_xlabel('Graph Order (N)')\n";
    py << "axs[0, 0].set_ylabel('Generations')\n";
    py << "axs[0, 0].legend()\n";
    py << "axs[0, 0].grid(True, linestyle='--', alpha=0.7)\n\n";

    py << "# 2. Time per Generation vs Order\n";
    py << "for g_type, group in df[df['param_name'] == 'Graph_Order'].groupby('graph_type'):\n";
    py << "    axs[0, 1].plot(group['param_value'], group['time_per_gen_ms'], marker='^', linewidth=2, label=g_type)\n";
    py << "axs[0, 1].set_title('Graph Order (N) vs Time per Generation (Fixed Pop=10)')\n";
    py << "axs[0, 1].set_xlabel('Graph Order (N)')\n";
    py << "axs[0, 1].set_ylabel('Time per Generation (ms)')\n";
    py << "axs[0, 1].legend()\n";
    py << "axs[0, 1].grid(True, linestyle='--', alpha=0.7)\n\n";

    py << "pop_df = df[df['param_name'] == 'Population_Size']\n\n";

    py << "# 3. Generations vs Population Size\n";
    py << "axs[1, 0].plot(pop_df['param_value'], pop_df['gens_to_solution'], marker='s', color='green', linewidth=2, label='Random Graph (N=30)')\n";
    py << "axs[1, 0].set_title('Population Size (P) vs Generations to Solution (Fixed N=30)')\n";
    py << "axs[1, 0].set_xlabel('Population Size (P)')\n";
    py << "axs[1, 0].set_ylabel('Generations')\n";
    py << "axs[1, 0].legend()\n";
    py << "axs[1, 0].grid(True, linestyle='--', alpha=0.7)\n\n";

    py << "# 4. Time per Generation vs Population Size\n";
    py << "axs[1, 1].plot(pop_df['param_value'], pop_df['time_per_gen_ms'], marker='d', color='purple', linewidth=2, label='Random Graph (N=30)')\n";
    py << "axs[1, 1].set_title('Population Size (P) vs Time per Generation (Fixed N=30)')\n";
    py << "axs[1, 1].set_xlabel('Population Size (P)')\n";
    py << "axs[1, 1].set_ylabel('Time per Generation (ms)')\n";
    py << "axs[1, 1].legend()\n";
    py << "axs[1, 1].grid(True, linestyle='--', alpha=0.7)\n\n";

    py << "plt.tight_layout(rect=[0, 0.03, 1, 0.95])\n";
    py << "plt.savefig('comprehensive_analysis.png', dpi=300)\n";
    py << "plt.close()\n";
    py << "print('Successfully generated: comprehensive_analysis.png')\n";
    py.close();
}

int main() {
    std::ofstream csv("benchmark_parameter_results.csv");
    csv << "param_name,param_value,graph_type,order,pop_size,runtime_ms,gens_to_solution,time_per_gen_ms,best_z\n";

    auto log_res = [&](const TrialResult& r) {
        csv << r.param_name << "," << r.param_value << "," << r.graph_type << ","
            << r.order << "," << r.pop_size << "," << r.runtime_ms << ","
            << r.gens_to_solution << "," << r.time_per_gen_ms << "," << r.best_z << "\n";
    };

    std::cout << "====================================================================\n";
    std::cout << "     GENETIC SOLVER PARAMETER SWEEP & AUTO-GRAPHING GENERATOR       \n";
    std::cout << "====================================================================\n\n";

    // Experiment 1: Sweeping Order (N). Variables controlled: Population Size is strictly fixed.
    std::cout << "--- Stage 1: Sweeping Graph Order (N) ---\n";
    std::vector<std::size_t> orders;
    for (std::size_t a = 8; a <= 100; a += 4) orders.push_back(a); // Stepped loop for cleaner graphs

    std::size_t default_pop = 10;

    for (std::size_t N : orders) {
        log_res(run_adaptive_trial("Graph_Order", N, "Path", GraphGenerator::path(N), default_pop));
        log_res(run_adaptive_trial("Graph_Order", N, "Cycle", GraphGenerator::cycle(N), default_pop));
        log_res(run_adaptive_trial("Graph_Order", N, "Complete", GraphGenerator::complete(N), default_pop));
        log_res(run_adaptive_trial("Graph_Order", N, "Random", GraphGenerator::random(N, 1, 0.4)[0], default_pop));
    }

    // Experiment 2: Sweeping Population Size (P). Variables controlled: Graph Order and Graph Instance are strictly fixed.
    std::cout << "\n--- Stage 2: Sweeping Population Size (P) ---\n";
    std::vector<std::size_t> pop_sizes;
    for (std::size_t a = 8; a <= 128; a += 4) pop_sizes.push_back(a);

    std::size_t fixed_order = 30;
    Graph fixed_rand_graph = GraphGenerator::random(fixed_order, 1, 0.4)[0];

    for (std::size_t P : pop_sizes) {
        log_res(run_adaptive_trial("Population_Size", P, "Random", fixed_rand_graph, P));
    }

    csv.close();
    std::cout << "\nBenchmark data saved to 'benchmark_parameter_results.csv'.\n";
    
    generate_python_plotter();
    std::cout << "Executing Python plotting script...\n";
    int ret = std::system("python3 plot_results.py");
    
    if (ret == 0) {
        std::cout << "Graph generated successfully: 'comprehensive_analysis.png'.\n";
    } else {
        std::cout << "Note: Python execution encountered an issue. Run 'python3 plot_results.py' manually.\n";
    }
    std::cout << "====================================================================\n";

    return 0;
}
