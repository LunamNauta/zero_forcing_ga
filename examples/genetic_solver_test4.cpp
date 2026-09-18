#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>

#include "../include/graph.hpp"
#include "../include/genetic_solver.hpp"

namespace fs = std::filesystem;

void generate_python_plotter(const fs::path& csv_path) {
    fs::path py_path = csv_path.parent_path() / "plot_trajectory.py";
    fs::path img_path = csv_path.parent_path() / "trajectory_analysis.png";

    std::ofstream py(py_path);
    py << "import pandas as pd\n";
    py << "import matplotlib.pyplot as plt\n\n";
    py << "df = pd.read_csv('" << csv_path.filename().string() << "')\n\n";
    
    py << "plt.figure(figsize=(15, 5))\n\n";
    
    // 1. Population Variance over Generations
    py << "plt.subplot(1, 3, 1)\n";
    py << "for g_type, group in df.groupby('graph_type'):\n";
    py << "    plt.plot(group['generation'], group['variance'], label=g_type, linewidth=2)\n";
    py << "plt.title('Population Variance vs Generation', fontsize=11, fontweight='bold')\n";
    py << "plt.xlabel('Generation')\n";
    py << "plt.ylabel('Population Variance')\n";
    py << "plt.legend()\n";
    py << "plt.grid(True, linestyle='--', alpha=0.7)\n\n";
    
    // 2. Best Score over Generations
    py << "plt.subplot(1, 3, 2)\n";
    py << "for g_type, group in df.groupby('graph_type'):\n";
    py << "    plt.plot(group['generation'], group['best_score'], label=g_type, linewidth=2)\n";
    py << "plt.title('Best Score vs Generation', fontsize=11, fontweight='bold')\n";
    py << "plt.xlabel('Generation')\n";
    py << "plt.ylabel('Best Score')\n";
    py << "plt.legend()\n";
    py << "plt.grid(True, linestyle='--', alpha=0.7)\n\n";
    
    // 3. Best Z(G) over Generations
    py << "plt.subplot(1, 3, 3)\n";
    py << "for g_type, group in df.groupby('graph_type'):\n";
    py << "    plt.plot(group['generation'], group['best_z'], label=g_type, linewidth=2, marker='.')\n";
    py << "plt.title('Best Z(G) vs Generation', fontsize=11, fontweight='bold')\n";
    py << "plt.xlabel('Generation')\n";
    py << "plt.ylabel('Best Zero Forcing Number Z(G)')\n";
    py << "plt.legend()\n";
    py << "plt.grid(True, linestyle='--', alpha=0.7)\n\n";
    
    py << "plt.tight_layout()\n";
    py << "plt.savefig('" << img_path.filename().string() << "', dpi=300)\n";
    py << "plt.close()\n";
    py << "print('Successfully generated: " << img_path.filename().string() << "')\n";
    py.close();
}

int main() {
    fs::path csv_path = "trajectory_results.csv";
    std::ofstream csv(csv_path);
    csv << "graph_type,generation,best_score,best_z,variance\n";

    std::size_t default_order = 64;

    std::vector<std::pair<std::string, Graph>> test_cases = {
        {"Path", GraphGenerator::path(default_order)},
        {"Cycle", GraphGenerator::cycle(default_order)},
        {"Complete", GraphGenerator::complete(default_order)},
        {"Random", GraphGenerator::random(default_order, 1, 0.4)[0]}
    };

    std::size_t total_gens = 128;
    std::size_t pop_size = 10;

    std::cout << "====================================================================\n";
    std::cout << "     GENETIC SOLVER TRAJECTORY & METRIC TRACKING TEST               \n";
    std::cout << "====================================================================\n\n";

    for (auto& [name, graph] : test_cases) {
        std::cout << "Tracking evolution metrics for: " << name << " ... " << std::flush;
        GeneticSolver solver(&graph, pop_size);

        // Record initial state (Generation 0)
        csv << name << ",0," << solver.best_score() << "," << solver.best_z() << "," << solver.variance() << "\n";

        std::cout << "Finding exact... ";
        std::size_t exact = zero_forcing_wavefront(graph);
        std::cout << "Found... ";

        for (std::size_t g = 1; g <= total_gens; ++g) {
            // solver.run(1);
            if (solver.best_z() == exact) break;
            csv << name << "," << g << "," << solver.best_score() << "," << solver.best_z() << "," << solver.variance() << "\n";
        }
        std::cout << "DONE\n";
    }

    csv.close();
    std::cout << "\nTrajectory metrics saved to " << fs::absolute(csv_path) << ".\n";
    
    // Generate and run Python plotting script
    generate_python_plotter(csv_path);
    std::cout << "Executing Python plotting script...\n";
    std::string cmd = "python3 " + (csv_path.parent_path() / "plot_trajectory.py").string();
    int ret = std::system(cmd.c_str());
    
    if (ret == 0) {
        std::cout << "Graph generated successfully: 'trajectory_analysis.png'.\n";
    } else {
        std::cout << "Note: Python execution encountered an issue. Ensure pandas and matplotlib are installed.\n";
    }
    std::cout << "====================================================================\n";

    return 0;
}
