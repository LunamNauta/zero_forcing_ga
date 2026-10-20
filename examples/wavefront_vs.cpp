#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cstdlib>

#include "../include/graph.hpp"

// Referencing verbatim to ensure linking with the exact implementation file provided
// Note: Depending on your build system, you may either link this during compilation
// or include it directly if doing a unity build.
#include "../include/zero_forcing.hpp"

int main() {
    std::cout << "====================================================================\n";
    std::cout << "     ZERO FORCING WAVEFRONT ALGORITHM COMPARISON BENCHMARK          \n";
    std::cout << "====================================================================\n\n";

    std::vector<std::size_t> orders = {2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128};
    std::vector<double> probabilities = {0.9};
    const std::size_t TRIALS_PER_CONFIG = 5;

    std::cout << std::left << std::setw(10) << "Order(N)" 
              << std::setw(10) << "Prob(p)" 
              << std::setw(15) << "New Time(ms)" 
              << std::setw(15) << "Old Time 1 (ms)" 
              << std::setw(15) << "Old Time 2 (ms)" 
              << std::setw(20) << "Speedup 1" 
              << std::setw(20) << "Speedup 2" 
              << std::setw(20) << "Status\n";
    std::cout << std::string(75, '-') << "\n";

    for (std::size_t N : orders) {
        for (double p : probabilities) {
            
            // Generate deterministic sets of random graphs for consistent comparison
            std::vector<Graph> graphs = GraphGenerator::random(N, TRIALS_PER_CONFIG, p);

            double total_time_new = 0.0;
            double total_time_old1 = 0.0;
            double total_time_old2 = 0.0;
            bool all_matched = true;

            for (const Graph& g : graphs) {
                // Benchmark new optimized algorithm
                auto start_new = std::chrono::high_resolution_clock::now();
                std::size_t z_new = zero_forcing_wavefront(g, g.order());
                auto end_new = std::chrono::high_resolution_clock::now();
                total_time_new += std::chrono::duration<double, std::milli>(end_new - start_new).count();

                // Benchmark old algorithm
                auto start_old1 = std::chrono::high_resolution_clock::now();
                std::size_t z_old1 = zero_forcing_wavefront_old1(g);
                auto end_old1 = std::chrono::high_resolution_clock::now();
                total_time_old1 += std::chrono::duration<double, std::milli>(end_old1 - start_old1).count();

                auto start_old2 = std::chrono::high_resolution_clock::now();
                std::size_t z_old2 = zero_forcing_wavefront_old2(g);
                auto end_old2 = std::chrono::high_resolution_clock::now();
                total_time_old2 += std::chrono::duration<double, std::milli>(end_old2 - start_old2).count();

                // Value Verification
                if (z_new != z_old1 || z_new != z_old2) {
                    std::cerr << "\n[FATAL ERROR] Value mismatch!\n"
                              << "Graph Order: " << N << " | Edge Prob: " << p << "\n"
                              << "zero_forcing_wavefront returned: " << z_new << "\n"
                              << "zero_forcing_wavefront_old1 returned: " << z_old1 << "\n"
                              << "zero_forcing_wavefront_old2 returned: " << z_old2 << "\n";
                    std::exit(EXIT_FAILURE);
                }
            }

            double avg_new = total_time_new / TRIALS_PER_CONFIG;
            double avg_old1 = total_time_old1 / TRIALS_PER_CONFIG;
            double avg_old2 = total_time_old2 / TRIALS_PER_CONFIG;
            double speedup1 = avg_old1 / (avg_new > 0 ? avg_new : 1e-9);
            double speedup2 = avg_old2 / (avg_new > 0 ? avg_new : 1e-9);

            std::cout << std::left << std::setw(10) << N 
                      << std::setw(10) << p 
                      << std::setw(15) << std::fixed << std::setprecision(4) << avg_new 
                      << std::setw(15) << std::fixed << std::setprecision(4) << avg_old1
                      << std::setw(15) << std::fixed << std::setprecision(4) << avg_old2
                      << "    " << std::setw(10) << std::setprecision(2) << (std::to_string(speedup1) + "x")
                      << "    " << std::setw(10) << std::setprecision(2) << (std::to_string(speedup2) + "x")
                      << " " << (all_matched ? "[PASS]" : "[FAIL]") << "\n";
        }
    }

    std::cout << "====================================================================\n";
    std::cout << "Benchmark complete. All return values exactly matched.\n";

    return 0;
}
