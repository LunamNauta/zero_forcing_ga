#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "genetic_solver.hpp"
#include "graph.hpp"
#include "zero_forcing.hpp"

struct PerGraphResults {
  std::size_t generations = 0;
  std::size_t distance = 0;
  std::size_t expected_z = 0;
  std::size_t best_z = 0;
};

struct SuiteResults {
  std::vector<PerGraphResults> results;
  double accuracy = 0.0;
  double error = 0.0;
  double fail_rate = 0.0;
  double avg_generations = 0.0;
  std::size_t max_generations = 0;
  std::size_t tested = 0;
  std::size_t skipped = 0;
};

struct TestParams {
  std::size_t n;
  double p;
  std::size_t graph_count;
  std::size_t population_size;
  std::size_t generations_before_quit;
  std::size_t max_generations;
  std::size_t r_min;
  std::size_t r_max;
  std::size_t seed;
};

static SuiteResults test_graph_suite(
  const std::vector<Graph>& graphs,
  std::size_t population_size,
  std::size_t generations_before_quit,
  std::size_t max_generations
) {
  SuiteResults suite;

  if (graphs.empty()) {
    return suite;
  }

  double total_error = 0.0;
  std::size_t num_failed = 0;
  double total_generations = 0.0;

  for (const Graph& graph : graphs) {
    if (!GraphMiscInfo::is_connected(graph)) {
      ++suite.skipped;
      continue;
    }

    // The configuration constructor takes psi.
    GeneticSolverConfig config(population_size);

    GeneticSolver solver(&graph, config);

    std::size_t generation = 0;
    std::size_t generation_found = 0;

    while (generation < max_generations) {
      solver.run(1);
      ++generation;

      const GeneticSolverInfo info = solver.info();

      // Record the generation in which best_z most recently improved.
      if (info.since_better_z() == 0) {
        generation_found = generation;
      }

      if (info.variance() < 0.5 / population_size ||
        info.since_better_score() > generations_before_quit) {
        break;
      }
    }

    const std::size_t expected = zero_forcing_wavefront(graph);
    const std::size_t actual = solver.info().best_z();

    // Avoid unsigned subtraction: it would wrap if actual < expected.
    if (actual < expected) {
      std::cerr
        << "ERROR: solver result is below the reference value!\n"
        << "Order: " << graph.order()
        << ", expected: " << expected
        << ", actual: " << actual << '\n';
      std::exit(EXIT_FAILURE);
    }

    const std::size_t distance = actual - expected;

    suite.results.push_back({
      generation_found,
      distance,
      expected,
      actual
    });

    total_error += static_cast<double>(distance) /
             static_cast<double>(graph.order());

    total_generations += static_cast<double>(generation_found);
    suite.max_generations =
      std::max(suite.max_generations, generation_found);

    if (distance != 0) {
      ++num_failed;
    }

    ++suite.tested;
  }

  if (suite.tested == 0) {
    return suite;
  }

  suite.error = total_error / static_cast<double>(suite.tested);
  suite.accuracy = 1.0 - suite.error;
  suite.fail_rate =
    static_cast<double>(num_failed) /
    static_cast<double>(suite.tested);
  suite.avg_generations =
    total_generations / static_cast<double>(suite.tested);

  return suite;
}

static void print_results(
  std::size_t r,
  const TestParams& params,
  const SuiteResults& result,
  double elapsed_seconds
) {
  std::cout << "\n====================================================\n";
  std::cout << "Results for r = " << r << '\n';
  std::cout << "====================================================\n";

  std::cout << "Order n:        " << params.n << '\n';
  std::cout << "Edge probability p:   " << params.p << '\n';
  std::cout << "Graphs requested:   " << params.graph_count << '\n';
  std::cout << "Connected tested:   " << result.tested << '\n';
  std::cout << "Disconnected skipped: " << result.skipped << '\n';

  if (result.tested == 0) {
    std::cout << "No connected graphs to evaluate.\n";
    return;
  }

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "Accuracy:       " << result.accuracy * 100.0
        << "%\n";
  std::cout << "Normalized error:   " << result.error << '\n';
  std::cout << "Fail rate:      " << result.fail_rate * 100.0
        << "%\n";
  std::cout << "Average generations:  " << result.avg_generations << '\n';
  std::cout << "Max generations found:" << result.max_generations << '\n';
  std::cout << "Elapsed seconds:    " << elapsed_seconds << '\n';

  std::cout << "\nPer-graph results:\n";
  std::cout << "Graph  Expected  Best  Distance  Gen-found\n";

  for (std::size_t i = 0; i < result.results.size(); ++i) {
    const auto& graph_result = result.results[i];

    std::cout << std::setw(5) << i
          << std::setw(10) << graph_result.expected_z
          << std::setw(6) << graph_result.best_z
          << std::setw(10) << graph_result.distance
          << std::setw(11) << graph_result.generations
          << '\n';
  }
}

int main(int argc, char* argv[]) {
  if (argc != 10) {
    std::cerr
      << "Usage: " << argv[0]
      << " <n> <p> <graph_count> <population_size>"
      << " <gens_before_quit> <max_generations>"
      << " <r_min> <r_max> <seed>\n\n"
      << "Example:\n"
      << argv[0]
      << " 12 0.5 20 30 100 1000 1 10 42\n";
    return EXIT_FAILURE;
  }

  TestParams params{};

  try {
    params.n = std::stoull(argv[1]);
    params.p = std::stod(argv[2]);
    params.graph_count = std::stoull(argv[3]);
    params.population_size = std::stoull(argv[4]);
    params.generations_before_quit = std::stoull(argv[5]);
    params.max_generations = std::stoull(argv[6]);
    params.r_min = std::stoull(argv[7]);
    params.r_max = std::stoull(argv[8]);
    params.seed = std::stoull(argv[9]);
  } catch (const std::exception& e) {
    std::cerr << "Invalid argument: " << e.what() << '\n';
    return EXIT_FAILURE;
  }

  if (params.n == 0 ||
    params.graph_count == 0 ||
    params.population_size == 0 ||
    params.max_generations == 0 ||
    params.r_min > params.r_max ||
    !(params.p >= 0.0 && params.p <= 1.0) ||
    !std::isfinite(params.p)) {
    std::cerr << "Invalid parameters.\n";
    return EXIT_FAILURE;
  }

  std::cout << "Generating random graphs...\n"
        << "n = " << params.n
        << ", p = " << params.p
        << ", count = " << params.graph_count
        << ", seed = " << params.seed << '\n';

  // Assumes GraphGenerator::random(order, count, p) is available.
  // If your generator has a seed parameter, pass params.seed to it.
  std::vector<Graph> graphs =
    GraphGenerator::random(params.n, params.graph_count, params.p);

  std::cout << "Generated " << graphs.size() << " graphs.\n";

  std::cout << "\nTesting every r in ["
        << params.r_min << ", " << params.r_max << "]\n";

  // Inclusive range; the break avoids overflow at SIZE_MAX.
  for (std::size_t r = params.r_min;; ++r) {
    const auto start = std::chrono::steady_clock::now();

    // IMPORTANT:
    // This assumes the configuration constructor's argument is psi.
    // If r is the genetic crossover parameter, set config.r = r instead.
    GeneticSolverConfig config(params.population_size);
    config.r = r;

    std::vector<PerGraphResults> per_graph;
    SuiteResults result;

    double total_error = 0.0;
    double total_generations = 0.0;
    std::size_t num_failed = 0;

    for (const Graph& graph : graphs) {
      if (!GraphMiscInfo::is_connected(graph)) {
        ++result.skipped;
        continue;
      }

      GeneticSolver solver(&graph, config);

      std::size_t generation = 0;
      std::size_t generation_found = 0;

      while (generation < params.max_generations) {
        solver.run(1);
        ++generation;

        const auto info = solver.info();

        if (info.since_better_z() == 0) {
          generation_found = generation;
        }

        if (info.variance() < 0.5 / params.population_size ||
          info.since_better_score() >
            params.generations_before_quit) {
          break;
        }
      }

      const std::size_t expected = zero_forcing_wavefront(graph);
      const std::size_t actual = solver.info().best_z();

      if (actual < expected) {
        std::cerr << "Reference mismatch: best_z=" << actual
              << ", expected=" << expected
              << ", r=" << r << '\n';
        return EXIT_FAILURE;
      }

      const std::size_t distance = actual - expected;

      per_graph.push_back({
        generation_found, distance, expected, actual
      });

      total_error += static_cast<double>(distance) /
               static_cast<double>(graph.order());
      total_generations += static_cast<double>(generation_found);

      if (distance != 0) {
        ++num_failed;
      }
    }

    result.results = std::move(per_graph);
    result.tested = result.results.size();

    if (result.tested != 0) {
      result.error =
        total_error / static_cast<double>(result.tested);
      result.accuracy = 1.0 - result.error;
      result.fail_rate =
        static_cast<double>(num_failed) /
        static_cast<double>(result.tested);
      result.avg_generations =
        total_generations / static_cast<double>(result.tested);

      for (const auto& item : result.results) {
        result.max_generations =
          std::max(result.max_generations, item.generations);
      }
    }

    const auto end = std::chrono::steady_clock::now();
    const double elapsed =
      std::chrono::duration<double>(end - start).count();

    print_results(r, params, result, elapsed);

    if (r == params.r_max) {
      break;
    }
  }

  std::cout << "\nAll r values tested.\n";
  return EXIT_SUCCESS;
}
