#include <limits>
#include <thread>

#include "genetic_algorithm.hpp"
#include "fort_cover.hpp"

int main() {
  Graph graph = Graph::generate_random(20, 1, 0.5).front();

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