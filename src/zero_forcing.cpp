#include "zero_forcing.hpp"

#include <limits>
#include <vector>
#include <limits>
#include <list>

#include "graph.hpp"

#include "gurobi_c++.h"

std::size_t zero_forcing_closure(const Graph &graph, VertexBitset &closure) {
  std::vector<std::size_t> uncolored_degree(graph.get_order(), 0);

  std::vector<Vertex> current_ready; 
  std::vector<Vertex> next_ready;
  current_ready.reserve(graph.get_order());
  next_ready.reserve(graph.get_order());

  for (Vertex u = 0; u < graph.get_order(); u++) {
    for (Vertex neighbor : graph.get_neighbors(u)) uncolored_degree[u] += !closure.test(neighbor);
    if (closure.test(u) && uncolored_degree[u] == 1) current_ready.push_back(u);
  }

  std::size_t pt = 0;
  while (!current_ready.empty()) {
    while (!current_ready.empty()) {
      Vertex u = current_ready.back();
      current_ready.pop_back();
      if (uncolored_degree[u] != 1) continue;

      Vertex v = graph.get_order(); 
      for (Vertex vi : graph.get_neighbors(u)) {
        if (closure.test(vi)) continue;
        v = vi;
        break;
      }
      if (v == graph.get_order()) continue; 

      closure.set(v);
      for (Vertex w : graph.get_neighbors(v)) {
        uncolored_degree[w]--;      
        if (closure.test(w) && uncolored_degree[w] == 1) next_ready.push_back(w);
      }
          
      if (uncolored_degree[v] == 1) next_ready.push_back(v);
    }

    pt++;
    
    current_ready = std::move(next_ready);
    next_ready.clear();
  }

  if (closure.count() != graph.get_order()) return std::numeric_limits<std::size_t>::max();
  return pt;
}

std::size_t zero_forcing_closure(const Graph &graph, VertexSet &filled) {
  VertexSet next;
  std::size_t pt;

  for (pt = 0; pt < graph.get_order(); pt++){
    next.clear();

    for (Vertex u : filled) {
      std::size_t white_count = 0;
      Vertex forced_vertex;

      for (Vertex v : graph.get_adjacent(u)) {
        if (filled.find(v) != filled.end()) continue;
        if (++white_count > 1) break;
        forced_vertex = v;
      }
      if (white_count == 1) next.insert(forced_vertex);
    }
    if (next.empty()) break;
    
    filled.insert(next.begin(), next.end());
  }

  if (filled.size() != graph.get_order()) return std::numeric_limits<std::size_t>::max();
  return pt;
}

std::size_t zero_forcing_wavefront(const Graph &graph, std::size_t upper_bound) {
  if (graph.get_order() == 0) return 0;
  upper_bound = std::min(graph.get_order(), upper_bound);

  std::list<std::pair<VertexBitset, std::size_t>> cl_pairs;
  cl_pairs.emplace_back(VertexBitset(graph.get_order(), false), 0);

  for (std::size_t R = 1; R <= upper_bound - 1; R++) {
    std::list<std::pair<VertexBitset, std::size_t>> next_cl_pairs;

    // Use structured binding for cleaner access
    for (const auto& [S, r] : cl_pairs) {
      for (Vertex v = 0; v < graph.get_order(); v++) {
        std::size_t r_new = r;
        if (!S.test(v)) r_new++;

        const VertexSet& neighbors = graph.get_adjacent(v);
        std::size_t neighbors_outside_S = 0;
        
        for (Vertex neighbor : neighbors) {
          if (!S.test(neighbor)) neighbors_outside_S++;
        }
        
        if (neighbors_outside_S > 1) r_new += (neighbors_outside_S - 1);

        // Prune paths early
        if (r_new != R) continue;

        VertexBitset S_new = S; 
        S_new.set(v);
        for (Vertex u : neighbors) S_new.set(u);

        zero_forcing_closure(graph, S_new);

        if (S_new.count() == graph.get_order()) return r_new;

        bool already_present = false;
        
        // Check the main list for duplicates or better paths
        for (const auto& existing : cl_pairs) {
          if (existing.first == S_new && existing.second <= r_new) {
            already_present = true;
            break;
          }
        }

        // Check the pending generation list to prevent duplicate branching
        if (!already_present) {
          for (const auto& pending : next_cl_pairs) {
            if (pending.first == S_new && pending.second <= r_new) {
              already_present = true;
              break;
            }
          }
        }

        if (!already_present) {
          next_cl_pairs.emplace_back(std::move(S_new), r_new);
        }
      }
    }
    cl_pairs.splice(cl_pairs.end(), next_cl_pairs);
  }

  return upper_bound;
}

std::size_t zero_forcing_wavefront_old(const Graph &graph) {
  if (graph.get_order() == 0) return 0;

  std::list<std::pair<VertexSet, std::size_t>> cl_pairs;
  cl_pairs.push_back(std::pair<VertexSet, std::size_t>({}, 0));

  for (std::size_t R = 1; R <= graph.get_order(); R++) {
    for (auto j = cl_pairs.cbegin(); j != cl_pairs.cend(); j++) {
      const VertexSet &S = j->first;
      std::size_t r = j->second;

      for (Vertex v = 0; v < graph.get_order(); v++) {
        VertexSet S_new(S.cbegin(), S.cend());

        S_new.insert(v);

        VertexSet neighbors = graph.get_adjacent(v);
        S_new.insert(neighbors.cbegin(), neighbors.cend());

        zero_forcing_closure(graph,S_new);

        std::size_t r_new = r;
        if (S.find(v) == S.cend()) r_new++;

        std::size_t neighbors_outside_S = 0;
        for (VertexSet::const_iterator it = neighbors.cbegin(); it != neighbors.cend(); it++) {
          if (S.find(*it) != S.cend()) continue;

          neighbors_outside_S++;
        }

        if (neighbors_outside_S > 1) r_new += (neighbors_outside_S - 1);

        if (r_new <= R) {
          bool already_present = false;

          for (auto it = cl_pairs.cbegin(); it != cl_pairs.cend(); it++) {
            if (it->first != S_new || it->second > r_new) continue;

            already_present = true;
            break;
          }

          if (!already_present) {
            cl_pairs.push_back(std::pair<VertexSet, int>(S_new, r_new));

            if (S_new.size() == graph.get_order()) return r_new;
          }
        }
      }
    }
  }

  return graph.get_order();
}

static GRBEnv& get_gurobi_base_env() {
  static GRBEnv env([]() {
    GRBEnv empty_env(true); // Create empty environment
    empty_env.set(GRB_IntParam_OutputFlag, 0);
    empty_env.set(GRB_IntParam_LogToConsole, 0);
    empty_env.set(GRB_DoubleParam_Heuristics, 0);
    empty_env.set(GRB_IntParam_Threads, 0);
    empty_env.set(GRB_IntParam_Cuts, 0);
    // empty_env.set(GRB_IntParam_Presolve, 0);
    empty_env.start();
    return empty_env;
  }());
  return env;
}
VertexBitset minimum_fort_ip_subgraph(const Graph &graph, const VertexBitset &filled) {
  VertexBitset closure(filled);
  std::size_t pt = zero_forcing_closure(graph, closure);

  VertexSet vertices;
  for (Vertex u = 0; u < graph.get_order(); u++) {
    if (closure.test(u)) continue;
    vertices.insert(u);
    for (Vertex v : graph.get_neighbors(u)) {
      if (!closure.test(v)) continue;
      vertices.insert(v);
    }
  }

  Graph induced = graph.subgraph(vertices);
  std::size_t order = induced.get_order();

  // Fetch the globally unique environment
  GRBEnv &env = get_gurobi_base_env();
  GRBModel model(env);
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

  std::vector<GRBVar> x;
  x.reserve(order);
  for (Vertex u = 0; u < order; u++) {
    double ub = closure.test(induced.get_label(u)) ? 0.0 : 1.0;
    x.push_back(model.addVar(0.0, ub, 1.0, GRB_BINARY));
  }

  GRBLinExpr expr;
  std::vector<double> ones(order, 1.0);
  expr.addTerms(ones.data(), x.data(), order);
  model.addConstr(expr, GRB_GREATER_EQUAL, 1.0);

  std::vector<double> coeffs;
  std::vector<GRBVar> vars;
  coeffs.reserve(order + 2);
  vars.reserve(order + 2);

  for (Vertex u = 0; u < order; u++) {
    for (Vertex v : induced.get_neighbors(u)) {
      coeffs.clear(); 
      vars.clear();

      // Add x[v]
      coeffs.push_back(1.0); 
      vars.push_back(x[v]);

      // Add -2 * x[u]
      coeffs.push_back(-2.0); 
      vars.push_back(x[u]);

      // Add sum(x[w])
      for (Vertex w : induced.get_neighbors(v)) {
        coeffs.push_back(1.0); 
        vars.push_back(x[w]);
      }

      expr.clear();
      expr.addTerms(coeffs.data(), vars.data(), coeffs.size());
      model.addConstr(expr, GRB_GREATER_EQUAL, 0.0);
    }
  }

  model.optimize();

  VertexBitset sol(graph.get_order());
  if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    for (Vertex u = 0; u < order; u++) {
      if (x[u].get(GRB_DoubleAttr_X) > 0.5) {
        sol.set(induced.get_label(u));
      }
    }
  }

  return sol;
}