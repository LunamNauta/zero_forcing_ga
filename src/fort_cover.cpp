#include "fort_cover.hpp"

#include <cassert>
#include <chrono>
#include <limits>
#include <memory>
#include <cmath>

#include "gurobi_c++.h"

#include "genetic_algorithm.hpp"
#include "zero_forcing.hpp"
#include "graph.hpp"

violated_fort_testing_v1::violated_fort_testing_v1(FortCoverSolver *fci, const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci) :
  fc_solver(fci),
  graph(gi),
  s(si),
  x(xi),
  model_sub(msi),
  sub_model_count(smci),
  sub_model_time(smti)
{}

void violated_fort_testing_v1::callback() {
  if (where != GRB_CB_MIPSOL) return;

  std::vector<double> sol(graph->get_order());
  for (Vertex u = 0; u < graph->get_order(); u++) {
    sol[u] = getSolution(s[u]);
  }

  double ub_to_set = -1.0;
  VertexBitset inc_data;
  bool apply_inc = false;

  if (fc_solver) {
    double lower_bound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
    if (lower_bound > 0 && lower_bound - fc_solver->ga_solver->lower_bound() >= 1) {
      std::cout << "Setting lower bound (FC -> GA): " << lower_bound << "\n";
      fc_solver->ga_solver->force_lower_bound(lower_bound);
    }
    
    VertexBitset incumbent(graph->get_order());
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (sol[u] <= 0.5) continue;
      incumbent.set(u);
    }

    bool updated = false;
    {
      std::lock_guard<std::mutex> lock(fc_solver->_update_mutex);
      if (incumbent.count() < fc_solver->_best_incumbent.count()) {
        fc_solver->_best_incumbent = incumbent;
        updated = true;
      }
    }
    if (updated) {
      std::cout << "Setting incumbent solution (FC -> GA): " << incumbent.count() << "\n";
      fc_solver->ga_solver->force_incumbent(incumbent);
    }

    {
      std::lock_guard<std::mutex> lock(fc_solver->_update_mutex);
      
      if (fc_solver->_has_new_incumbent) {
        inc_data = fc_solver->_pending_incumbent;
        apply_inc = true;
        fc_solver->_has_new_incumbent = false;
      }
    }

    if (apply_inc) {
      std::vector<double> vals(graph->get_order());
      for(Vertex u = 0; u < graph->get_order(); u++) {
        vals[u] = inc_data.test(u) ? 1.0 : 0.0;
      }
      setSolution(s, vals.data(), graph->get_order());
      std::cout << "Inserted incumbent solution: " << inc_data.count() << "\n";
    }
  }

  try {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    std::chrono::microseconds duration;
    start = std::chrono::high_resolution_clock::now();

    model_sub->remove(model_sub->getConstrByName("dynamic_constr"));
    model_sub->update();
    VertexSet filled;

    for (Vertex u = 0; u < graph->get_order(); u++){
      if(sol[u] <= 0.5) continue;
      filled.insert(u);
    }

    zero_forcing_closure(*graph, filled);

    GRBLinExpr expr = 0;
    for (Vertex u : filled) {
      expr += x[u];
    }

    model_sub->addConstr(expr, GRB_LESS_EQUAL, 0, "dynamic_constr");
    model_sub->update();
    model_sub->optimize();

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    *sub_model_time += duration.count() * 1E-6;
    *sub_model_count += 1;

    if (model_sub->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      expr = 0;
      for(Vertex u = 0; u < graph->get_order(); u++){
        if (x[u].get(GRB_DoubleAttr_X) <= 0.5) continue;
        expr += s[u];
      }
      addLazy(expr >= 1);
    }
  }
  catch(GRBException &e){
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...){
    std::cout << "Error during violated_fort_testingv1 callback" << "\n";
  }
}

violated_fort_testing_v2::violated_fort_testing_v2(FortCoverSolver *fci, const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci) : 
  fc_solver(fci),
  graph(gi),
  s(si),
  x(xi),
  model_sub(msi),
  sub_model_time(smti),
  sub_model_count(smci)
{}

void violated_fort_testing_v2::callback() {
  if (where != GRB_CB_MIPSOL) return;

  std::vector<double> sol(graph->get_order());
  for (Vertex u = 0; u < graph->get_order(); u++) {
    sol[u] = getSolution(s[u]);
  }

  double ub_to_set = -1.0;
  VertexBitset inc_data;
  bool apply_inc = false;

  if (fc_solver) {
    double lower_bound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
    if (lower_bound > 0 && lower_bound - fc_solver->ga_solver->lower_bound() >= 1) {
      std::cout << "Setting lower bound (FC -> GA): " << lower_bound << "\n";
      fc_solver->ga_solver->force_lower_bound(lower_bound);
    }
    
    VertexBitset incumbent(graph->get_order());
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (sol[u] <= 0.5) continue;
      incumbent.set(u);
    }

    bool updated = false;
    {
      std::lock_guard<std::mutex> lock(fc_solver->_update_mutex);
      if (incumbent.count() < fc_solver->_best_incumbent.count()) {
        fc_solver->_best_incumbent = incumbent;
        updated = true;
      }
    }
    if (updated) {
      std::cout << "Setting incumbent solution (FC -> GA): " << incumbent.count() << "\n";
      fc_solver->ga_solver->force_incumbent(incumbent);
    }

    {
      std::lock_guard<std::mutex> lock(fc_solver->_update_mutex);
      
      if (fc_solver->_has_new_incumbent) {
        inc_data = fc_solver->_pending_incumbent;
        apply_inc = true;
        fc_solver->_has_new_incumbent = false;
      }
    }

    if (apply_inc) {
      std::vector<double> vals(graph->get_order());
      for(Vertex u = 0; u < graph->get_order(); u++) {
        vals[u] = inc_data.test(u) ? 1.0 : 0.0;
      }
      setSolution(s, vals.data(), graph->get_order());
      std::cout << "Inserted incumbent solution: " << inc_data.count() << "\n";
    }
  }

  try {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    std::chrono::microseconds duration;
    start = std::chrono::high_resolution_clock::now();

    model_sub->remove(model_sub->getConstrByName("dynamic_constr"));
    model_sub->update();

    GRBLinExpr expr = 0;
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (getSolution(s[u]) <= 0.5) continue;
      expr += x[u];
    }

    model_sub->addConstr(expr,GRB_LESS_EQUAL,0,"dynamic_constr");
    model_sub->update();
    model_sub->optimize();

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    *sub_model_time += duration.count()*1E-6;
    *sub_model_count += 1;

    if (model_sub->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      expr = 0;
      for (Vertex u = 0; u < graph->get_order(); u++) {
        if(x[u].get(GRB_DoubleAttr_X) <= 0.5) continue;
        expr += s[u];
      }
      addLazy(expr >= 1);
    }
  }
  catch(GRBException &e){
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...){
    std::cout << "Error during violated_fort_testingv2 callback" << "\n";
  } 
}

violated_fort_testing_v3::violated_fort_testing_v3(FortCoverSolver *fci, const Graph *gi, GRBVar *si, double *smti, std::size_t *smci) : 
  fc_solver(fci),
  graph(gi),
  s(si),
  sub_model_time(smti),
  sub_model_count(smci)
{}

void violated_fort_testing_v3::callback() {
  if (where != GRB_CB_MIPSOL) return;

  std::vector<double> sol(graph->get_order());
  for (Vertex u = 0; u < graph->get_order(); u++) {
    sol[u] = getSolution(s[u]);
  }

  double ub_to_set = -1.0;
  VertexBitset inc_data;
  bool apply_inc = false;

  if (fc_solver) {
    double lower_bound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
    if (lower_bound > 0 && lower_bound - fc_solver->ga_solver->lower_bound() >= 1) {
      std::cout << "Setting lower bound (FC -> GA): " << lower_bound << "\n";
      fc_solver->ga_solver->force_lower_bound(lower_bound);
    }
    
    VertexBitset incumbent(graph->get_order());
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if (sol[u] <= 0.5) continue;
      incumbent.set(u);
    }

    bool updated = false;
    {
      std::lock_guard<std::mutex> lock(fc_solver->_update_mutex);
      if (incumbent.count() < fc_solver->_best_incumbent.count()) {
        fc_solver->_best_incumbent = incumbent;
        updated = true;
      }
    }
    if (updated) {
      std::cout << "Setting incumbent solution (FC -> GA): " << incumbent.count() << "\n";
      fc_solver->ga_solver->force_incumbent(incumbent);
    }

    {
      std::lock_guard<std::mutex> lock(fc_solver->_update_mutex);
      
      if (fc_solver->_has_new_incumbent) {
        inc_data = fc_solver->_pending_incumbent;
        apply_inc = true;
        fc_solver->_has_new_incumbent = false;
      }
    }

    if (apply_inc) {
      std::vector<double> vals(graph->get_order());
      for(Vertex u = 0; u < graph->get_order(); u++) {
        vals[u] = inc_data.test(u) ? 1.0 : 0.0;
      }
      setSolution(s, vals.data(), graph->get_order());
      std::cout << "Inserted incumbent solution: " << inc_data.count() << "\n";
    }
  }
  
  try{
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    std::chrono::microseconds duration;
    start = std::chrono::high_resolution_clock::now();

    VertexSet filled;
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if(getSolution(s[u]) <= 0.5) continue;
      filled.insert(u);
    }

    std::size_t pt = zero_forcing_closure(*graph, filled);

    if (pt == std::numeric_limits<std::size_t>::max()) {
      VertexSet vert;
      for (Vertex u = 0; u < graph->get_order(); u++) {
        if (filled.find(u) == filled.cend()) {
          vert.insert(u);
          for (Vertex v : graph->get_adjacent(u)) {
            if (filled.find(v) == filled.cend()) continue;
            vert.insert(v);
          }
        }
      }

      Graph induced = graph->subgraph(vert);

      GRBEnv env_sub(true);
      env_sub.set(GRB_IntParam_OutputFlag, 0);
      env_sub.set(GRB_IntParam_LogToConsole, 0);
      env_sub.start();
      GRBModel model_sub(env_sub);

      std::vector<GRBVar> x(induced.get_order());
      for (Vertex u = 0; u < induced.get_order(); u++) {
        x[u] = model_sub.addVar(0, 1, 1, GRB_BINARY);
      }

      model_sub.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

      GRBLinExpr expr = 0;
      for (Vertex u = 0; u < induced.get_order(); u++) {
        expr += x[u];
      }
      model_sub.addConstr(expr, GRB_GREATER_EQUAL, 1);

      for (Vertex u = 0; u < induced.get_order(); u++) {
        if (filled.find(induced.get_label(u)) != filled.cend()) {
          model_sub.addConstr(x[u], GRB_EQUAL, 0);
        }
      }
        
      for (Vertex u = 0; u < induced.get_order(); u++) {
        for (Vertex v : induced.get_adjacent(u)) {
          expr = x[v] - 2 * x[u];
          for (Vertex w : induced.get_adjacent(v)) {
            expr += x[w];
          }
          model_sub.addConstr(expr, GRB_GREATER_EQUAL, 0);
        }
      }

      model_sub.update();
      model_sub.optimize();

      stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      *sub_model_time += duration.count()*1E-6;
      *sub_model_count += 1;

      if (model_sub.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        expr = 0;
        for (Vertex u = 0; u < induced.get_order(); u++) {
          if (x[u].get(GRB_DoubleAttr_X) <= 0.5) continue;
          expr += s[induced.get_label(u)];
        }
        addLazy(expr >= 1);
      }
    }
  }
  catch(GRBException &e){
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...){
    std::cout << "Error during violated_fort_testingv3 callback" << "\n";
  }
}

FortCoverSolver::FortCoverSolver(const Graph *gi, int type, double t_limit) : 
  graph(gi),
  call_type(type),
  time_limit(t_limit)
{
  _best_incumbent.resize(graph->get_order());
  _best_incumbent.reset(); // Initialized safely to empty
}

void FortCoverSolver::set_ga_solver(GeneticSolver *ga) {
  ga_solver = ga;
}

VertexBitset FortCoverSolver::best_incumbent() const {
  std::lock_guard<std::mutex> lock(_update_mutex);
  return _best_incumbent;
}

void FortCoverSolver::force_incumbent(const VertexBitset& solution) {
  std::lock_guard<std::mutex> lock(_update_mutex);
  _pending_incumbent = solution;
  _has_new_incumbent = true;
}

void FortCoverSolver::run(FortCoverData& data) {
  try {
    // Main model environment
    GRBEnv env_main(true);
    env_main.set(GRB_IntParam_OutputFlag, 1);
    env_main.set(GRB_IntParam_LogToConsole, 1);
    env_main.set(GRB_DoubleParam_TimeLimit, time_limit);
    env_main.set(GRB_IntParam_LazyConstraints, 1);
    env_main.set(GRB_IntParam_Threads, 0);
    env_main.set(GRB_DoubleParam_Heuristics, 0);
    env_main.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);
    env_main.set(GRB_IntParam_Cuts, 2);
    env_main.set(GRB_IntParam_Presolve, 2);
    env_main.start();
    GRBModel model_main(env_main);

    // Main model variables
    std::vector<GRBVar> s(graph->get_order());
    for (Vertex u = 0; u < graph->get_order(); u++) {
      s[u] = model_main.addVar(0, 1, 1, GRB_BINARY);
    }

    // Main model objective
    model_main.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Base constraint ensuring the model has at least 1 row to process
    GRBLinExpr base_expr = 0;
    for (Vertex u = 0; u < graph->get_order(); u++) {
      base_expr += s[u];
    }
    model_main.addConstr(base_expr >= 1, "base_at_least_one");
    model_main.update();

    // Sub model environment
    GRBEnv env_sub(true);
    env_sub.set(GRB_IntParam_OutputFlag, 0);
    env_sub.set(GRB_IntParam_LogToConsole, 0);
    env_sub.start();
    GRBModel model_sub(env_sub);

    // Sub model variables
    std::vector<GRBVar> x(graph->get_order());
    for (Vertex u = 0; u < graph->get_order(); u++) {
      x[u] = model_sub.addVar(0, 1, 1, GRB_BINARY);
    }

    // Sub model objective
    model_sub.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Sub model static constraints
    GRBLinExpr expr = 0;
    for (Vertex u = 0; u < graph->get_order(); ++u) {
      expr += x[u];
    }
    model_sub.addConstr(expr, GRB_GREATER_EQUAL, 1);

    for (Vertex u = 0; u < graph->get_order(); u++) {
      for (Vertex v : graph->get_adjacent(u)) {
        expr = x[v] - 2 * x[u];
        for (Vertex w : graph->get_adjacent(v)) {
          expr += x[w];
        }
        model_sub.addConstr(expr, GRB_GREATER_EQUAL, 0);
      }
    }

    // Sub model dynamic constraint placeholder
    expr = 0;
    model_sub.addConstr(expr, GRB_LESS_EQUAL, 0, "dynamic_constr");
    model_sub.update();

    // Callback setup
    double sub_model_time = 0;
    std::size_t sub_model_count = 0;
    std::unique_ptr<GRBCallback> callback;

    if (call_type == 1) callback = std::make_unique<violated_fort_testing_v1>(this, graph, s.data(), x.data(), &model_sub, &sub_model_time, &sub_model_count);
    else if (call_type == 2) callback = std::make_unique<violated_fort_testing_v2>(this, graph, s.data(), x.data(), &model_sub, &sub_model_time, &sub_model_count);
    else if (call_type == 3) callback = std::make_unique<violated_fort_testing_v3>(this, graph, s.data(), &sub_model_time, &sub_model_count);

    if (callback) {
      model_main.setCallback(callback.get());
      model_main.optimize();
    }

    // Extract results
    data.status = model_main.get(GRB_IntAttr_Status);
    data.sub_model_time = sub_model_time;
    data.sub_model_count = sub_model_count;

    if (data.status == GRB_OPTIMAL) {
      data.val = static_cast<int>(std::round(model_main.get(GRB_DoubleAttr_ObjVal)));
      for (Vertex u = 0; u < graph->get_order(); u++) {
        if (s[u].get(GRB_DoubleAttr_X) <= 0.5) continue;
        data.zf_set.insert(u);
      }
    }
  }
  catch(GRBException &e) {
    std::cerr << "Gurobi Error: " << e.getErrorCode() << ": " << e.getMessage() << std::endl;
  }
  catch(...) {
    std::cerr << "Unknown error during fort cover model execution." << std::endl;
  }
}

void fort_cover_ip(const Graph &graph, FortCoverData &data, const int call_type) {
  try {
    // Main model
    GRBEnv env_main(true);
    env_main.set(GRB_IntParam_OutputFlag, 1);
    env_main.set(GRB_IntParam_LogToConsole, 1);
    env_main.set(GRB_DoubleParam_TimeLimit, IP_MAX_TIME);
    env_main.set(GRB_IntParam_LazyConstraints, 1);
    env_main.set(GRB_IntParam_Threads, 1);
    env_main.set(GRB_DoubleParam_TimeLimit, IP_MAX_TIME);
    env_main.start();
    GRBModel model_main(env_main);

    // Main model variables
    std::vector<GRBVar> s(graph.get_order());
    for (Vertex u = 0; u < graph.get_order(); u++) {
      s[u] = model_main.addVar(0, 1, 1, GRB_BINARY);
    }

    // Main model objective
    model_main.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Base constraint ensuring the model has at least 1 row to process
    GRBLinExpr base_expr = 0;
    for (Vertex u = 0; u < graph.get_order(); u++) {
      base_expr += s[u];
    }
    model_main.addConstr(base_expr >= 1, "base_at_least_one");
    model_main.update();

    // Sub model
    GRBEnv env_sub(true);
    env_sub.set(GRB_IntParam_OutputFlag, 0);
    env_sub.set(GRB_IntParam_LogToConsole, 0);
    env_sub.start();
    GRBModel model_sub(env_sub);

    // Sub model variables
    std::vector<GRBVar> x(graph.get_order());
    for (Vertex u = 0; u < graph.get_order(); u++) {
      x[u] = model_sub.addVar(0, 1, 1, GRB_BINARY);
    }

    // Sub model objective
    model_sub.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Sub model static constraints
    GRBLinExpr expr = 0;
    for (Vertex u = 0; u < graph.get_order(); ++u) {
      expr += x[u];
    }

    model_sub.addConstr(expr, GRB_GREATER_EQUAL, 1);
    for (Vertex u = 0; u < graph.get_order(); u++) {
      for (Vertex v : graph.get_adjacent(u)) {
        expr = x[v] - 2 * x[u];
        for (Vertex w : graph.get_adjacent(v)) {
          expr += x[w];
        }
        model_sub.addConstr(expr, GRB_GREATER_EQUAL, 0);
      }
    }

    // Sub model dynamic constraint placeholder
    expr = 0;
    model_sub.addConstr(expr, GRB_LESS_EQUAL, 0, "dynamic_constr");
    model_sub.update();

    // Callback
    double sub_model_time = 0;
    std::size_t sub_model_count = 0;
    std::unique_ptr<GRBCallback> callback;

    if (call_type == 1) callback = std::make_unique<violated_fort_testing_v1>(nullptr, &graph, s.data(), x.data(), &model_sub, &sub_model_time, &sub_model_count);
    else if (call_type == 2) callback = std::make_unique<violated_fort_testing_v2>(nullptr, &graph, s.data(), x.data(), &model_sub, &sub_model_time, &sub_model_count);
    else if (call_type == 3) callback = std::make_unique<violated_fort_testing_v3>(nullptr, &graph, s.data(), &sub_model_time, &sub_model_count);
    if (callback) {
      model_main.setCallback(callback.get());

      // Optimize
      model_main.optimize();
    }

    // Extract solution only if optimal
    data.status = model_main.get(GRB_IntAttr_Status);
    data.sub_model_time = sub_model_time;
    data.sub_model_count = sub_model_count;

    if (data.status == GRB_OPTIMAL) {
      data.val = std::round(model_main.get(GRB_DoubleAttr_ObjVal));
      for (Vertex u = 0; u < graph.get_order(); u++) {
        if (s[u].get(GRB_DoubleAttr_X) <= 0.5) continue;
        data.zf_set.insert(u);
      }
    }
  }
  catch(GRBException &e){
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...){
    std::cout << "Error during fort cover model" << "\n";
  }
}