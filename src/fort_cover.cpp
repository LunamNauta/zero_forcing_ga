#include "fort_cover.hpp"

#include <cassert>
#include <chrono>
#include <limits>
#include <memory>
#include <cmath>

#include "gurobi_c++.h"

#include "zero_forcing.hpp"
#include "graph.hpp"

violated_fort_testing_v1::violated_fort_testing_v1(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci) :
  graph(gi),
  s(si),
  x(xi),
  model_sub(msi),
  sub_model_count(smci),
  sub_model_time(smti)
{}

void violated_fort_testing_v1::callback() {
  if (where != GRB_CB_MIPSOL) return;
  try {
    // Found an integer solution, use closure to update sub model to see if it violates a fort constraint
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    std::chrono::microseconds duration;
    start = std::chrono::high_resolution_clock::now();

    model_sub->remove(model_sub->getConstrByName("dynamic_constr"));
    model_sub->update();
    VertexSet filled;

    for (Vertex u = 0; u < graph->get_order(); u++){
      if(getSolution(s[u]) <= 0.5) continue;
      filled.insert(u);
    }

    zero_forcing_closure(*graph, filled);

    GRBLinExpr expr = 0;
    for (Vertex u : filled) {
      expr += x[u];
    }

    model_sub->addConstr(expr,GRB_LESS_EQUAL,0,"dynamic_constr");
    model_sub->update();
    model_sub->optimize();

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    *sub_model_time += duration.count()*1E-6;
    *sub_model_count += 1;

    // If sub model is feasible then there is a minimum fort violation, which we add back to main model
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

violated_fort_testing_v2::violated_fort_testing_v2(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci) : 
  graph(gi),
  s(si),
  x(xi),
  model_sub(msi),
  sub_model_time(smti),
  sub_model_count(smci)
{}

void violated_fort_testing_v2::callback() {
  if (where != GRB_CB_MIPSOL) return;
  try {
    // Found an integer solution, use initial set to update sub model to see if it violates a fort constraint
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

    // If sub model is feasible then there is a minimum fort violation, which we add back to main model
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

violated_fort_testing_v3::violated_fort_testing_v3(const Graph *gi, GRBVar *si, double *smti, std::size_t *smci) : 
  graph(gi),
  s(si),
  sub_model_time(smti),
  sub_model_count(smci)
{}

void violated_fort_testing_v3::callback() {
  if (where != GRB_CB_MIPSOL) return;
  try{
    // Found an integer solution, use closure to construct sub model to see if it violates a fort constraint
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    std::chrono::microseconds duration;
    start = std::chrono::high_resolution_clock::now();

    VertexSet filled;
    for (Vertex u = 0; u < graph->get_order(); u++) {
      if(getSolution(s[u]) <= 0.5) continue;
      filled.insert(u);
    }

    std::size_t pt = zero_forcing_closure(*graph, filled);

    // If filled is not entire vertex set, add fort constraint to main model
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

      // If sub model is feasible then there is a minimum fort violation, which we add back to main model
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

    if (call_type == 1) callback = std::make_unique<violated_fort_testing_v1>(&graph, s.data(), x.data(), &model_sub, &sub_model_time, &sub_model_count);
    else if (call_type == 2) callback = std::make_unique<violated_fort_testing_v2>(&graph, s.data(), x.data(), &model_sub, &sub_model_time, &sub_model_count);
    else if (call_type == 3) callback = std::make_unique<violated_fort_testing_v3>(&graph, s.data(), &sub_model_time, &sub_model_count);
    if (callback) {
      model_main.setCallback(callback.get());

      // Optimize
      model_main.optimize();
    };

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