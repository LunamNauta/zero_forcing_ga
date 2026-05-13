#include "fort_cover.hpp"

#include <chrono>
#include <cmath>

#include "zero_forcing.hpp"
#include "graph.hpp"

violated_fort_testing_v1::violated_fort_testing_v1(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *mi, double *smi) : 
  graph(gi), s(si), x(xi), model_sub(mi), sub_model_time(smi)
{}

void violated_fort_testing_v1::callback() {
  try {
    if (where == GRB_CB_MIPSOL) {
      // Found an integer solution, use closure to update sub model to see if it violates a fort constraint
      std::chrono::time_point<std::chrono::high_resolution_clock> start;
      std::chrono::time_point<std::chrono::high_resolution_clock> stop;
      std::chrono::microseconds duration;

      start = std::chrono::high_resolution_clock::now();
      model_sub->remove(model_sub->getConstrByName("dynamic_constr"));
      model_sub->update();

      VertexSet filled;
      for (std::size_t a = 0; a < graph->get_order(); a++){
        if (getSolution(s[a]) <= 0.5) continue;
        filled.insert(a);
      }

      filled = zero_forcing_closure(*graph, filled);
        
      GRBLinExpr expr = 0;
      for (VertexSet::const_iterator it = filled.cbegin(); it != filled.cend(); it++){
        expr += x[*it];
      }

      model_sub->addConstr(expr, GRB_LESS_EQUAL, 0, "dynamic_constr");
      model_sub->update();
      model_sub->optimize();

      stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      *sub_model_time += duration.count()*1E-6;

      // If sub model is feasible then there is a minimum fort violation, which we add back to main model
      if (model_sub->get(GRB_IntAttr_Status) != GRB_OPTIMAL) return;

      expr = 0;
      for(std::size_t a = 0; a < graph->get_order(); a++) {
        if (x[a].get(GRB_DoubleAttr_X) <= 0.5) continue;
        expr += s[a];
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

violated_fort_testing_v2::violated_fort_testing_v2(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *mi, double *smi) :
  graph(gi), s(si), x(xi), model_sub(mi), sub_model_time(smi) {}

void violated_fort_testing_v2::callback() {
  try {
    if (where == GRB_CB_MIPSOL) {
      // Found an integer solution, use initial set to update sub model to see if it violates a fort constraint
      std::chrono::time_point<std::chrono::high_resolution_clock> start;
      std::chrono::time_point<std::chrono::high_resolution_clock> stop;
      std::chrono::microseconds duration;

      start = std::chrono::high_resolution_clock::now();
      model_sub->remove(model_sub->getConstrByName("dynamic_constr"));
      model_sub->update();

      GRBLinExpr expr = 0;
      for (std::size_t a = 0; a < graph->get_order(); a++){
        if (getSolution(s[a]) <= 0.5) continue;
        expr += x[a];
      }

      model_sub->addConstr(expr, GRB_LESS_EQUAL, 0, "dynamic_constr");
      model_sub->update();
      model_sub->optimize();

      stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      *sub_model_time += duration.count()*1E-6;

      // If sub model is feasible then there is a minimum fort violation, which we add back to main model
      if (model_sub->get(GRB_IntAttr_Status) != GRB_OPTIMAL) return;

      expr = 0;
      for (std::size_t a = 0; a < graph->get_order(); a++) {
        if (x[a].get(GRB_DoubleAttr_X) <= 0.5) continue;
        expr += s[a];
      }
      addLazy(expr >= 1);
      }
  }
  catch(GRBException &e) {
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...){
    std::cout << "Error during violated_fort_testingv2 callback" << "\n";
  }
}

violated_fort_testing_v3::violated_fort_testing_v3(const Graph *gi, GRBVar *si, double *smi) :
  graph(gi), s(si), sub_model_time(smi) {}

void violated_fort_testing_v3::callback() {
  try {
    if (where == GRB_CB_MIPSOL){
      // Found an integer solution, use closure to construct sub model to see if it violates a fort constraint
      std::chrono::time_point<std::chrono::high_resolution_clock> start;
      std::chrono::time_point<std::chrono::high_resolution_clock> stop;
      std::chrono::microseconds duration;
      start = std::chrono::high_resolution_clock::now();

      VertexSet filled;
      for (std::size_t a = 0; a < graph->get_order(); a++) {
        if (getSolution(s[a]) <= 0.5) continue;
        filled.insert(a);
      }

      filled = zero_forcing_closure(*graph, filled);

      // If filled is not entire vertex set, add fort constraint to main model
      if (filled.size() >= graph->get_order()) return;

      VertexSet verts;
      for (std::size_t a = 0; a < graph->get_order(); a++) {
        if (filled.find(a) != filled.cend()) continue;
        verts.insert(a);

        VertexSet neighbors = graph->get_adjacent(a);
        for (VertexSet::const_iterator v = neighbors.cbegin(); v != neighbors.cend(); ++v) {
          if (filled.find(*v) == filled.cend()) continue;
          verts.insert(*v);
        }
      }

      Graph induced = graph->subgraph(verts);

      GRBEnv env_sub(true);
      env_sub.set(GRB_IntParam_OutputFlag, 0);
      env_sub.set(GRB_IntParam_LogToConsole, 0);
      env_sub.start();
      GRBModel model_sub(env_sub);

      std::vector<GRBVar> x(induced.get_order());
      for (std::size_t a = 0; a < induced.get_order(); a++) {
        x[a] = model_sub.addVar(0, 1, 1, GRB_BINARY);
      }

      model_sub.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

      GRBLinExpr expr = 0;
      for (std::size_t a = 0; a < induced.get_order(); a++) {
        expr += x[a];
      }
      model_sub.addConstr(expr, GRB_GREATER_EQUAL, 1);

      for (std::size_t a = 0; a < induced.get_order(); a++) {
        if (filled.find(induced.get_label(a)) == filled.cend()) continue;
        model_sub.addConstr(x[a], GRB_EQUAL, 0);
      }

      for (std::size_t a = 0; a < induced.get_order(); a++) {
        VertexSet neighbors1 = induced.get_adjacent(a);
        for (VertexSet::const_iterator it1 = neighbors1.cbegin(); it1 != neighbors1.cend(); it1++) {
          VertexSet neighbors2 = induced.get_adjacent(*it1);
          expr = x[*it1] - 2*x[a];
          for (VertexSet::const_iterator it2 = neighbors2.cbegin(); it2 != neighbors2.cend(); it2++) {
            expr += x[*it2];
          }
          model_sub.addConstr(expr, GRB_GREATER_EQUAL, 0);
        }
      }

      model_sub.update();
      model_sub.optimize();

      stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      *sub_model_time += duration.count()*1E-6;

      // If sub model is feasible then there is a minimum fort violation, which we add back to main model
      if (model_sub.get(GRB_IntAttr_Status) != GRB_OPTIMAL) return;

      expr = 0;
      for (std::size_t a = 0; a < induced.get_order(); a++) {
        if (x[a].get(GRB_DoubleAttr_X) <= 0.5) continue;
        expr += s[induced.get_label(a)];
      }
      addLazy(expr >= 1);
    }
  }
  catch(GRBException &e) {
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...) {
    std::cout << "Error during violated_fort_testingv3 callback" << "\n";
  }
}

void fort_cover_ip(const Graph &graph, fort_cover_data &data, int call_type) {
  try {
    // Main model ----------------------------------------------------------------
    GRBEnv env_main(true);
    env_main.set(GRB_IntParam_OutputFlag, 1);
    env_main.set(GRB_IntParam_LogToConsole, 1);
    env_main.set(GRB_IntParam_LazyConstraints, 1);
    env_main.set(GRB_IntParam_Threads, 1);
    env_main.set(GRB_DoubleParam_TimeLimit, IP_MAX_TIME);
    env_main.start();
    GRBModel model_main(env_main);
    // Main model variables
    std::vector<GRBVar> s(graph.get_order());
    for (std::size_t a = 0; a < graph.get_order(); a++) {
      s[a] = model_main.addVar(0, 1, 1, GRB_BINARY);
    }
    // Main model objective
    model_main.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Sub Model ----------------------------------------------------------------
    GRBEnv env_sub(true);
    env_sub.set(GRB_IntParam_OutputFlag, 0);
    env_sub.set(GRB_IntParam_LogToConsole, 0);
    env_sub.start();
    GRBModel model_sub(env_sub);

    // Sub model Variables
    std::vector<GRBVar> x(graph.get_order());
    for (std::size_t a = 0; a < graph.get_order(); a++) {
      x[a] = model_sub.addVar(0, 1, 1, GRB_BINARY);
    }

    // Sub model objective
    model_sub.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Sub model static constraints
    GRBLinExpr expr = 0;
    for (std::size_t a = 0; a < graph.get_order(); a++) {
      expr += x[a];
    }
    model_sub.addConstr(expr, GRB_GREATER_EQUAL, 1);

    for (std::size_t a = 0; a < graph.get_order(); a++) {
      VertexSet neighbors1 = graph.get_adjacent(a);
      for (VertexSet::const_iterator it1 = neighbors1.cbegin(); it1 != neighbors1.cend(); it1++) {
        expr = x[*it1] - 2 * x[a];

        VertexSet neighbors2 = graph.get_adjacent(*it1);
        for (VertexSet::const_iterator it2 = neighbors2.cbegin(); it2 != neighbors2.cend(); it2++) {
          expr += x[*it2];
        }

        model_sub.addConstr(expr, GRB_GREATER_EQUAL, 0);
      }
    }

    // Sub model dynamic constraint placeholder
    expr = 0;
    model_sub.addConstr(expr, GRB_LESS_EQUAL, 0, "dynamic_constr");
    model_sub.update();

    // Callback ----------------------------------------------------------------
    double sub_model_time = 0;
    GRBCallback* callback = nullptr;

    if (call_type == 1) {
      callback = new violated_fort_testing_v1(&graph, s.data(), x.data(), &model_sub, &sub_model_time);
    }
    else if (call_type == 2) {
      callback = new violated_fort_testing_v2(&graph, s.data(), x.data(), &model_sub, &sub_model_time);
    }
    else if (call_type == 3) {
      callback = new violated_fort_testing_v3(&graph, s.data(), &sub_model_time);
    }

    if (callback != nullptr) {
      model_main.setCallback(callback);
    }

    // Optimize
    model_main.optimize();

    // Extract solution only if optimal
    data.status = model_main.get(GRB_IntAttr_Status);
    data.zf_set = VertexSet();
    data.val = -1;

    if (data.status == GRB_OPTIMAL) {
      data.val = std::round(model_main.get(GRB_DoubleAttr_ObjVal));

      for (std::size_t a = 0; a < graph.get_order(); a++) {
        if (s[a].get(GRB_DoubleAttr_X) <= 0.5) continue;
        data.zf_set.insert(a);
      }
    }

    data.sub_model_time = sub_model_time;
  }
  catch(GRBException &e){
    std::cout << "Error number: " << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  }
  catch(...){
    std::cout << "Error during fort cover model" << "\n";
  }
}