#ifndef FORT_COVER_HEADER
#define FORT_COVER_HEADER

#include <mutex>

#include "gurobi_c++.h"

#include "graph.hpp"

class GeneticSolver;

const double IP_MAX_TIME = 7200;

struct FortCoverData {
  std::size_t status;
  std::size_t val;
  std::size_t sub_model_count;
  double sub_model_time;
  VertexSet zf_set;
};

class FortCoverSolver;

class violated_fort_testing_v1: public GRBCallback {
public:
  FortCoverSolver *fc_solver;
  const Graph *graph;
  GRBVar *s;
  GRBVar *x;    
  GRBModel *model_sub;
  double *sub_model_time;
  std::size_t *sub_model_count;
    int cuts_added = 0;

  violated_fort_testing_v1(FortCoverSolver *fc, const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci);

private:
  void callback();
};

class violated_fort_testing_v2: public GRBCallback {
public:
  FortCoverSolver *fc_solver;
  const Graph *graph;
  GRBVar *s;
  GRBVar *x;
  GRBModel *model_sub;
  double *sub_model_time;
  std::size_t *sub_model_count;

  violated_fort_testing_v2(FortCoverSolver *fc, const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci);

private:
  void callback();
};

class violated_fort_testing_v3: public GRBCallback{
public:
  FortCoverSolver *fc_solver;
  const Graph *graph;
  GRBVar *s;
  double *sub_model_time;
  std::size_t *sub_model_count;

  violated_fort_testing_v3(FortCoverSolver *fc, const Graph *gi, GRBVar *si, double *smti, std::size_t *smci);

private:
  void callback();
};

class FortCoverSolver {
private:
  mutable std::mutex _update_mutex;

  VertexBitset _best_incumbent;
  VertexBitset _pending_incumbent;
  bool _has_new_incumbent = false;

  GeneticSolver *ga_solver;
  const Graph *graph;
  int call_type;
  double time_limit;

public:
  FortCoverSolver(const Graph *gi, int type, double t_limit = 60.0);

  void set_ga_solver(GeneticSolver *ga);

  VertexBitset best_incumbent() const;
  void force_incumbent(const VertexBitset& solution);
  void force_lower_bound(double lower_bound);

  void run(FortCoverData& data);

  friend class violated_fort_testing_v1;
  friend class violated_fort_testing_v2;
  friend class violated_fort_testing_v3;
};

void fort_cover_ip(const Graph &graph, FortCoverData &data, int call_type);

#endif