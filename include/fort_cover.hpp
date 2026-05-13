#ifndef FORT_COVER_HEADER
#define FORT_COVER_HEADER

#include "gurobi_c++.h"

#include "graph.hpp"

const double IP_MAX_TIME = 7200;

struct fort_cover_data {
  int status;
  int val;
  double sub_model_time;
  std::unordered_set<VertexIndex> zf_set;
};

class violated_fort_testing_v1: public GRBCallback {
public:
  const Graph *graph;
  GRBVar *s;
  GRBVar *x;
  GRBModel *model_sub;
  double *sub_model_time;

  violated_fort_testing_v1(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *mi, double *smi);

private:
  void callback();
};

class violated_fort_testing_v2: public GRBCallback {
public:
  const Graph *graph;
  GRBVar *s;
  GRBVar *x;
  GRBModel *model_sub;
  double *sub_model_time;

  violated_fort_testing_v2(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *mi, double *smi);

private:
  void callback();
};

class violated_fort_testing_v3: public GRBCallback {
public:
  const Graph *graph;
  GRBVar *s;
  double *sub_model_time;

  violated_fort_testing_v3(const Graph *gi, GRBVar *si, double *smi);

private:
  void callback();
};

void fort_cover_ip(const Graph &graph, fort_cover_data &data, int call_type);

#endif