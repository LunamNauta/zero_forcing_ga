#ifndef FORT_COVER_HEADER
#define FORT_COVER_HEADER

#include "gurobi_c++.h"

#include "graph.hpp"

const double IP_MAX_TIME = 7200;

struct FortCoverData {
  std::size_t status;
  std::size_t val;
  std::size_t sub_model_count;
  double sub_model_time;
  VertexSet zf_set;
};

class violated_fort_testing_v1: public GRBCallback {
public:
  const Graph *graph;
  GRBVar *s;
  GRBVar *x;    
  GRBModel *model_sub;
  double *sub_model_time;
  std::size_t *sub_model_count;

  violated_fort_testing_v1(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci);

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
  std::size_t *sub_model_count;

  violated_fort_testing_v2(const Graph *gi, GRBVar *si, GRBVar *xi, GRBModel *msi, double *smti, std::size_t *smci);

private:
  void callback();
};

class violated_fort_testing_v3: public GRBCallback{
public:
  const Graph *graph;
  GRBVar *s;
  double *sub_model_time;
  std::size_t *sub_model_count;

  violated_fort_testing_v3(const Graph *gi, GRBVar *si, double *smti, std::size_t *smci);

private:
  void callback();
};

void fort_cover_ip(const Graph &graph, FortCoverData &data, int call_type);

#endif