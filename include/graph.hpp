#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include <unordered_set>
#include <fstream>
#include <limits>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

typedef std::size_t Vertex;
typedef std::pair<Vertex, Vertex> Edge;
typedef std::unordered_set<Vertex> VertexSet;
typedef boost::dynamic_bitset<> VertexBitset;

constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();
constexpr Vertex INVALID_VERTEX = std::numeric_limits<Vertex>::max();

class Graph {
private:
  std::vector<VertexSet> _adjacent;
  std::vector<Vertex> labels;
  std::size_t _vert_count;
  std::size_t _edge_count;

public:
  Graph(std::size_t order, bool label = true);

  void from_edge_list(const std::ifstream &file);
  void from_edge_list(const std::string &str);

  void from_sparse6(const std::ifstream &file);
  void from_sparse6(const std::string &str);

  void from_graph6(const std::ifstream &file);
  void from_graph6(const std::string &str);

  Graph subgraph(const VertexBitset &vertices);
  Graph subgraph(const VertexSet &vertices);

  VertexBitset adjacent_bitset(Vertex u) const;
  const VertexSet& adjacent(Vertex u) const;
  Vertex label(Vertex u) const;
  std::size_t order() const;
  std::size_t size() const;

  std::size_t degree(Vertex u) const;
  std::size_t max_degree() const;
  std::size_t min_degree() const;

  bool has_edge(Vertex u, Vertex v) const;

  void insert_edge(Vertex u, Vertex v);

  void erase_edge(Vertex u, Vertex v);

  friend class GraphGenerator;
  friend class GraphTreeInfo;
  friend class GraphMiscInfo;
};

class GraphGenerator {
public:
  static std::vector<Graph> random(std::size_t order, std::size_t count, double probability);
  static std::vector<Graph> cubic(std::size_t order, std::size_t count);
  static Graph path(std::size_t order);
  static Graph cycle(std::size_t order);
  static Graph complete(std::size_t order);
};

class GraphTreeInfo {
public:
  static bool is_forest(const Graph &graph);
  static bool is_tree(const Graph &graph);

  static std::size_t diameter(const Graph &graph);

  static std::pair<Vertex, Vertex> center(const Graph &graph);
};

class GraphMiscInfo {
public:
  static bool is_connected(const Graph &graph);

  static VertexBitset pendants_bs(const Graph &graph);
  static VertexSet pendants(const Graph &graph);

  static VertexBitset complement_bs(const Graph &graph, const VertexBitset &vertices);
  static VertexSet complement(const Graph &graph, const VertexSet &vertices);
};

#endif
