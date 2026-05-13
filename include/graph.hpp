#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include <unordered_set>
#include <iosfwd>
#include <limits>
#include <vector>
#include <string>

typedef std::size_t VertexIndex;
typedef std::unordered_set<VertexIndex> AdjacencySet;
typedef std::unordered_set<VertexIndex> VertexSet;
typedef std::pair<VertexIndex, VertexIndex> EdgePair;

inline constexpr VertexIndex INVALID_VERTEX = std::numeric_limits<VertexIndex>::max();
inline constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

class Graph { 
private:
  std::size_t vert_count;
  std::size_t edge_count;
  std::vector<VertexIndex> labels;
  std::vector<AdjacencySet> adj;

public:
  // Initialization ----------------------------------------------------------------
  Graph();
  Graph(std::size_t ord, bool label = true);

  void clear();
  void clear(std::size_t ord, bool label = true);

  void from_graph6(const std::string &str);
  void from_sparse6(const std::string &str);
  void from_edge_list(const std::string &str);

  void from_graph6(const std::ifstream &file);  
  void from_sparse6(const std::ifstream &file);
  void from_edge_list(const std::ifstream &file);

  Graph subgraph(const std::unordered_set<VertexIndex> &vertices) const;

  // Element Access ----------------------------------------------------------------
  std::size_t get_order() const;
  std::size_t get_vertex_count() const;
  std::size_t get_edge_count() const;

  std::size_t get_degree(VertexIndex vert) const;
  std::size_t get_max_degree() const;
  VertexIndex get_label(VertexIndex vert) const;
  AdjacencySet get_adjacent(VertexIndex vert) const;

  // Query ----------------------------------------------------------------
  bool has_edge(VertexIndex vert1, VertexIndex vert2);

  std::vector<VertexIndex> get_vertices() const;
  std::vector<EdgePair> get_edges() const;

  // Insertion ----------------------------------------------------------------
  void insert_edge(VertexIndex vert1, VertexIndex vert2);
  void insert_vertex(std::size_t count = 1);

  // Erasure ----------------------------------------------------------------
  void erase_edge(VertexIndex vert1, VertexIndex vert2);
  void erase_vertex(VertexIndex vert);

  // Tree ----------------------------------------------------------------
  bool is_valid_forest() const;
  bool is_valid_tree() const;
  std::size_t tree_diameter() const;
  std::pair<VertexIndex, VertexIndex> tree_center() const;

  // Zero Forcing ----------------------------------------------------------------
  bool is_valid_zf(const VertexSet &filled) const;
  VertexSet zf_closure(VertexSet filled, std::size_t *pt = nullptr) const;
  std::size_t zf_wavefront() const;

  // Misc ----------------------------------------------------------------
  VertexSet get_pendants() const;
  bool is_connected() const;

  // Output ----------------------------------------------------------------
  friend std::ostream& operator<<(std::ostream &os, const Graph &graph);
};

#endif