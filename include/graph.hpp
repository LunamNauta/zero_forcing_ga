#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include <unordered_set>
#include <iosfwd>
#include <limits>
#include <vector>
#include <string>

typedef std::size_t Vertex;
typedef std::unordered_set<Vertex> VertexSet;
typedef std::pair<Vertex, Vertex> Edge;

inline constexpr Vertex INVALID_VERTEX = std::numeric_limits<Vertex>::max();
inline constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

class Graph { 
private:
  std::size_t vert_count;
  std::size_t edge_count;
  std::vector<Vertex> labels;
  std::vector<VertexSet> adj;

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

  Graph subgraph(const std::unordered_set<Vertex> &vertices) const;

  // Element Access ----------------------------------------------------------------
  std::size_t get_order() const;
  std::size_t get_vertex_count() const;
  std::size_t get_edge_count() const;

  std::size_t get_degree(Vertex u) const;
  std::size_t get_max_degree() const;
  Vertex get_label(Vertex u) const;
  VertexSet get_adjacent(Vertex u) const;

  // Query ----------------------------------------------------------------
  bool has_edge(Vertex u, Vertex v);

  std::vector<Vertex> get_vertices() const;
  std::vector<Edge> get_edges() const;

  // Insertion ----------------------------------------------------------------
  void insert_edge(Vertex u, Vertex v);
  void insert_vertex(std::size_t count = 1);

  // Erasure ----------------------------------------------------------------
  void erase_edge(Vertex u, Vertex v);
  void erase_vertex(Vertex u);

  // Tree ----------------------------------------------------------------
  bool is_valid_forest() const;
  bool is_valid_tree() const;
  std::size_t tree_diameter() const;
  std::pair<Vertex, Vertex> tree_center() const;

  // Misc ----------------------------------------------------------------
  VertexSet get_pendants() const;
  VertexSet get_complement(const VertexSet &verts) const;
  bool is_connected() const;

  // Output ----------------------------------------------------------------
  friend std::ostream& operator<<(std::ostream &os, const Graph &graph);
};

#endif