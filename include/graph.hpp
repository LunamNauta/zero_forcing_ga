#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include <unordered_set>
#include <iosfwd>
#include <limits>
#include <vector>
#include <string>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

/** @brief Type alias for vertex indices. */
typedef std::size_t Vertex;

/** @brief Collection of unique vertices. */
typedef std::unordered_set<Vertex> VertexSet;

/** @brief Efficient bit representation of a set of vertices. */
typedef boost::dynamic_bitset<> VertexBitset;

/** @brief Represents an undirected edge as a pair of vertices. */
typedef std::pair<Vertex, Vertex> Edge;

/** @brief Sentinel value for an uninitialized or non-existent vertex. */
inline constexpr Vertex INVALID_VERTEX = std::numeric_limits<Vertex>::max();

/** @brief Sentinel value for an invalid index or size. */
inline constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

/**
 * @class Graph
 * @brief An undirected graph representation using an adjacency set.
 * 
 * Provides utilities for graph construction, property querying, 
 * and and standard serialization formats (graph6/sparse6/edge lists).
 */
class Graph { 
private:
  std::size_t vert_count;     //!< The number of vertices (order) in the graph.
  std::size_t edge_count;     //!< The number of edges (size) in the graph.
  std::vector<Vertex> labels; //!< Mapping of internal indices to external labels.
  std::vector<VertexSet> adj; //!< Adjacency list using sets.

public:
  // Initialization ----------------------------------------------------------------
  /** @brief Constructs an empty graph. */
  Graph();

  /** 
   * @brief Constructs a graph with a fixed number of vertices. 
   * @param ord Initial number of vertices.
   * @param label Whether to initialize default labels.
   */
  Graph(std::size_t ord, bool label = true);

  /** @brief Resets the graph, removing all vertices and edges. */
  void clear();

  /** @brief Clears the graph and re-initializes with a specific order. */
  void clear(std::size_t ord, bool label = true);

  /** @name Serialization
   *  Methods for loading graphs from standard compressed formats.
   *  @{ */
  void from_graph6(const std::string &str);
  void from_sparse6(const std::string &str);
  void from_edge_list(const std::string &str);

  void from_graph6(const std::ifstream &file);  
  void from_sparse6(const std::ifstream &file);
  void from_edge_list(const std::ifstream &file);
  /** @} */

  /**
   * @brief Creates a new Graph object containing only the specified vertices.
   * @param vertices The set of vertices to include in the subgraph.
   * @return An induced subgraph.
   */
  Graph subgraph(const std::unordered_set<Vertex> &vertices) const;

  // Generation ----------------------------------------------------------------
  /** @name Generation
  *  Methods for generating common types of graphs
  *  @{ */
  static std::vector<Graph> generate_random(std::size_t ord, std::size_t count, double edge_prob);
  static Graph generate_path(std::size_t ord);
  static Graph generate_cycle(std::size_t ord);
  static Graph generate_complete(std::size_t ord);
  static Graph generate_cubic(std::size_t ord);
  /** @} */

  // Element Access ----------------------------------------------------------------
  /** @brief Returns the number of vertices (Order). */
  std::size_t get_order() const;
  std::size_t get_vertex_count() const;

  /** @brief Returns the number of edges (Size). */
  std::size_t get_size() const;
  std::size_t get_edge_count() const;

  /** @brief Returns the number of edges connected to vertex u. */
  std::size_t get_degree(Vertex u) const;

  /** @brief Returns the highest degree present in the graph. */
  std::size_t get_max_degree() const;

  /** @brief Retrieves the external label for a given internal vertex index. */
  Vertex get_label(Vertex u) const;

  /** @brief Returns the set of vertices directly connected to vertex u. */
  VertexSet get_neighbors(Vertex u) const;
  VertexSet get_adjacent(Vertex u) const;

  // Query ----------------------------------------------------------------
  /** @brief Checks if an edge exists between u and v. */
  bool has_edge(Vertex u, Vertex v);

  /** @brief Returns a list of all vertex indices. */
  std::vector<Vertex> get_vertices() const;

  /** @brief Returns a list of all edges in the graph. */
  std::vector<Edge> get_edges() const;

  // Insertion ----------------------------------------------------------------
  /** @brief Adds an undirected edge between u and v. */
  void insert_edge(Vertex u, Vertex v);

  /** @brief Adds new vertices to the graph. */
  void insert_vertex(std::size_t count = 1);

  // Erasure ----------------------------------------------------------------
  /** @brief Removes the edge between u and v if it exists. */
  void erase_edge(Vertex u, Vertex v);

  /** @brief Removes a vertex and all its incident edges. */
  void erase_vertex(Vertex u);

  // Tree ----------------------------------------------------------------
  /** @brief Checks if the graph is acyclic. */
  bool is_valid_forest() const;

  /** @brief Checks if the graph is both connected and acyclic. */
  bool is_valid_tree() const;

  /** @brief Finds the longest path between any two nodes in a tree. */
  std::size_t tree_diameter() const;

  /** @brief Finds the vertex (or pair) that minimizes the distance to all other nodes. */
  std::pair<Vertex, Vertex> tree_center() const;

  // Misc ----------------------------------------------------------------
  /** @brief Returns all vertices with degree 1. */
  VertexSet get_pendants() const;

  /** @brief Given a set, returns all vertices in the graph NOT in that set. */
  VertexSet get_complement(const VertexSet &verts) const;

  /** @brief Checks if there is a path between every pair of vertices. */
  bool is_connected() const;

  // Output ----------------------------------------------------------------
  /** @brief Formats the graph for stream output (e.g., for debugging). */
  friend std::ostream& operator<<(std::ostream &os, const Graph &graph);
};

#endif