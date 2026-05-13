#include "graph.hpp"

#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cassert>
#include <sstream>
#include <queue>

std::string _file_to_string(const std::ifstream &file) {
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

size_t _get_header_offset(const char *data, std::size_t &order) {
  if (data[0] <= 62) {
    order = data[0];
    return 1;
  }

  if (data[1] <= 62) {
    order = (data[1] << 12) + (data[2] << 6) + data[3];
    return 4;
  }

  order = (static_cast<size_t>(data[2]) << 30) + (data[3] << 24) + (data[4] << 18) + (data[5] << 12) + (data[6] << 6) + data[7];
  return 8;
}

std::vector<bool> _parse_graph6_data(const char *edges, std::size_t count) {
  std::vector<bool> output(count*6);
	for (std::size_t a = 0; a < count; a++) {
		for (std::ptrdiff_t b = 5; b >= 0; b--) {
			output[a*6 + (5 - b)] = (edges[a] >> b) & 1;
		}
	}
	return output;
}

void _parse_sparse6_data(std::size_t order, std::size_t count, const char *edges, std::vector<bool> &is_vertex, std::vector<Vertex> &neighbors) {
  if (order == 0) return;

  std::size_t bits_per = 0;
  while ((1ULL << bits_per) < order) bits_per++;

  std::size_t current_byte_idx = 0;
  std::size_t bits_left_in_byte = 0;
  std::size_t current_byte = 0;

  auto get_next_bit = [&]() -> bool {
    if (bits_left_in_byte == 0) {
      current_byte = static_cast<unsigned char>(edges[current_byte_idx++]);
      bits_left_in_byte = 6;
    }
    bits_left_in_byte--;
    return (current_byte >> bits_left_in_byte) & 1;
  };

  auto get_next_k_bits = [&](std::size_t bit_count) -> Vertex {
    Vertex val = 0;
    for (std::size_t i = 0; i < bit_count; ++i) {
      val = (val << 1) | (get_next_bit() ? 1 : 0);
    }
    return val;
  };

  std::size_t total_bits = count * 6;
  std::size_t bits_consumed = 0;

  while (bits_consumed + 1 + bits_per <= total_bits) {
    bool bit = get_next_bit();
    Vertex val = get_next_k_bits(bits_per);
    bits_consumed += (1 + bits_per);

    is_vertex.push_back(bit);
    neighbors.push_back(val);
  }
}

Vertex _bfs_farthest(const Graph &graph, Vertex start, std::size_t &max_dist) {
  std::vector<bool> visited(graph.get_order(), false);
  std::vector<std::size_t> dist(graph.get_order(), 0);

  std::queue<std::size_t> queue;
	queue.push(start);
	visited[start] = true;

	Vertex farthest = start;
	max_dist = 0;

	while (!queue.empty()) {
    Vertex vert1 = queue.front();
		queue.pop();

		if (dist[vert1] > max_dist) {
			max_dist = dist[vert1];
			farthest = vert1;
		}

    VertexSet neighbors = graph.get_adjacent(vert1);
		for (VertexSet::const_iterator it = neighbors.cbegin(); it != neighbors.cend(); it++) {
			Vertex vert2 = *it;
			if (visited[vert2]) continue;
			visited[vert2] = true;
			dist[vert2] = dist[vert1] + 1;
			queue.push(vert2);
		}
	}

	return farthest;
}

// Initialization ----------------------------------------------------------------
Graph::Graph() : 
  vert_count(0), edge_count(0)
{}

Graph::Graph(std::size_t ord, bool label) :
  vert_count(ord), edge_count(0)
{ 
  assert(ord > 0);
  labels.resize(ord);
  if (label) std::iota(labels.begin(), labels.end(), 0);
  adj.resize(ord);
}

void Graph::clear() {
  vert_count = 0;
  edge_count = 0;
  labels.resize(0);
  adj.resize(0);
}

void Graph::clear(std::size_t ord, bool label) {
  assert(ord > 0);
  vert_count = ord;
  edge_count = 0;
  labels.resize(ord);
  if (label) std::iota(labels.begin(), labels.end(), 0);
  adj.resize(0);
}

void Graph::from_graph6(const std::string &str) {
  if (str.empty()) {
    clear();
    return;
  }
    
  std::string str_tmp(str);
  for (char &c : str_tmp) c -= 63;

  clear();
  size_t header_offset = _get_header_offset(str_tmp.data(), vert_count);
    
  size_t triangle_bits = (vert_count * (vert_count - 1)) / 2;
    
  const char* data_start = str_tmp.data() + header_offset;
  size_t data_char_count = str_tmp.size() - header_offset;
  std::vector<bool> is_edge = _parse_graph6_data(data_start, data_char_count);

  adj.resize(vert_count);
  labels.resize(vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  size_t bit_pos = 0;
  for (Vertex u = 0; u < vert_count; u++) {
    for (Vertex v = 0; v < u; v++) {
      if (bit_pos < triangle_bits && bit_pos < is_edge.size() && is_edge[bit_pos]) {
        insert_edge(v, u);
      }
      bit_pos++;
    }
  }
}

void Graph::from_sparse6(const std::string &str) {
  if (str.empty() || str[0] != ':') {
    clear();
    return;
  }

  std::string str_tmp = str.substr(1);
  for (char &c : str_tmp) c -= 63;

  clear();
  size_t header_offset = _get_header_offset(str_tmp.data(), vert_count);
    
  std::vector<Vertex> neighbors;
  std::vector<bool> is_new_vertex_bit;
  const char* data_ptr = str_tmp.data() + header_offset;
  size_t data_len = str_tmp.size() - header_offset;
    
  _parse_sparse6_data(vert_count, data_len, data_ptr, is_new_vertex_bit, neighbors);

  adj.resize(vert_count);
  labels.resize(vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  Vertex u = 0;
  for (size_t a = 0; a < neighbors.size(); a++) {
    if (is_new_vertex_bit[a]) u++;
        
    Vertex v = neighbors[a];
    if (v >= vert_count || u >= vert_count) break;

    if (v > u) u = v;
    else insert_edge(u, v);
  }
}

void Graph::from_edge_list(const std::string &str) {
  if (str.empty()){
    clear();
    return;
  }

  std::stringstream ss(str);

  Vertex u;
  Vertex v;
  std::size_t ec_tmp;
  std::size_t vc_tmp;

  clear();

  if (!(ss >> vc_tmp >> ec_tmp)) return;

  vert_count = vc_tmp;
  edge_count = 0;
  adj.resize(vert_count);
  labels.resize(vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  for (std::size_t a = 0; a < ec_tmp; a++) {
    if (ss >> u >> v) insert_edge(u, v);
    else break;
  }
}

void Graph::from_graph6(const std::ifstream &file) {
  return from_graph6(_file_to_string(file));
}
void Graph::from_sparse6(const std::ifstream &file) {
  return from_sparse6(_file_to_string(file));
}
void Graph::from_edge_list(const std::ifstream &file) {
  return from_edge_list(_file_to_string(file));
}

Graph Graph::subgraph(const VertexSet &vertices) const {
  Graph induced(vertices.size(), false);
  std::unordered_map<Vertex, Vertex> old_to_new;

  std::size_t idx = 0;
  for (VertexSet::const_iterator it = vertices.cbegin(); it != vertices.cend(); it++) {
    Vertex old_v = *it;
    old_to_new[old_v] = idx;
    induced.labels[idx] = old_v;
    idx++;
  }

  for (VertexSet::const_iterator it_u = vertices.cbegin(); it_u != vertices.cend(); it_u++) {
    Vertex old_u = *it_u;
    Vertex new_u = old_to_new[old_u];

    for (VertexSet::const_iterator it_v = adj[old_u].cbegin(); it_v != adj[old_u].cend(); it_v++) {
      Vertex old_v = *it_v;
      if (vertices.find(old_v) == vertices.cend()) continue;
      Vertex new_v = old_to_new[old_v];
      if (new_u < new_v) induced.insert_edge(new_u, new_v);
    }
  }

  return induced;
}

// Element Access ----------------------------------------------------------------
std::size_t Graph::get_order() const {
  return get_vertex_count();
}

std::size_t Graph::get_vertex_count() const {
  return vert_count;
}

std::size_t Graph::get_edge_count() const {
  return edge_count;
}

std::size_t Graph::get_degree(Vertex u) const {
  assert(u < vert_count);
  return adj[u].size();
}

std::size_t Graph::get_max_degree() const {
  if (vert_count == 0) return INVALID_INDEX;

  std::size_t max = get_degree(0);
  for (Vertex u = 1; u < vert_count; u++) {
    max = std::max(max, get_degree(u));
  }

  return max;
}

Vertex Graph::get_label(Vertex u) const {
  assert(u < vert_count);
  return labels[u];
}

VertexSet Graph::get_adjacent(Vertex u) const {
  assert(u < vert_count);
  return adj[u];
}

// Query ----------------------------------------------------------------
bool Graph::has_edge(Vertex u, Vertex v) {
  if (u >= vert_count) return false;
  return adj[u].find(v) != adj[u].cend();
}

std::vector<Vertex> Graph::get_vertices() const {
  return labels;
}

std::vector<Edge> Graph::get_edges() const {
  std::vector<Edge> edges;
  edges.reserve(edge_count);

  for (Vertex u = 0; u < adj.size(); u++) {
    for (Vertex v : adj[u]) {
      if (u >= v) continue;
      edges.emplace_back(u, v);
    }
  }

  return edges;
}

// Insertion ----------------------------------------------------------------
void Graph::insert_edge(Vertex u, Vertex v) {
  assert(u < vert_count && v < vert_count);
  if (has_edge(u, v)) return;
	adj[u].insert(v);
	adj[v].insert(u);
  edge_count++;
}

void Graph::insert_vertex(std::size_t count) {
  if (count == 0) return;
  if (vert_count == 0) {
    clear(count);
    return;
  }

  std::size_t max = *std::max_element(labels.begin(), labels.end());
  labels.reserve(labels.size() + count);
  adj.resize(adj.size() + count);
  for (Vertex u = 0; u < count; u++) labels.push_back(++max);
  vert_count += count;
}

// Erasure ----------------------------------------------------------------
void Graph::erase_edge(Vertex u, Vertex v) {
  assert(u < vert_count && v < vert_count);
  if (!has_edge(u, v)) return;
	adj[u].erase(v);
	adj[v].erase(u);
	edge_count--;
}

void Graph::erase_vertex(Vertex u) {
  assert(u < vert_count);

  for (Vertex neighbor : adj[u]) {
    adj[neighbor].erase(u);
  }

  adj.erase(adj.begin() + u);
  labels.erase(labels.begin() + u);
  vert_count--;

  for (VertexSet &neighbors : adj) {
    VertexSet updated_set;
    for (Vertex neighbor_idx : neighbors) {
      if (neighbor_idx > u) updated_set.insert(neighbor_idx - 1);
      else updated_set.insert(neighbor_idx);
    }
    neighbors = std::move(updated_set);
  }
}

// Tree ----------------------------------------------------------------
bool Graph::is_valid_forest() const {
  if (vert_count == 0) return true;

  std::vector<bool> visited(vert_count, false);

  for (std::size_t u = 0; u < vert_count; u++) {
    if (visited[u]) continue;

    std::vector<std::pair<Vertex, Vertex>> stack;
    stack.emplace_back(u, u);
    visited[u] = true;

    while (!stack.empty()) {
      auto [curr, parent] = stack.back();
      stack.pop_back();

      for (Vertex neighbor : adj[curr]) {
        if (visited[neighbor]) {
          if (neighbor != parent) return false;
          continue;
        }
        
        stack.emplace_back(neighbor, curr);
        visited[neighbor] = true;
      }
    }
  }

  return true;
}

bool Graph::is_valid_tree() const {
  if (vert_count == 0) return false;
  if (edge_count != vert_count - 1) return false;
  return is_valid_forest();
}

std::size_t Graph::tree_diameter() const {
  if (vert_count == 0) return 0;

  std::size_t max_dist = 0;
	Vertex farthest = _bfs_farthest(*this, 0, max_dist);

	max_dist = 0;
	_bfs_farthest(*this, farthest, max_dist);

	return max_dist;
}

std::pair<Vertex, Vertex> Graph::tree_center() const {
  if (vert_count == 0) return {INVALID_INDEX, INVALID_INDEX};
  if (vert_count == 1) return {0, INVALID_INDEX};
  if (vert_count == 2) return {0, 1};

  std::vector<std::size_t> degrees(vert_count);
  std::vector<Vertex> leaves;

  for (Vertex u = 0; u < vert_count; u++) {
    degrees[u] = adj[u].size();
    if (degrees[u] == 1) leaves.push_back(u);
  }

  std::size_t removed_count = leaves.size();
  while (removed_count < vert_count) {
    std::vector<Vertex> next_leaves;
        
    for (Vertex leaf : leaves) {
      for (Vertex neighbor : adj[leaf]) {
        if (degrees[neighbor] <= 1) continue;
        if (--degrees[neighbor] == 1) next_leaves.push_back(neighbor);
      }
      degrees[leaf] = 0;
    }
        
    if (next_leaves.empty()) break;

    removed_count += next_leaves.size();
    leaves = std::move(next_leaves);
  }

  if (leaves.size() == 1) return {leaves[0], INVALID_INDEX};
  if (leaves.size() == 2) return {leaves[0], leaves[1]};
  return {INVALID_INDEX, INVALID_INDEX};
}

// Misc ----------------------------------------------------------------
VertexSet Graph::get_pendants() const {
  VertexSet pendants;

  for (Vertex u = 0; u < adj.size(); u++) {
    if (adj[u].size() != 1) continue;
    pendants.insert(u);
  }

  return pendants; 
}

VertexSet Graph::get_complement(const VertexSet &verts) const {
	VertexSet comp;
	for (Vertex u = 0; u < vert_count; u++){
		if(verts.find(u) != verts.cend()) continue;
		comp.insert(u);
	}
	return comp;
}

bool Graph::is_connected() const {
	if (vert_count == 0) return true;

  std::vector<bool> visited(vert_count, false);
  std::vector<Vertex> stack;

  visited[0] = true;
	stack.push_back(0);
  std::size_t count = 1;

	while (!stack.empty()) {
		Vertex u = stack.back();
		stack.pop_back();

    for (Vertex v : adj[u]) {
      if (visited[v]) continue;
      visited[v] = true;
      stack.push_back(v);
      count++;
    }
  }

  return count == vert_count;
}

// Output ----------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const Graph &graph) {
  os << "Order: " << graph.vert_count << ", #Edges: " << graph.edge_count << "\n";
  for (std::size_t u = 0; u < graph.vert_count; u++){
    os << graph.labels[u] << ": ";
    for (VertexSet::const_iterator it = graph.adj[u].cbegin(); it != graph.adj[u].cend(); it++) {
      os << graph.labels[*it] << " ";
    }
    os << "\n";
  }
  return os;
}