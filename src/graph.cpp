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

void _parse_sparse6_data(std::size_t order, std::size_t count, const char *edges, std::vector<bool> &is_vertex, std::vector<VertexIndex> &neighbors) {
  if (order == 0) return;

  std::size_t bits_per = 0;
  while ((1ULL << bits_per) < order) bits_per++;

  std::size_t current_byte_idx = 0;
  std::size_t bits_left_in_byte = 0;
  std::size_t current_byte = 0;

  auto get_next_bit = [&]() mutable -> bool {
    if (bits_left_in_byte == 0) {
      current_byte = static_cast<unsigned char>(edges[current_byte_idx++]);
      bits_left_in_byte = 6;
    }
    bits_left_in_byte--;
    return (current_byte >> bits_left_in_byte) & 1;
  };

  auto get_next_k_bits = [&](std::size_t bit_count) mutable -> VertexIndex {
    VertexIndex val = 0;
    for (std::size_t i = 0; i < bit_count; ++i) {
      val = (val << 1) | (get_next_bit() ? 1 : 0);
    }
    return val;
  };

  std::size_t total_bits = count * 6;
  std::size_t bits_consumed = 0;

  while (bits_consumed + 1 + bits_per <= total_bits) {
    bool bit = get_next_bit();
    VertexIndex val = get_next_k_bits(bits_per);
    bits_consumed += (1 + bits_per);

    is_vertex.push_back(bit);
    neighbors.push_back(val);
  }
}

VertexIndex _bfs_farthest(const Graph &graph, VertexIndex start, std::size_t &max_dist) {
  std::vector<bool> visited(graph.get_order(), false);
  std::vector<std::size_t> dist(graph.get_order(), 0);

  std::queue<std::size_t> queue;
	queue.push(start);
	visited[start] = true;

	VertexIndex farthest = start;
	max_dist = 0;

	while (!queue.empty()) {
    VertexIndex vert1 = queue.front();
		queue.pop();

		if (dist[vert1] > max_dist) {
			max_dist = dist[vert1];
			farthest = vert1;
		}

    AdjacencySet neighbors = graph.get_adjacent(vert1);
		for (AdjacencySet::const_iterator it = neighbors.cbegin(); it != neighbors.cend(); it++) {
			VertexIndex vert2 = *it;
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
  for (std::size_t a = 0; a < vert_count; a++) {
    for (std::size_t b = 0; b < a; b++) {
      if (bit_pos < triangle_bits && bit_pos < is_edge.size() && is_edge[bit_pos]) {
        insert_edge(b, a);
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
    
  std::vector<VertexIndex> neighbors;
  std::vector<bool> is_new_vertex_bit;
  const char* data_ptr = str_tmp.data() + header_offset;
  size_t data_len = str_tmp.size() - header_offset;
    
  _parse_sparse6_data(vert_count, data_len, data_ptr, is_new_vertex_bit, neighbors);

  adj.resize(vert_count);
  labels.resize(vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  VertexIndex u = 0;
  for (size_t i = 0; i < neighbors.size(); ++i) {
    if (is_new_vertex_bit[i]) u++;
        
    VertexIndex v = neighbors[i];
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

  VertexIndex u;
  VertexIndex v;
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
  std::unordered_map<VertexIndex, VertexIndex> old_to_new;

  std::size_t idx = 0;
  for (VertexSet::const_iterator it = vertices.cbegin(); it != vertices.cend(); it++) {
    VertexIndex old_v = *it;
    old_to_new[old_v] = idx;
    induced.labels[idx] = old_v;
    idx++;
  }

  for (VertexSet::const_iterator it1 = vertices.cbegin(); it1 != vertices.cend(); it1++) {
    VertexIndex old_u = *it1;
    VertexIndex new_u = old_to_new[old_u];

    for (VertexSet::const_iterator it2 = adj[old_u].cbegin(); it2 != adj[old_u].cend(); it2++) {
      VertexIndex old_v = *it2;
      if (vertices.find(old_v) == vertices.cend()) continue;
      VertexIndex new_v = old_to_new[old_v];
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

std::size_t Graph::get_degree(VertexIndex vert) const {
  assert(vert < vert_count);
  return adj[vert].size();
}

std::size_t Graph::get_max_degree() const {
  if (vert_count == 0) return INVALID_INDEX;

  std::size_t max = get_degree(0);
  for (std::size_t a = 1; a < get_order(); a++) {
    max = std::max(max, get_degree(a));
  }

  return max;
}

VertexIndex Graph::get_label(VertexIndex vert) const {
  assert(vert < vert_count);
  return labels[vert];
}

AdjacencySet Graph::get_adjacent(VertexIndex vert) const {
  assert(vert < vert_count);
  return adj[vert];
}

// Query ----------------------------------------------------------------
bool Graph::has_edge(VertexIndex vert1, VertexIndex vert2) {
  if (vert1 >= vert_count) return false;
  return adj[vert1].find(vert2) != adj[vert1].cend();
}

std::vector<VertexIndex> Graph::get_vertices() const {
  return labels;
}

std::vector<EdgePair> Graph::get_edges() const {
  std::vector<EdgePair> edges;
  edges.reserve(edge_count);

  for (VertexIndex a = 0; a < adj.size(); a++) {
    for (VertexIndex b : adj[a]) {
      if (a >= b) continue;
      edges.emplace_back(a, b);
    }
  }

  return edges;
}

// Insertion ----------------------------------------------------------------
void Graph::insert_edge(VertexIndex vert1, VertexIndex vert2) {
  assert(vert1 < vert_count && vert2 < vert_count);
  if (has_edge(vert1, vert2)) return;
	adj[vert1].insert(vert2);
	adj[vert2].insert(vert1);
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
  for (size_t a = 0; a < count; a++) labels.push_back(++max);
  vert_count += count;
}

// Erasure ----------------------------------------------------------------
void Graph::erase_edge(VertexIndex vert1, VertexIndex vert2) {
  assert(vert1 < vert_count && vert2 < vert_count);
  if (!has_edge(vert1, vert2)) return;
	adj[vert1].erase(vert2);
	adj[vert2].erase(vert1);
	edge_count--;
}

void Graph::erase_vertex(VertexIndex vert) {
  assert(vert < vert_count);

  for (VertexIndex neighbor : adj[vert]) {
    adj[neighbor].erase(vert);
  }

  adj.erase(adj.begin() + vert);
  labels.erase(labels.begin() + vert);
  vert_count--;

  for (AdjacencySet &neighbors : adj) {
    AdjacencySet updated_set;
    for (VertexIndex neighbor_idx : neighbors) {
      if (neighbor_idx > vert) updated_set.insert(neighbor_idx - 1);
      else updated_set.insert(neighbor_idx);
    }
    neighbors = std::move(updated_set);
  }
}

// Tree ----------------------------------------------------------------
bool Graph::is_valid_forest() const {
  if (vert_count == 0) return true;

  std::vector<bool> visited(vert_count, false);

  for (std::size_t a = 0; a < vert_count; a++) {
    if (visited[a]) continue;

    std::vector<std::pair<VertexIndex, VertexIndex>> stack;
    stack.push_back({a, a});
    visited[a] = true;

    while (!stack.empty()) {
      auto [curr, parent] = stack.back();
      stack.pop_back();

      for (VertexIndex neighbor : adj[curr]) {
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
	VertexIndex farthest = _bfs_farthest(*this, 0, max_dist);

	max_dist = 0;
	_bfs_farthest(*this, farthest, max_dist);

	return max_dist;
}

std::pair<VertexIndex, VertexIndex> Graph::tree_center() const {
  if (vert_count == 0) return {INVALID_INDEX, INVALID_INDEX};
  if (vert_count == 1) return {0, INVALID_INDEX};
  if (vert_count == 2) return {0, 1};

  std::vector<std::size_t> degrees(vert_count);
  std::vector<VertexIndex> leaves;

  for (VertexIndex a = 0; a < vert_count; a++) {
    degrees[a] = adj[a].size();
    if (degrees[a] == 1) leaves.push_back(a);
  }

  std::size_t removed_count = leaves.size();
  while (removed_count < vert_count) {
    std::vector<VertexIndex> next_leaves;
        
    for (VertexIndex leaf : leaves) {
      for (VertexIndex neighbor : adj[leaf]) {
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

// Zero Forcing ----------------------------------------------------------------
bool Graph::is_valid_zf(const VertexSet &filled) const {
  std::size_t pt;
  zf_closure(filled, &pt);
  return pt != INVALID_INDEX;
}

VertexSet Graph::zf_closure(VertexSet filled, std::size_t *pt) const {
  std::size_t propagation_time;
  VertexIndex vert;

  for (propagation_time = 0; propagation_time < vert_count; propagation_time++) {
    VertexSet active;

    for (VertexSet::const_iterator it_u = filled.cbegin(); it_u != filled.cend(); it_u++) {
      std::size_t count = 0;
      AdjacencySet neighbors = get_adjacent(*it_u);
      
      for (AdjacencySet::const_iterator it_v = neighbors.cbegin(); it_v != neighbors.cend(); it_v++) {
        if (filled.find(*it_v) != filled.cend()) continue;
        count++;
        vert = *it_v;
      }
      if (count == 1) active.insert(vert);
    }

    if (active.empty()) break;
    filled.insert(active.cbegin(), active.cend()); 
  }

  if (pt) {
    if (filled.size() == vert_count) *pt = propagation_time;
    else *pt = INVALID_INDEX;
  }

  return filled;
}

std::size_t Graph::zf_wavefront() const {
  if (vert_count == 0) return 0;

  typedef std::pair<VertexSet, std::size_t> ZFSolution;
  std::vector<ZFSolution> cl_pairs;
  cl_pairs.emplace_back(VertexSet(), 0);

  for (std::size_t rank_limit = 1; rank_limit <= vert_count; rank_limit++) {
    for (std::vector<ZFSolution>::const_iterator it = cl_pairs.cbegin(); it != cl_pairs.cend(); it++) {
      auto [source, source_rank] = *it;

      for (VertexIndex a = 0; a < vert_count; a++) {
        VertexSet candidate(source.cbegin(), source.cend());

        candidate.insert(a);
        AdjacencySet neighbors = get_adjacent(a);
        candidate.insert(neighbors.cbegin(), neighbors.cend());

        zf_closure(candidate);

        std::size_t candidate_rank = source_rank;
        if (source.find(a) == source.cend()) candidate_rank++;

        int neighbors_outside_source = 0;
        for (AdjacencySet::const_iterator it = neighbors.cbegin(); it != neighbors.cend(); it++) {
          if (source.find(*it) != source.cend()) continue;
          neighbors_outside_source++;
        }

        if (neighbors_outside_source > 1) candidate_rank += (neighbors_outside_source - 1);

        if (candidate_rank <= rank_limit) {
          bool already_present = false;

          for (std::vector<ZFSolution>::const_iterator it = cl_pairs.cbegin(); it != cl_pairs.cend(); it++) {
            if (it->first != candidate || it->second > candidate_rank) continue;
            already_present = true;
            break;
          }

          if (!already_present) {
            cl_pairs.emplace_back(candidate, candidate_rank);
            if (candidate.size() == vert_count) return candidate_rank;
          }
        }
      }
    }
  }

  return vert_count;
}

// Misc ----------------------------------------------------------------
VertexSet Graph::get_pendants() const {
  VertexSet pendants;

  for (VertexIndex a = 0; a < adj.size(); a++) {
    if (adj[a].size() != 1) continue;
    pendants.insert(a);
  }

  return pendants; 
}

bool Graph::is_connected() const {
	if (vert_count == 0) return true;

  std::vector<bool> visited(vert_count, false);
  std::vector<VertexIndex> stack;

  visited[0] = true;
	stack.push_back(0);
  std::size_t count = 1;

	while (!stack.empty()) {
		VertexIndex vert1 = stack.back();
		stack.pop_back();

    for (VertexIndex vert2 : adj[vert1]) {
      if (visited[vert2]) continue;
      visited[vert2] = true;
      stack.push_back(vert2);
      count++;
    }
  }

  return count == vert_count;
}

// Output ----------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const Graph &graph) {
  os << "Order: " << graph.vert_count << ", #Edges: " << graph.edge_count << "\n";
  for (std::size_t a = 0; a < graph.vert_count; a++){
    os << graph.labels[a] << ": ";
    for (AdjacencySet::const_iterator it = graph.adj[a].cbegin(); it != graph.adj[a].cend(); it++) {
      os << graph.labels[*it] << " ";
    }
    os << "\n";
  }
  return os;
}