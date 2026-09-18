#include "../include/graph.hpp"

#include <algorithm>
#include <numeric>
#include <random>
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
    for (std::size_t a = 0; a < bit_count; a++) {
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
  std::vector<bool> visited(graph.order(), false);
  std::vector<std::size_t> dist(graph.order(), 0);

  std::queue<std::size_t> queue;
	queue.push(start);
	visited[start] = true;

	Vertex farthest = start;
	max_dist = 0;

	while (!queue.empty()) {
    Vertex u = queue.front();
		queue.pop();

		if (dist[u] > max_dist) {
			max_dist = dist[u];
			farthest = u;
		}

    VertexSet neighbors = graph.adjacent(u);
		for (VertexSet::const_iterator it_v = neighbors.cbegin(); it_v != neighbors.cend(); it_v++) {
			Vertex v = *it_v;
			if (visited[v]) continue;
			visited[v] = true;
			dist[v] = dist[u] + 1;
			queue.push(v);
		}
	}

	return farthest;
}

Graph::Graph(std::size_t order, bool label) :
  _vert_count(order),
  _edge_count(0)
{
  assert(order != 0);
  labels.resize(order);
  if (label) std::iota(labels.begin(), labels.end(), 0);
  _adjacent.resize(order);
}

void Graph::from_edge_list(const std::ifstream &file) {
  return from_edge_list(_file_to_string(file));
}

void Graph::from_edge_list(const std::string &str) {
  assert(!str.empty());

  std::stringstream ss(str);

  Vertex u;
  Vertex v;
  std::size_t ec_tmp;
  std::size_t vc_tmp;

  if (!(ss >> vc_tmp >> ec_tmp)) return;

  _vert_count = vc_tmp;
  _edge_count = 0;

  _adjacent.clear();
  _adjacent.resize(_vert_count);

  labels.clear();
  labels.resize(_vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  for (std::size_t a = 0; a < ec_tmp; a++) {
    if (ss >> u >> v) insert_edge(u, v);
    else break;
  }
}

void Graph::from_sparse6(const std::ifstream &file) {
  return from_sparse6(_file_to_string(file));
}

void Graph::from_sparse6(const std::string &str) {
  assert(!str.empty());

  std::string str_tmp = str.substr(1);
  for (char &c : str_tmp) c -= 63;

  size_t header_offset = _get_header_offset(str_tmp.data(), _vert_count);
  const char* data_ptr = str_tmp.data() + header_offset;
  size_t data_len = str_tmp.size() - header_offset;
    
  std::vector<bool> is_new_vertex_bit;
  std::vector<Vertex> neighbors;
  _parse_sparse6_data(_vert_count, data_len, data_ptr, is_new_vertex_bit, neighbors);

  _adjacent.clear();
  _adjacent.resize(_vert_count);

  labels.clear();
  labels.resize(_vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  Vertex u = 0;
  for (size_t a = 0; a < neighbors.size(); a++) {
    if (is_new_vertex_bit[a]) u++;
        
    Vertex v = neighbors[a];
    if (v >= _vert_count || u >= _vert_count) break;

    if (v > u) u = v;
    else insert_edge(u, v);
  }
}

void Graph::from_graph6(const std::ifstream &file) {
  return from_graph6(_file_to_string(file));
}

void Graph::from_graph6(const std::string &str) {
  assert(!str.empty());
    
  std::string str_tmp(str);
  for (char &c : str_tmp) c -= 63;

  size_t header_offset = _get_header_offset(str_tmp.data(), _vert_count);
  size_t triangle_bits = (_vert_count * (_vert_count - 1)) / 2;
    
  const char* data_start = str_tmp.data() + header_offset;
  size_t data_char_count = str_tmp.size() - header_offset;
  std::vector<bool> is_edge = _parse_graph6_data(data_start, data_char_count);

  _adjacent.clear();
  _adjacent.resize(_vert_count);

  labels.clear();
  labels.resize(_vert_count);
  std::iota(labels.begin(), labels.end(), 0);

  size_t bit_pos = 0;
  for (Vertex u = 0; u < _vert_count; u++) {
    for (Vertex v = 0; v < u; v++) {
      if (bit_pos < triangle_bits && bit_pos < is_edge.size() && is_edge[bit_pos]) {
        insert_edge(v, u);
      }
      bit_pos++;
    }
  }
}

Graph Graph::subgraph(const VertexBitset &vertices) {
  std::unordered_map<Vertex, Vertex> old_to_new;
  Graph induced(vertices.count(), false);

  std::size_t idx = 0;
  for (Vertex v = 0; v < _vert_count; v++) {
    if (!vertices.test(v)) continue;
    Vertex old_v = v;
    old_to_new[old_v] = idx;
    induced.labels[idx] = old_v;
    idx++;
  }
  

  for (Vertex u = 0; u < _vert_count; u++) {
    Vertex old_u = u;
    Vertex new_u = old_to_new[old_u];

    for (VertexSet::const_iterator it_v = _adjacent[old_u].cbegin(); it_v != _adjacent[old_u].cend(); it_v++) {
      Vertex old_v = *it_v;
      if (!vertices.test(old_v)) continue;
      Vertex new_v = old_to_new[old_v];
      if (new_u < new_v) induced.insert_edge(new_u, new_v);
    }
  }

  return induced;
}

Graph Graph::subgraph(const VertexSet &vertices) {
  std::unordered_map<Vertex, Vertex> old_to_new;
  Graph induced(vertices.size(), false);

  std::size_t idx = 0;
  for (VertexSet::const_iterator it_v = vertices.cbegin(); it_v != vertices.cend(); it_v++) {
    Vertex old_v = *it_v;
    old_to_new[old_v] = idx;
    induced.labels[idx] = old_v;
    idx++;
  }

  for (VertexSet::const_iterator it_u = vertices.cbegin(); it_u != vertices.cend(); it_u++) {
    Vertex old_u = *it_u;
    Vertex new_u = old_to_new[old_u];

    for (VertexSet::const_iterator it_v = _adjacent[old_u].cbegin(); it_v != _adjacent[old_u].cend(); it_v++) {
      Vertex old_v = *it_v;
      if (vertices.find(old_v) == vertices.cend()) continue;
      Vertex new_v = old_to_new[old_v];
      if (new_u < new_v) induced.insert_edge(new_u, new_v);
    }
  }

  return induced;
}

VertexBitset Graph::adjacent_bitset(Vertex u) const {
  VertexBitset adjacent_bs(_vert_count);
  VertexSet adjacent = _adjacent[u];
  for (Vertex u : adjacent) {
    adjacent_bs.set(u);
  }
  return adjacent_bs;
}

const VertexSet& Graph::adjacent(Vertex u) const {
  return _adjacent[u];
}

Vertex Graph::label(Vertex u) const {
  return labels[u];
}

std::size_t Graph::order() const {
  return _vert_count;
}

std::size_t Graph::size() const {
  return _edge_count;
}

std::size_t Graph::degree(Vertex u) const {
  return _adjacent[u].size();
}

std::size_t Graph::max_degree() const {
  std::size_t max = degree(0);
  for (Vertex u = 1; u < _vert_count; u++) {
    max = std::max(max, degree(u));
  }
  return max;
}

std::size_t Graph::min_degree() const {
  std::size_t min = degree(0);
  for (Vertex u = 1; u < _vert_count; u++) {
    min = std::min(min, degree(u));
  }
  return min;
}

bool Graph::has_edge(Vertex u, Vertex v) const {
  return _adjacent[u].find(v) != _adjacent[u].end();
}

void Graph::insert_edge(Vertex u, Vertex v) {
  if (has_edge(u, v)) return;
  _adjacent[u].insert(v);
  _adjacent[v].insert(u);
  _edge_count++;
}

void Graph::erase_edge(Vertex u, Vertex v) {
  if (!has_edge(u, v)) return;
  _adjacent[u].erase(v);
  _adjacent[v].erase(u);
  _edge_count--;
}

std::vector<Graph> GraphGenerator::random(std::size_t order, std::size_t count, double edge_probability) {
  std::vector<Graph> output;
  std::bernoulli_distribution dist(edge_probability);
  std::mt19937 gen(std::random_device{}());
  output.reserve(count);

  for (std::size_t a = 0; a < count; a++) {
    Graph graph(order);
    for (Vertex u = 0; u < order; u++) {
      for (Vertex v = u + 1; v < order; v++) {
        if (dist(gen)) graph.insert_edge(u, v);
      }
    }
    output.push_back(std::move(graph));
  }

  return output;
}

std::vector<Graph> GraphGenerator::cubic(std::size_t order, std::size_t count) {
  assert(order % 2 == 0 && order >= 4);
  std::mt19937 gen(std::random_device{}());

  std::vector<Graph> output;

  while (output.size() < count) {
    while (true) {
      Graph graph(order);
      std::vector<std::size_t> points;
      points.reserve(order * 3);
      for (std::size_t a = 0; a < order; a++) {
        points.push_back(a);
        points.push_back(a);
        points.push_back(a);
      }

      std::shuffle(points.begin(), points.end(), gen);

      bool valid = true;
      for (std::size_t a = 0; a < points.size(); a += 2) {
        Vertex u = points[a];
        Vertex v = points[a + 1];

        VertexSet neighbors(u);
        if (u == v || neighbors.find(v) != neighbors.cend()) { 
          valid = false;
          break;
        }

        graph.insert_edge(u, v);
      }

      if (valid) output.push_back(graph);
    }
  }

  return output;
}

Graph GraphGenerator::path(std::size_t order) {
  Graph graph(order);
  for (Vertex u = 0; u < order - 1; u++) {
    graph.insert_edge(u, u + 1);
  }
  return graph;
}

Graph GraphGenerator::cycle(std::size_t order) {
  Graph graph = GraphGenerator::path(order);
  if (order != 1) graph.insert_edge(order - 1, 0);
  return graph;
}

Graph GraphGenerator::complete(std::size_t order) {
  Graph graph(order);
  for (Vertex u = 0; u < order; u++) {
    for (Vertex v = u + 1; v < order; v++) {
      graph.insert_edge(u, v);
    }
  }
  return graph;
}

bool GraphTreeInfo::is_forest(const Graph &graph) {
  std::vector<bool> visited(graph._vert_count, false);

  for (std::size_t u = 0; u < graph._vert_count; u++) {
    if (visited[u]) continue;

    std::vector<std::pair<Vertex, Vertex>> stack;
    stack.emplace_back(u, u);
    visited[u] = true;

    while (!stack.empty()) {
      auto [curr, parent] = stack.back();
      stack.pop_back();

      for (Vertex neighbor : graph._adjacent[curr]) {
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

bool GraphTreeInfo::is_tree(const Graph &graph) {
  if (graph._edge_count != graph._vert_count - 1) return false;
  return GraphTreeInfo::is_forest(graph);
}

std::size_t GraphTreeInfo::diameter(const Graph &graph) {
  if (graph._vert_count == 0) return 0;

  std::size_t max_dist = 0;
	Vertex farthest = _bfs_farthest(graph, 0, max_dist);

	max_dist = 0;
	_bfs_farthest(graph, farthest, max_dist);

	return max_dist;
}

std::pair<Vertex, Vertex> GraphTreeInfo::center(const Graph &graph) {
  if (graph._vert_count == 1) return {0, INVALID_INDEX};
  if (graph._vert_count == 2) return {0, 1};

  std::vector<std::size_t> degrees(graph._vert_count);
  std::vector<Vertex> leaves;

  for (Vertex u = 0; u < graph._vert_count; u++) {
    degrees[u] = graph._adjacent[u].size();
    if (degrees[u] == 1) leaves.push_back(u);
  }

  std::size_t removed_count = leaves.size();
  while (removed_count < graph._vert_count) {
    std::vector<Vertex> next_leaves;
        
    for (Vertex leaf : leaves) {
      for (Vertex neighbor : graph._adjacent[leaf]) {
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

bool GraphMiscInfo::is_connected(const Graph &graph) {
	if (graph._vert_count == 0) return true;

  std::vector<bool> visited(graph._vert_count, false);
  std::vector<Vertex> stack;

  visited[0] = true;
	stack.push_back(0);
  std::size_t count = 1;

	while (!stack.empty()) {
		Vertex u = stack.back();
		stack.pop_back();

    for (Vertex v : graph._adjacent[u]) {
      if (visited[v]) continue;
      visited[v] = true;
      stack.push_back(v);
      count++;
    }
  }

  return count == graph._vert_count;
}

VertexBitset GraphMiscInfo::pendants_bs(const Graph &graph) {
  VertexBitset pendants(graph._vert_count);
  for (Vertex u = 0; u < graph._adjacent.size(); u++) {
    if (graph._adjacent[u].size() != 1) continue;
    pendants.set(u);
  }
  return pendants; 
}

VertexSet GraphMiscInfo::pendants(const Graph &graph) {
  VertexSet pendants;
  for (Vertex u = 0; u < graph._adjacent.size(); u++) {
    if (graph._adjacent[u].size() != 1) continue;
    pendants.insert(u);
  }
  return pendants; 
}

VertexBitset GraphMiscInfo::complement_bs(const Graph &graph, const VertexBitset &vertices) {
	VertexBitset complement(graph._vert_count);
	for (Vertex u = 0; u < graph._vert_count; u++){
		if (vertices.test(u)) continue;
		complement.set(u);
	}
	return complement;
}

VertexSet GraphMiscInfo::complement(const Graph &graph, const VertexSet &vertices) {
	VertexSet complement;
	for (Vertex u = 0; u < graph._vert_count; u++){
		if (vertices.find(u) != vertices.end()) continue;
		complement.insert(u);
	}
	return complement;
}
