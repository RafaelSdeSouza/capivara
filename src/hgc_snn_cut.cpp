#include <Rcpp.h>
#include <queue>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

struct QueueEdge {
  int u;
  int v;
  double w;
  double score;
};

struct EdgeCompare {
  bool operator()(const QueueEdge& a, const QueueEdge& b) const {
    return a.score < b.score;
  }
};

int hgc_find_root(std::vector<int>& parent, int x) {
  int r = x;
  while (parent[r] != r) r = parent[r];

  while (parent[x] != x) {
    int p = parent[x];
    parent[x] = r;
    x = p;
  }

  return r;
}

// [[Rcpp::export]]
IntegerVector capivara_hgc_snn_cut_cpp(IntegerMatrix nn_idx, int target_k) {
  int n = nn_idx.nrow();
  int k = nn_idx.ncol();

  if (target_k < 1 || target_k > n) {
    stop("target_k must be between 1 and n.");
  }

  std::vector< std::unordered_map<int, double> > adj(n);
  std::vector<int> mark(n, 0);

  // Build shared-nearest-neighbor graph.
  for (int i = 0; i < n; ++i) {
    int stamp = i + 1;

    for (int a = 0; a < k; ++a) {
      int nb = nn_idx(i, a) - 1;
      if (nb >= 0 && nb < n && nb != i) {
        mark[nb] = stamp;
      }
    }

    for (int a = 0; a < k; ++a) {
      int j = nn_idx(i, a) - 1;
      if (j < 0 || j >= n || j == i) continue;

      int shared = 0;

      for (int b = 0; b < k; ++b) {
        int nbj = nn_idx(j, b) - 1;
        if (nbj >= 0 && nbj < n && mark[nbj] == stamp) {
          ++shared;
        }
      }

      if (shared <= 0) continue;

      double w = static_cast<double>(shared);
      adj[i][j] += w;
      adj[j][i] += w;
    }
  }

  std::vector<double> degree(n, 0.0);

  for (int i = 0; i < n; ++i) {
    for (std::unordered_map<int, double>::const_iterator it = adj[i].begin();
         it != adj[i].end(); ++it) {
      degree[i] += it->second;
    }
  }

  std::priority_queue<QueueEdge, std::vector<QueueEdge>, EdgeCompare> pq;

  for (int i = 0; i < n; ++i) {
    for (std::unordered_map<int, double>::const_iterator it = adj[i].begin();
         it != adj[i].end(); ++it) {
      int j = it->first;

      if (i < j && degree[i] > 0 && degree[j] > 0) {
        double w = it->second;
        pq.push(QueueEdge{i, j, w, w / (degree[i] * degree[j])});
      }
    }
  }

  std::vector<int> parent(n);
  std::vector<char> active(n, 1);

  for (int i = 0; i < n; ++i) {
    parent[i] = i;
  }

  int active_count = n;

  while (active_count > target_k && !pq.empty()) {
    QueueEdge e = pq.top();
    pq.pop();

    int u = hgc_find_root(parent, e.u);
    int v = hgc_find_root(parent, e.v);

    if (u == v || !active[u] || !active[v]) continue;

    std::unordered_map<int, double>::iterator uv_it = adj[u].find(v);
    if (uv_it == adj[u].end()) continue;

    double current_w = uv_it->second;
    if (std::fabs(current_w - e.w) > 1e-8) continue;

    if (degree[u] <= 0 || degree[v] <= 0) continue;

    if (adj[u].size() < adj[v].size()) {
      std::swap(u, v);
    }

    double internal_w = 0.0;
    std::unordered_map<int, double>::iterator int_it = adj[u].find(v);
    if (int_it != adj[u].end()) {
      internal_w = int_it->second;
    }

    std::vector< std::pair<int, double> > neigh_v;
    neigh_v.reserve(adj[v].size());

    for (std::unordered_map<int, double>::const_iterator it = adj[v].begin();
         it != adj[v].end(); ++it) {
      neigh_v.push_back(*it);
    }

    adj[u].erase(v);
    adj[v].erase(u);

    degree[u] = degree[u] + degree[v] - 2.0 * internal_w;
    degree[v] = 0.0;

    active[v] = 0;
    parent[v] = u;
    --active_count;

    for (size_t a = 0; a < neigh_v.size(); ++a) {
      int t = hgc_find_root(parent, neigh_v[a].first);
      double wvt = neigh_v[a].second;

      if (t == u || !active[t]) continue;

      adj[t].erase(v);

      double new_w = wvt;
      std::unordered_map<int, double>::iterator old = adj[u].find(t);

      if (old != adj[u].end()) {
        new_w += old->second;
      }

      adj[u][t] = new_w;
      adj[t][u] = new_w;

      if (degree[u] > 0 && degree[t] > 0) {
        pq.push(QueueEdge{u, t, new_w, new_w / (degree[u] * degree[t])});
      }
    }

    adj[v].clear();
  }

  IntegerVector out(n);
  std::unordered_map<int, int> remap;
  int next_label = 1;

  for (int i = 0; i < n; ++i) {
    int r = hgc_find_root(parent, i);

    if (remap.find(r) == remap.end()) {
      remap[r] = next_label++;
    }

    out[i] = remap[r];
  }

  return out;
}
