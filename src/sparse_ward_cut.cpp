#include <Rcpp.h>
#include <queue>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

using namespace Rcpp;

struct SparseWardEdge {
  int u;
  int v;
  double score;
};

struct SparseWardEdgeCompare {
  bool operator()(const SparseWardEdge& a, const SparseWardEdge& b) const {
    return a.score > b.score;
  }
};

int sparse_ward_find_root(std::vector<int>& parent, int x) {
  int r = x;
  while (parent[r] != r) r = parent[r];

  while (parent[x] != x) {
    int p = parent[x];
    parent[x] = r;
    x = p;
  }

  return r;
}

double sparse_ward_score(const NumericMatrix& center,
                         const std::vector<double>& size,
                         int u,
                         int v) {
  double ss = 0.0;
  int p = center.ncol();

  for (int j = 0; j < p; ++j) {
    double d = center(u, j) - center(v, j);
    ss += d * d;
  }

  return (size[u] * size[v] / (size[u] + size[v])) * ss;
}

// [[Rcpp::export]]
IntegerVector capivara_sparse_ward_cut_cpp(NumericMatrix features,
                                           IntegerMatrix nn_idx,
                                           int target_k) {
  int n = features.nrow();
  int p = features.ncol();
  int k = nn_idx.ncol();

  if (target_k < 1 || target_k > n) {
    stop("target_k must be between 1 and n.");
  }

  NumericMatrix center(clone(features));
  std::vector<double> size(n, 1.0);
  std::vector<int> parent(n);
  std::vector<char> active(n, 1);
  std::vector< std::unordered_set<int> > nbr(n);

  for (int i = 0; i < n; ++i) {
    parent[i] = i;
  }

  for (int i = 0; i < n; ++i) {
    for (int a = 0; a < k; ++a) {
      int j = nn_idx(i, a) - 1;
      if (j >= 0 && j < n && j != i) {
        nbr[i].insert(j);
        nbr[j].insert(i);
      }
    }
  }

  std::priority_queue<
    SparseWardEdge,
    std::vector<SparseWardEdge>,
    SparseWardEdgeCompare
  > pq;

  for (int i = 0; i < n; ++i) {
    for (std::unordered_set<int>::const_iterator it = nbr[i].begin();
         it != nbr[i].end(); ++it) {
      int j = *it;
      if (i < j) {
        pq.push(SparseWardEdge{i, j, sparse_ward_score(center, size, i, j)});
      }
    }
  }

  int active_count = n;

  while (active_count > target_k && !pq.empty()) {
    SparseWardEdge e = pq.top();
    pq.pop();

    int u = sparse_ward_find_root(parent, e.u);
    int v = sparse_ward_find_root(parent, e.v);

    if (u == v || !active[u] || !active[v]) continue;
    if (nbr[u].find(v) == nbr[u].end()) continue;

    double current = sparse_ward_score(center, size, u, v);
    if (std::fabs(current - e.score) > 1e-8 * (1.0 + std::fabs(current))) {
      continue;
    }

    if (nbr[u].size() < nbr[v].size()) {
      std::swap(u, v);
    }

    double new_size = size[u] + size[v];
    for (int j = 0; j < p; ++j) {
      center(u, j) = (size[u] * center(u, j) + size[v] * center(v, j)) / new_size;
    }
    size[u] = new_size;
    size[v] = 0.0;
    active[v] = 0;
    parent[v] = u;
    --active_count;

    nbr[u].erase(v);
    nbr[v].erase(u);

    std::vector<int> v_neigh;
    v_neigh.reserve(nbr[v].size());
    for (std::unordered_set<int>::const_iterator it = nbr[v].begin();
         it != nbr[v].end(); ++it) {
      v_neigh.push_back(*it);
    }

    for (size_t a = 0; a < v_neigh.size(); ++a) {
      int t = sparse_ward_find_root(parent, v_neigh[a]);
      if (t == u || !active[t]) continue;

      nbr[t].erase(v);
      nbr[t].insert(u);
      nbr[u].insert(t);
    }
    nbr[v].clear();

    for (std::unordered_set<int>::const_iterator it = nbr[u].begin();
         it != nbr[u].end(); ++it) {
      int t = sparse_ward_find_root(parent, *it);
      if (t == u || !active[t]) continue;
      pq.push(SparseWardEdge{u, t, sparse_ward_score(center, size, u, t)});
    }
  }

  IntegerVector out(n);
  std::map<int, int> remap;
  int next_label = 1;

  for (int i = 0; i < n; ++i) {
    int r = sparse_ward_find_root(parent, i);
    if (remap.find(r) == remap.end()) {
      remap[r] = next_label++;
    }
    out[i] = remap[r];
  }

  return out;
}
