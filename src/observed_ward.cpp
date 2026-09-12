#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

const double UNIT_ROUNDOFF = 0x1p-53;

double gamma_n(int k) {
  if (k <= 0) return 0.0;
  double ku = static_cast<double>(k) * UNIT_ROUNDOFF;
  return ku < 1.0 ? ku / (1.0 - ku) : R_PosInf;
}

struct Candidate {
  double cost;
  double bound;
  int key1;
  int key2;
  int proposal;
  int node_a;
  int node_b;
  bool exact_demonstrated;
};

bool canonical_less(const Candidate& a, const Candidate& b) {
  return a.key1 < b.key1 || (a.key1 == b.key1 && a.key2 < b.key2);
}

bool strict_less(const Candidate& a, const Candidate& b) {
  if (a.cost < b.cost) return true;
  if (a.cost > b.cost) return false;
  return canonical_less(a, b);
}

bool equivalent(const Candidate& a, const Candidate& b) {
  return std::fabs(a.cost - b.cost) <= a.bound + b.bound;
}

int root(std::vector<int>& parent, int x) {
  int r = x;
  while (parent[r] != r) r = parent[r];
  while (parent[x] != x) {
    int p = parent[x]; parent[x] = r; x = p;
  }
  return r;
}

} // namespace

// [[Rcpp::export]]
List capivara_observed_ward_cpp(NumericMatrix features,
                                LogicalMatrix validity,
                                IntegerMatrix edges,
                                IntegerVector original_ids,
                                IntegerVector feature_group,
                                double min_shared_fraction,
                                double min_contributor_fraction,
                                double min_feature_fraction,
                                double input_error_gamma,
                                CharacterVector leaf_signatures,
                                int target_k) {
  const int n = features.nrow();
  const int p = features.ncol();
  if (n < 1 || p < 1 || validity.nrow() != n || validity.ncol() != p) {
    stop("Features and validity must be non-empty matrices with equal dimensions.");
  }
  if (original_ids.size() != n || leaf_signatures.size() != n || feature_group.size() != p) {
    stop("Original IDs, signatures and feature groups do not match the features.");
  }
  if (target_k < 1 || target_k > n) stop("target_k must be between 1 and the number of leaves.");
  if (!R_finite(min_shared_fraction) || !R_finite(min_contributor_fraction) ||
      !R_finite(min_feature_fraction) || min_shared_fraction < 0.0 ||
      min_shared_fraction > 1.0 || min_contributor_fraction < 0.0 ||
      min_contributor_fraction > 1.0 || min_feature_fraction < 0.0 ||
      min_feature_fraction > 1.0) {
    stop("Overlap and contributor fractions must lie in [0, 1].");
  }
  if (!R_finite(input_error_gamma) || input_error_gamma < 0.0) {
    stop("The input forward-error gamma must be finite and non-negative.");
  }
  std::set<int> unique_ids;
  for (int i = 0; i < n; ++i) {
    if (original_ids[i] == NA_INTEGER || original_ids[i] < 1 ||
        !unique_ids.insert(original_ids[i]).second) {
      stop("Original spaxel IDs must be unique positive integers.");
    }
  }
  int groups = 0;
  for (int j = 0; j < p; ++j) {
    if (feature_group[j] < 0) stop("Feature-group identifiers cannot be negative.");
    groups = std::max(groups, feature_group[j]);
  }
  std::vector<int> group_total(groups + 1, 0);
  for (int j = 0; j < p; ++j) if (feature_group[j] > 0) ++group_total[feature_group[j]];

  std::vector<int> counts(static_cast<size_t>(n) * p, 0);
  std::vector<double> sums(static_cast<size_t>(n) * p, 0.0);
  std::vector<double> sumabs(static_cast<size_t>(n) * p, 0.0);
  std::vector<int> size(n, 1), parent(n), min_id(n), node_id(n);
  std::vector<char> active(n, 1), signature_known(n, 1);
  std::vector<std::string> signature(n);
  std::vector< std::set<int> > neighbours(n);
  for (int i = 0; i < n; ++i) {
    parent[i] = i; min_id[i] = original_ids[i]; node_id[i] = i + 1;
    signature[i] = as<std::string>(leaf_signatures[i]);
    if (signature[i].empty()) signature_known[i] = 0;
    for (int j = 0; j < p; ++j) {
      const size_t at = static_cast<size_t>(i) * p + j;
      const bool ok = validity(i, j) == TRUE;
      if (ok && !R_finite(features(i, j))) stop("A valid feature entry is non-finite.");
      if (ok) { counts[at] = 1; sums[at] = features(i, j); sumabs[at] = std::fabs(features(i, j)); }
    }
  }
  for (int e = 0; e < edges.nrow(); ++e) {
    int a = edges(e, 0) - 1, b = edges(e, 1) - 1;
    if (a < 0 || b < 0 || a >= n || b >= n || a == b) stop("Invalid spatial edge.");
    neighbours[a].insert(b); neighbours[b].insert(a);
  }

  std::map< std::pair<int,int>, Candidate > candidates;
  int proposals = 0, unsupported = 0, admitted = 0;
  auto propose = [&](int aa, int bb) {
    int a = std::min(aa, bb), b = std::max(aa, bb);
    std::pair<int,int> pair(a, b);
    int reliable_count = 0, common_count = 0;
    std::vector<int> reliable_group(groups + 1, 0);
    for (int j = 0; j < p; ++j) {
      const int ca = counts[static_cast<size_t>(a) * p + j];
      const int cb = counts[static_cast<size_t>(b) * p + j];
      if (ca > 0 && cb > 0) ++common_count;
      const bool reliable = ca > 0 && cb > 0 &&
        static_cast<double>(ca) >= min_contributor_fraction * size[a] &&
        static_cast<double>(cb) >= min_contributor_fraction * size[b];
      if (reliable) {
        ++reliable_count;
        if (feature_group[j] > 0) ++reliable_group[feature_group[j]];
      }
    }
    bool allowed = common_count > 0 &&
      static_cast<double>(reliable_count) / p >= min_shared_fraction;
    for (int g = 1; allowed && g <= groups; ++g) {
      allowed = group_total[g] > 0 &&
        static_cast<double>(reliable_group[g]) / group_total[g] >= min_feature_fraction;
    }
    if (!allowed) { candidates.erase(pair); ++unsupported; return; }

    double cost = 0.0, bound = 0.0, absolute_cost = 0.0;
    int represented = 0;
    for (int j = 0; j < p; ++j) {
      const size_t ia = static_cast<size_t>(a) * p + j;
      const size_t ib = static_cast<size_t>(b) * p + j;
      const int ca = counts[ia], cb = counts[ib];
      if (ca <= 0 || cb <= 0) continue;
      ++represented;
      const double ma = sums[ia] / ca, mb = sums[ib] / cb;
      const double ema = (input_error_gamma * sumabs[ia] + gamma_n(ca - 1) * sumabs[ia] + UNIT_ROUNDOFF * std::fabs(sums[ia])) / ca;
      const double emb = (input_error_gamma * sumabs[ib] + gamma_n(cb - 1) * sumabs[ib] + UNIT_ROUNDOFF * std::fabs(sums[ib])) / cb;
      const double d = ma - mb;
      const double ed = ema + emb + UNIT_ROUNDOFF * (std::fabs(ma) + std::fabs(mb));
      const double coefficient = static_cast<double>(ca) * cb / (ca + cb);
      const double term = coefficient * d * d;
      const double term_bound = coefficient * (2.0 * std::fabs(d) * ed + ed * ed) +
        gamma_n(4) * coefficient * std::pow(std::fabs(d) + ed, 2.0);
      if (!R_finite(term) || !R_finite(term_bound)) {
        stop("A Ward merge cost or its binary64 bound is non-finite.");
      }
      cost += term; absolute_cost += std::fabs(term); bound += term_bound;
    }
    bound += gamma_n(represented - 1) * absolute_cost;
    bound = std::max(bound, std::nextafter(0.0, 1.0));
    Candidate candidate{cost, bound, std::min(min_id[a], min_id[b]),
                        std::max(min_id[a], min_id[b]), proposals,
                        node_id[a], node_id[b],
                        cost == 0.0 && signature_known[a] && signature_known[b] && signature[a] == signature[b]};
    candidates[pair] = candidate; ++proposals; ++admitted;
  };
  for (int i = 0; i < n; ++i) for (int j : neighbours[i]) if (i < j) propose(i, j);

  std::vector<int> child_a, child_b, selected_key_a, selected_key_b, tie_group_size;
  std::vector<double> costs, bounds, strict_costs, maximum_group_cost, maximum_group_bound;
  std::vector<int> tie_steps, tie_exact, tie_strict_a, tie_strict_b, tie_selected_a, tie_selected_b;
  std::set<int> tie_proposal_ids;
  int active_count = n, exact_groups = 0, numerical_groups = 0, selected_not_strict = 0;
  IntegerVector target_labels;

  auto capture_labels = [&]() {
    IntegerVector labels(n); std::map<int,int> ordered_roots;
    for (int i = 0; i < n; ++i) {
      int r = root(parent, i); ordered_roots[min_id[r]] = r;
    }
    std::map<int,int> label_for_slot; int label = 1;
    for (const auto& item : ordered_roots) label_for_slot[item.second] = label++;
    for (int i = 0; i < n; ++i) labels[i] = label_for_slot[root(parent, i)];
    return labels;
  };
  if (active_count == target_k) target_labels = capture_labels();

  while (!candidates.empty()) {
    auto strict_it = candidates.begin();
    for (auto it = candidates.begin(); it != candidates.end(); ++it)
      if (strict_less(it->second, strict_it->second)) strict_it = it;
    const Candidate strict = strict_it->second;
    std::vector< std::map<std::pair<int,int>, Candidate>::iterator > group;
    for (auto it = candidates.begin(); it != candidates.end(); ++it)
      if (equivalent(it->second, strict)) group.push_back(it);
    auto selected_it = group.front();
    for (auto it : group) if (canonical_less(it->second, selected_it->second)) selected_it = it;
    const Candidate selected = selected_it->second;
    for (auto it : group) tie_proposal_ids.insert(it->second.proposal);
    if (group.size() > 1) {
      bool exact = true; double max_cost = -R_PosInf, max_bound = 0.0;
      for (auto it : group) {
        exact = exact && it->second.exact_demonstrated;
        max_cost = std::max(max_cost, it->second.cost);
        max_bound = std::max(max_bound, it->second.bound);
      }
      exact_groups += exact; numerical_groups += !exact;
      selected_not_strict += selected.proposal != strict.proposal;
      tie_steps.push_back(costs.size() + 1); tie_group_size.push_back(group.size()); tie_exact.push_back(exact);
      tie_strict_a.push_back(strict.key1); tie_strict_b.push_back(strict.key2);
      tie_selected_a.push_back(selected.key1); tie_selected_b.push_back(selected.key2);
      strict_costs.push_back(strict.cost); maximum_group_cost.push_back(max_cost); maximum_group_bound.push_back(max_bound);
    }
    int a = selected_it->first.first, b = selected_it->first.second;
    if (min_id[b] < min_id[a]) std::swap(a, b);
    const int kept = std::min(a, b), removed = std::max(a, b);
    const int left_node = min_id[a] <= min_id[b] ? node_id[a] : node_id[b];
    const int right_node = min_id[a] <= min_id[b] ? node_id[b] : node_id[a];
    child_a.push_back(left_node); child_b.push_back(right_node);
    costs.push_back(selected.cost); bounds.push_back(selected.bound);
    selected_key_a.push_back(selected.key1); selected_key_b.push_back(selected.key2);

    std::set<int> targets = neighbours[a];
    targets.insert(neighbours[b].begin(), neighbours[b].end());
    targets.erase(a); targets.erase(b);
    for (auto it = candidates.begin(); it != candidates.end(); ) {
      if (it->first.first == a || it->first.second == a || it->first.first == b || it->first.second == b) it = candidates.erase(it);
      else ++it;
    }
    for (int j = 0; j < p; ++j) {
      const size_t ik = static_cast<size_t>(kept) * p + j;
      const size_t ir = static_cast<size_t>(removed) * p + j;
      counts[ik] += counts[ir]; sums[ik] += sums[ir]; sumabs[ik] += sumabs[ir];
      counts[ir] = 0; sums[ir] = 0.0; sumabs[ir] = 0.0;
    }
    size[kept] = size[a] + size[b]; size[removed] = 0;
    min_id[kept] = std::min(min_id[a], min_id[b]);
    const bool same_signature = signature_known[a] && signature_known[b] && signature[a] == signature[b];
    signature[kept] = same_signature ? signature[a] : std::string();
    signature_known[kept] = same_signature; signature_known[removed] = 0;
    node_id[kept] = n + costs.size(); node_id[removed] = 0;
    active[kept] = 1; active[removed] = 0; parent[removed] = kept; parent[kept] = kept;
    neighbours[kept] = targets; neighbours[removed].clear();
    for (int t0 : targets) {
      int t = root(parent, t0);
      if (t == kept || !active[t]) continue;
      neighbours[t].erase(a); neighbours[t].erase(b); neighbours[t].insert(kept);
    }
    --active_count;
    for (int t0 : targets) {
      int t = root(parent, t0);
      if (t != kept && active[t]) propose(kept, t);
    }
    if (active_count == target_k && target_labels.size() == 0) target_labels = capture_labels();
  }
  if (target_labels.size() == 0) target_labels = capture_labels();
  int actual_k = 0;
  for (int x : target_labels) actual_k = std::max(actual_k, x);

  IntegerMatrix children(costs.size(), 2);
  for (int i = 0; i < children.nrow(); ++i) { children(i, 0) = child_a[i]; children(i, 1) = child_b[i]; }
  int inversions = 0;
  std::vector<double> node_cost(n + costs.size(), NA_REAL);
  for (size_t i = 0; i < costs.size(); ++i) {
    node_cost[n + i] = costs[i];
    for (int child : {child_a[i], child_b[i]})
      if (child > n && costs[i] < node_cost[child - 1]) ++inversions;
  }
  IntegerVector final_roots;
  std::vector<int> root_nodes;
  for (int i = 0; i < n; ++i) if (active[i]) root_nodes.push_back(node_id[i]);
  final_roots = wrap(root_nodes);

  DataFrame tie_groups = DataFrame::create(
    _["step"] = tie_steps, _["group_size"] = tie_group_size,
    _["exact_demonstrated"] = tie_exact,
    _["strict_key_a"] = tie_strict_a, _["strict_key_b"] = tie_strict_b,
    _["selected_key_a"] = tie_selected_a, _["selected_key_b"] = tie_selected_b,
    _["strict_cost"] = strict_costs, _["maximum_cost"] = maximum_group_cost,
    _["maximum_bound"] = maximum_group_bound
  );
  return List::create(
    _["children"] = children, _["cost"] = costs, _["cost_bound"] = bounds,
    _["selected_key_a"] = selected_key_a, _["selected_key_b"] = selected_key_b,
    _["labels"] = target_labels, _["requested_k"] = target_k,
    _["actual_k"] = actual_k, _["roots"] = final_roots,
    _["tie_groups"] = tie_groups,
    _["qc"] = List::create(
      _["objective"] = "observed-entry within-cluster SSE",
      _["comparability_rule"] = "separate reliable-overlap admission",
      _["tie_rule"] = "overlapping binary64 forward-error intervals; lexicographic minimum original IDs",
      _["unit_roundoff"] = UNIT_ROUNDOFF, _["proposals"] = proposals,
      _["admitted_proposals"] = admitted, _["unsupported_comparisons"] = unsupported,
      _["tie_proposals"] = tie_proposal_ids.size(),
      _["decisions_with_ties"] = tie_steps.size(),
      _["exact_demonstrated_tie_groups"] = exact_groups,
      _["numerical_tie_groups"] = numerical_groups,
      _["equivalent_selected_not_strict_minimum"] = selected_not_strict,
      _["raw_parent_child_inversions"] = inversions,
      _["imputed_samples"] = 0, _["costs_modified"] = false,
      _["height_monotonicization"] = false,
      _["spatial_connectivity"] = "four-neighbour supplied edges"
    )
  );
}
