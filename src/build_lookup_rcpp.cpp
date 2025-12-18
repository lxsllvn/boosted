// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// split(idx, leaf_native[idx], drop=TRUE) for integer groups
static inline List split_idx_by_leaf_cpp(const IntegerVector& idx,
                                         const IntegerVector& leaf_native) {
  const int n = leaf_native.size();
  const int k = idx.size();
  if (k == 0) return List::create();

  std::unordered_map<int, int> counts;
  counts.reserve((size_t)k * 2);

  std::vector<int> groups;
  groups.reserve((size_t)k);

  // Pass 1: count per leaf id + collect seen groups (ignore NA groups)
  for (int ii = 0; ii < k; ++ii) {
    int j = idx[ii];
    if (j == NA_INTEGER) continue;
    if (j < 1 || j > n) continue;

    int g = leaf_native[j - 1];
    if (g == NA_INTEGER) continue;

    auto it = counts.find(g);
    if (it == counts.end()) {
      counts.emplace(g, 1);
      groups.push_back(g);
    } else {
      it->second += 1;
    }
  }

  if (groups.empty()) return List::create();

  // Deterministic order by leaf id
  std::sort(groups.begin(), groups.end());
  groups.erase(std::unique(groups.begin(), groups.end()), groups.end());

  const int G = (int)groups.size();
  List out(G);

  // leaf id -> list position
  std::unordered_map<int, int> pos;
  pos.reserve((size_t)G * 2);
  for (int gi = 0; gi < G; ++gi) pos[groups[gi]] = gi;

  std::vector<int> cursor((size_t)G, 0);

  // allocate exact-sized buckets
  for (int gi = 0; gi < G; ++gi) {
    int g = groups[gi];
    out[gi] = IntegerVector(counts[g]);
  }

  // Pass 2: fill buckets
  for (int ii = 0; ii < k; ++ii) {
    int j = idx[ii];
    if (j == NA_INTEGER) continue;
    if (j < 1 || j > n) continue;

    int g = leaf_native[j - 1];
    if (g == NA_INTEGER) continue;

    int gi = pos[g];
    IntegerVector bucket = out[gi];
    bucket[cursor[gi]++] = j;
    out[gi] = bucket;
  }

  // names = as.character(leaf id)
  CharacterVector nms(G);
  for (int gi = 0; gi < G; ++gi) nms[gi] = std::to_string(groups[gi]);
  out.attr("names") = nms;

  return out;
}

// [[Rcpp::export(name = ".build_lookup_rcpp")]]
List build_lookup_rcpp(List dense_leaf_ids,
                       List native_leaf_ids,
                       SEXP extr_idx_sexp = R_NilValue,
                       SEXP bg_idx_sexp   = R_NilValue,
                       bool all = false) {

  const int Tm = dense_leaf_ids.size();

  const bool have_ext = !Rf_isNull(extr_idx_sexp) && Rf_length(extr_idx_sexp) > 0;
  const bool have_bg  = !Rf_isNull(bg_idx_sexp)   && Rf_length(bg_idx_sexp)   > 0;
  const bool have_all = all;

  List snps_ext_by_leaf = have_ext ? List(Tm) : List::create();
  List snps_bg_by_leaf  = have_bg  ? List(Tm) : List::create();

  // Keep as a real List when all==true; otherwise return NULL
  List snps_all_by_leaf;
  if (have_all) snps_all_by_leaf = List(Tm);

  IntegerVector extr_idx;
  IntegerVector bg_idx;
  if (have_ext) extr_idx = as<IntegerVector>(extr_idx_sexp);
  if (have_bg)  bg_idx   = as<IntegerVector>(bg_idx_sexp);

  for (int tt = 0; tt < Tm; ++tt) {
    IntegerVector inv_t      = dense_leaf_ids[tt];
    IntegerVector native_ids = native_leaf_ids[tt];

    if (inv_t.size() == 0 || native_ids.size() == 0) {
      if (have_ext) snps_ext_by_leaf[tt] = List::create();
      if (have_bg)  snps_bg_by_leaf[tt]  = List::create();
      if (have_all) snps_all_by_leaf[tt] = List::create();
      continue;
    }

    // Defensive: ensure inv_t indexes native_ids (like R would error if invalid)
    const int n_leaves = native_ids.size();
    for (int i = 0; i < inv_t.size(); ++i) {
      int di = inv_t[i];
      if (di == NA_INTEGER || di < 1 || di > n_leaves) {
        stop(".build_lookup_rcpp: dense_leaf_ids[[%d]] contains invalid index (must be 1..%d).",
             tt + 1, n_leaves);
      }
    }

    // leaf_native = native_ids[inv_t]
    IntegerVector leaf_native(inv_t.size());
    for (int i = 0; i < inv_t.size(); ++i) {
      leaf_native[i] = native_ids[inv_t[i] - 1];
    }

    if (have_ext) snps_ext_by_leaf[tt] = split_idx_by_leaf_cpp(extr_idx, leaf_native);
    if (have_bg)  snps_bg_by_leaf[tt]  = split_idx_by_leaf_cpp(bg_idx,  leaf_native);

    if (have_all) {
      IntegerVector all_idx = seq_len(leaf_native.size()); // 1..n
      snps_all_by_leaf[tt] = split_idx_by_leaf_cpp(all_idx, leaf_native);
    }
  }

  return List::create(
    _["snps_ext_by_leaf"] = (have_ext ? (SEXP)snps_ext_by_leaf : R_NilValue),
    _["snps_bg_by_leaf"]  = (have_bg  ? (SEXP)snps_bg_by_leaf  : R_NilValue),
    _["snps_all_by_leaf"] = (have_all ? (SEXP)snps_all_by_leaf : R_NilValue)
  );
}
