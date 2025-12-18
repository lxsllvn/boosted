// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>

using namespace Rcpp;

static inline Rcpp::IntegerVector dedup_first_seen(const std::vector<int>& x) {
  if (x.empty()) return Rcpp::IntegerVector(0);

  std::unordered_set<int> seen;
  seen.reserve(x.size() * 2);

  std::vector<int> out;
  out.reserve(x.size());

  for (int v : x) {
    if (seen.insert(v).second) out.push_back(v);
  }

  Rcpp::IntegerVector res(out.size());
  std::copy(out.begin(), out.end(), res.begin());
  return res;
}

// [[Rcpp::export(name = ".snp_lookup_rcpp")]]
SEXP snp_lookup_rcpp(DataFrame pairs,
                     SEXP snps_ext_by_leaf = R_NilValue,
                     SEXP snps_bg_by_leaf  = R_NilValue,
                     SEXP snps_all_by_leaf = R_NilValue) {

  const int n_pairs = pairs.nrows();
  if (n_pairs == 0) {
    return List::create(
      _["bucket_ext"] = IntegerVector(0),
      _["bucket_bg"]  = IntegerVector(0),
      _["bucket_all"] = IntegerVector(0)
    );
  }

  const bool have_ext = !Rf_isNull(snps_ext_by_leaf);
  const bool have_bg  = !Rf_isNull(snps_bg_by_leaf);
  const bool have_all = !Rf_isNull(snps_all_by_leaf);

  if (!have_ext && !have_bg && !have_all) {
    return List::create(
      _["bucket_ext"] = IntegerVector(0),
      _["bucket_bg"]  = IntegerVector(0),
      _["bucket_all"] = IntegerVector(0)
    );
  }

  List ext_list, bg_list, all_list;
  int Tm = 0;

  if (have_ext) {
    ext_list = as<List>(snps_ext_by_leaf);
    Tm = ext_list.size();
  } else if (have_bg) {
    bg_list = as<List>(snps_bg_by_leaf);
    Tm = bg_list.size();
  } else {
    all_list = as<List>(snps_all_by_leaf);
    Tm = all_list.size();
  }

  if (have_bg && bg_list.size() == 0) bg_list = as<List>(snps_bg_by_leaf);
  if (have_all && all_list.size() == 0) all_list = as<List>(snps_all_by_leaf);

  IntegerVector tree_col = pairs["Tree"];     // 0-based
  IntegerVector leaf_col = pairs["leaf_id"];  // native leaf id (int)

  std::vector<int> out_ext;
  std::vector<int> out_bg;
  std::vector<int> out_all;

  for (int j = 0; j < n_pairs; ++j) {
    const int tr0  = tree_col[j];
    const int leaf = leaf_col[j];
    if (tr0 < 0 || tr0 >= Tm) continue;

    const std::string lid = std::to_string(leaf);

    if (have_ext) {
      List ext_tree = ext_list[tr0];
      if (ext_tree.containsElementNamed(lid.c_str())) {
        IntegerVector v = ext_tree[lid];
        out_ext.insert(out_ext.end(), v.begin(), v.end());
      }
    }

    if (have_bg) {
      List bg_tree = bg_list[tr0];
      if (bg_tree.containsElementNamed(lid.c_str())) {
        IntegerVector v = bg_tree[lid];
        out_bg.insert(out_bg.end(), v.begin(), v.end());
      }
    }

    if (have_all) {
      List all_tree = all_list[tr0];
      if (all_tree.containsElementNamed(lid.c_str())) {
        IntegerVector v = all_tree[lid];
        out_all.insert(out_all.end(), v.begin(), v.end());
      }
    }
  }

  return List::create(
    _["bucket_ext"] = dedup_first_seen(out_ext),
    _["bucket_bg"]  = dedup_first_seen(out_bg),
    _["bucket_all"] = dedup_first_seen(out_all)
  );
}
