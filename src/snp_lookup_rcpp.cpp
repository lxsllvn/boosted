// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>

static inline Rcpp::IntegerVector dedup_first_seen(
    const std::vector<int>& x) {

  if (x.empty()) return Rcpp::IntegerVector(0);

  std::unordered_set<int> seen;
  seen.reserve(x.size() * 2);

  std::vector<int> out;
  out.reserve(x.size());

  for (int v : x) {
    if (seen.insert(v).second) {
      out.push_back(v);
    }
  }

  Rcpp::IntegerVector res(out.size());
  std::copy(out.begin(), out.end(), res.begin());
  return res;
}


using namespace Rcpp;

// [[Rcpp::export(name = ".snp_lookup_rcpp")]]
SEXP snp_lookup_rcpp(DataFrame pairs,
                     List snps_ext_by_leaf,
                     List snps_bg_by_leaf,
                     SEXP snps_all_by_leaf = R_NilValue) {

  int n_pairs = pairs.nrows();
  if (n_pairs == 0) {
    return List::create(
      _["bucket_ext"] = IntegerVector(0),
      _["bucket_bg"]  = IntegerVector(0),
      _["bucket_all"] = IntegerVector(0)
    );
  }

  IntegerVector tree_col = pairs["Tree"];     // 0-based
  IntegerVector leaf_col = pairs["leaf_id"];  // native leaf id
  int Tm = snps_ext_by_leaf.size();

  bool have_all = !Rf_isNull(snps_all_by_leaf);
  List all_list;
  if (have_all) all_list = List(snps_all_by_leaf);

  std::vector<int> out_ext;
  std::vector<int> out_bg;
  std::vector<int> out_all;

  for (int j = 0; j < n_pairs; ++j) {
    int tr0  = tree_col[j];
    int leaf = leaf_col[j];
    if (tr0 < 0 || tr0 >= Tm) continue;

    std::string lid = std::to_string(leaf);

    // extremes
    {
      List ext_tree = snps_ext_by_leaf[tr0];
      if (ext_tree.containsElementNamed(lid.c_str())) {
        IntegerVector v = ext_tree[lid];
        out_ext.insert(out_ext.end(), v.begin(), v.end());
      }
    }

    // background
    {
      List bg_tree = snps_bg_by_leaf[tr0];
      if (bg_tree.containsElementNamed(lid.c_str())) {
        IntegerVector v = bg_tree[lid];
        out_bg.insert(out_bg.end(), v.begin(), v.end());
      }
    }

    // all
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
