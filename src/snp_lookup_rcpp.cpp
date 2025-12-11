// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
using namespace Rcpp;

/*
 .snp_lookup_rcpp

 R equivalent:

 .snp_lookup_R <- function(pairs,
 snps_ext_by_leaf,
 snps_bg_by_leaf,
 snps_all_by_leaf = NULL) { ... }

 Inputs:
 pairs: data.frame/data.table with columns:
 - Tree    : integer, 0-based tree index (0..Tm-1)
 - leaf_id : integer, native leaf id

 snps_ext_by_leaf, snps_bg_by_leaf, snps_all_by_leaf:
 List of length Tm; each element is a named list where
 names are leaf_id as character (e.g. "17") and each entry
 is an integer vector of SNP indices.

 Output:
 list(
 bucket_ext = integer(),
 bucket_bg  = integer(),
 bucket_all = integer()
 )
 where each vector is deduplicated.
 */

// helper: build a map from leaf_id string -> index in the per-tree list
static std::unordered_map<std::string,int>
  build_leaf_index_map(const List &leaf_list) {
    std::unordered_map<std::string,int> m;
    CharacterVector nm = leaf_list.names();
    int n = nm.size();
    for (int i = 0; i < n; ++i) {
      if (nm[i] == NA_STRING) continue;
      std::string key = Rcpp::as<std::string>(nm[i]);
      // store position i (0-based for C++)
      m.emplace(key, i);
    }
    return m;
  }

// [[Rcpp::export(name = ".snp_lookup_rcpp")]]
SEXP snp_lookup_rcpp(DataFrame pairs,
                     List snps_ext_by_leaf,
                     List snps_bg_by_leaf,
                     SEXP snps_all_by_leaf = R_NilValue) {

  // empty output if no pairs
  int n_pairs = pairs.nrows();
  if (n_pairs == 0) {
    return List::create(
      _["bucket_ext"] = IntegerVector(0),
      _["bucket_bg"]  = IntegerVector(0),
      _["bucket_all"] = IntegerVector(0)
    );
  }

  IntegerVector tree_col = pairs["Tree"];     // 0-based tree index
  IntegerVector leaf_col = pairs["leaf_id"];  // native leaf id

  int Tm = snps_ext_by_leaf.size();

  // caches: for each tree, a map leaf_id_string -> index in that tree's list
  std::vector< std::unordered_map<std::string,int> > index_maps(Tm);
  std::vector<bool> map_built(Tm, false);

  bool have_all = !Rf_isNull(snps_all_by_leaf);
  List all_list;
  if (have_all) {
    all_list = List(snps_all_by_leaf);
  }

  // use sets to deduplicate SNP indices as we go
  std::unordered_set<int> set_ext;
  std::unordered_set<int> set_bg;
  std::unordered_set<int> set_all;

  // main loop over (Tree, leaf_id) pairs
  for (int j = 0; j < n_pairs; ++j) {
    int tr0   = tree_col[j];  // 0..Tm-1
    int leaf  = leaf_col[j];

    if (tr0 < 0 || tr0 >= Tm) continue;

    // lazy-build index map for this tree
    if (!map_built[tr0]) {
      List ext_tree = snps_ext_by_leaf[tr0];
      index_maps[tr0] = build_leaf_index_map(ext_tree);
      map_built[tr0]  = true;
    }

    std::string lid = std::to_string(leaf);
    auto &m = index_maps[tr0];
    auto it = m.find(lid);
    if (it == m.end()) {
      // no such leaf_id for this tree
      continue;
    }
    int pos = it->second; // index in per-tree list

    // extremes
    {
      List ext_tree = snps_ext_by_leaf[tr0];
      if (pos >= 0 && pos < ext_tree.size()) {
        SEXP re = ext_tree[pos];
        if (re != R_NilValue) {
          IntegerVector v(re);
          for (int k = 0; k < v.size(); ++k) {
            // SNP indices are 1-based; we keep them as-is
            set_ext.insert(v[k]);
          }
        }
      }
    }

    // background
    {
      List bg_tree = snps_bg_by_leaf[tr0];
      if (pos >= 0 && pos < bg_tree.size()) {
        SEXP rb = bg_tree[pos];
        if (rb != R_NilValue) {
          IntegerVector v(rb);
          for (int k = 0; k < v.size(); ++k) {
            set_bg.insert(v[k]);
          }
        }
      }
    }

    // all SNPs (optional)
    if (have_all) {
      List all_tree = all_list[tr0];
      if (pos >= 0 && pos < all_tree.size()) {
        SEXP ra = all_tree[pos];
        if (ra != R_NilValue) {
          IntegerVector v(ra);
          for (int k = 0; k < v.size(); ++k) {
            set_all.insert(v[k]);
          }
        }
      }
    }
  }

  // convert sets back to integer vectors
  auto set_to_iv = [](const std::unordered_set<int> &s) {
    IntegerVector out(s.size());
    int i = 0;
    for (auto v : s) {
      out[i++] = v;
    }
    return out;
  };

  IntegerVector bucket_ext = set_to_iv(set_ext);
  IntegerVector bucket_bg  = set_to_iv(set_bg);
  IntegerVector bucket_all = set_to_iv(set_all);

  return List::create(
    _["bucket_ext"] = bucket_ext,
    _["bucket_bg"]  = bucket_bg,
    _["bucket_all"] = bucket_all
  );
}
