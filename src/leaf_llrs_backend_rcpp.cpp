#include <Rcpp.h>
using namespace Rcpp;

/*
 .leaf_llrs_backend_rcpp

 Unified selected-index backend for resampling code. The dense leaf map is the
 same object used by the readable R reference implementation:

 dense_leaf_ids[[t]][i] = dense leaf ID (1..L_t) for training SNP i in tree t
 n_leaves[t]            = number of dense leaves in tree t

 The backend supports four modes:
 1) explicit extr_idx + bg_idx
 2) explicit bg_idx + fixed_ce_all
 3) explicit extr_idx + fixed_cb_all
 4) fixed relabelling pool + exactly one varying side

 It can also return stacked counts directly for one explicit index set, which
 lets the R resampling code precompute fixed counts once and keep the hot loop
 conceptually simple.
 */

// [[Rcpp::export(name = ".leaf_llrs_backend_rcpp")]]
SEXP leaf_llrs_backend_rcpp(List dense_leaf_ids,
                            IntegerVector n_leaves,
                            int N_extr,
                            int N_bg,
                            double alpha,
                            SEXP extr_idx_sexp = R_NilValue,
                            SEXP bg_idx_sexp = R_NilValue,
                            SEXP fixed_ce_all_sexp = R_NilValue,
                            SEXP fixed_cb_all_sexp = R_NilValue,
                            SEXP pool_counts_all_sexp = R_NilValue,
                            SEXP tree_idx_sexp = R_NilValue,
                            bool return_counts = false) {
  const int total_trees = dense_leaf_ids.size();
  IntegerVector tree_idx;

  if (Rf_isNull(tree_idx_sexp)) {
    tree_idx = seq_len(total_trees);
  } else {
    tree_idx = as<IntegerVector>(tree_idx_sexp);
  }

  const int n_tree_subset = tree_idx.size();
  IntegerVector tree_zero_based(n_tree_subset);
  IntegerVector offsets(n_tree_subset + 1, 0);

  for (int j = 0; j < n_tree_subset; ++j) {
    const int t = tree_idx[j] - 1;
    if (t < 0 || t >= total_trees) {
      stop("[.leaf_llrs_backend_rcpp] tree_idx contains an out-of-range tree.");
    }
    tree_zero_based[j] = t;
    offsets[j + 1] = offsets[j] + n_leaves[t];
  }

  const int total_leaves = offsets[n_tree_subset];
  const bool has_extr = !Rf_isNull(extr_idx_sexp);
  const bool has_bg = !Rf_isNull(bg_idx_sexp);
  const bool has_fixed_ce = !Rf_isNull(fixed_ce_all_sexp);
  const bool has_fixed_cb = !Rf_isNull(fixed_cb_all_sexp);
  const bool has_pool = !Rf_isNull(pool_counts_all_sexp);

  auto count_selected = [&](SEXP idx_sexp) {
    IntegerVector idx = as<IntegerVector>(idx_sexp);
    NumericVector counts(total_leaves, 0.0);

    for (int j = 0; j < n_tree_subset; ++j) {
      IntegerVector dense_ids = dense_leaf_ids[tree_zero_based[j]];
      const int n_train = dense_ids.size();
      const int offset = offsets[j];

      for (int k = 0; k < idx.size(); ++k) {
        const int row = idx[k] - 1;
        if (row < 0 || row >= n_train) {
          stop("[.leaf_llrs_backend_rcpp] Training index out of range.");
        }

        const int leaf = dense_ids[row];
        if (leaf > 0) {
          counts[offset + leaf - 1] += 1.0;
        }
      }
    }

    return counts;
  };

  if (return_counts) {
    if (has_extr == has_bg) {
      stop("[.leaf_llrs_backend_rcpp] return_counts requires exactly one explicit index vector.");
    }
    return count_selected(has_extr ? extr_idx_sexp : bg_idx_sexp);
  }

  NumericVector ce_all(total_leaves, 0.0);
  NumericVector cb_all(total_leaves, 0.0);

  if (has_fixed_ce) {
    NumericVector fixed_ce = as<NumericVector>(fixed_ce_all_sexp);
    if (fixed_ce.size() != total_leaves) {
      stop("[.leaf_llrs_backend_rcpp] fixed_ce_all has the wrong length for tree_idx.");
    }
    ce_all = clone(fixed_ce);
    cb_all = count_selected(bg_idx_sexp);
  } else if (has_fixed_cb) {
    NumericVector fixed_cb = as<NumericVector>(fixed_cb_all_sexp);
    if (fixed_cb.size() != total_leaves) {
      stop("[.leaf_llrs_backend_rcpp] fixed_cb_all has the wrong length for tree_idx.");
    }
    ce_all = count_selected(extr_idx_sexp);
    cb_all = clone(fixed_cb);
  } else if (has_pool) {
    NumericVector pool_counts = as<NumericVector>(pool_counts_all_sexp);
    if (pool_counts.size() != total_leaves) {
      stop("[.leaf_llrs_backend_rcpp] pool_counts_all has the wrong length for tree_idx.");
    }

    if (has_extr && !has_bg) {
      ce_all = count_selected(extr_idx_sexp);
      cb_all = clone(pool_counts);
      for (int i = 0; i < total_leaves; ++i) {
        cb_all[i] -= ce_all[i];
      }
    } else if (has_bg && !has_extr) {
      cb_all = count_selected(bg_idx_sexp);
      ce_all = clone(pool_counts);
      for (int i = 0; i < total_leaves; ++i) {
        ce_all[i] -= cb_all[i];
      }
    } else {
      stop("[.leaf_llrs_backend_rcpp] pool_counts_all requires exactly one explicit side.");
    }
  } else {
    if (!has_extr || !has_bg) {
      stop("[.leaf_llrs_backend_rcpp] explicit mode requires extr_idx and bg_idx.");
    }
    ce_all = count_selected(extr_idx_sexp);
    cb_all = count_selected(bg_idx_sexp);
  }

  const double eps = 1e-12;
  List leaf_llrs_by_tree(n_tree_subset);

  for (int j = 0; j < n_tree_subset; ++j) {
    const int L = n_leaves[tree_zero_based[j]];
    const int offset = offsets[j];
    NumericVector llrs(L, NA_REAL);

    if (alpha == 0.0) {
      for (int leaf = 0; leaf < L; ++leaf) {
        const double ce = ce_all[offset + leaf];
        const double cb = cb_all[offset + leaf];

        if ((ce + cb) > 0.0) {
          double pE = ce / static_cast<double>(N_extr);
          double pB = cb / static_cast<double>(N_bg);
          if (pE < eps) pE = eps;
          if (pB < eps) pB = eps;
          llrs[leaf] = std::log(pE / pB);
        }
      }
    } else {
      const double denom_E = static_cast<double>(N_extr) + alpha * static_cast<double>(L);
      const double denom_B = static_cast<double>(N_bg) + alpha * static_cast<double>(L);

      for (int leaf = 0; leaf < L; ++leaf) {
        const double ce = ce_all[offset + leaf];
        const double cb = cb_all[offset + leaf];

        if ((ce + cb) > 0.0) {
          const double pE = (ce + alpha) / denom_E;
          const double pB = (cb + alpha) / denom_B;
          llrs[leaf] = std::log(pE / pB);
        }
      }
    }

    leaf_llrs_by_tree[j] = llrs;
  }

  return List::create(_["leaf_llrs_by_tree"] = leaf_llrs_by_tree);
}
