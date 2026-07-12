#include <Rcpp.h>
using namespace Rcpp;

/*
 Leaf scoring primitives

 The R layer keeps the scoring algebra explicit:
   1) stack per-tree leaf counts or weighted masses,
   2) convert two stacked vectors into per-leaf LLRs,
   3) score SNPs from the LLR lookup.

 These compiled helpers accelerate the two mechanical hot spots without
 reintroducing a multi-mode "do everything" backend.
 */

struct TreeSubset {
  IntegerVector tree_zero_based;
  IntegerVector offsets;
};

static TreeSubset make_tree_subset(List dense_leaf_ids,
                                   IntegerVector n_leaves,
                                   IntegerVector tree_idx) {
  const int total_trees = dense_leaf_ids.size();
  const int n_tree_subset = tree_idx.size();
  TreeSubset subset;
  subset.tree_zero_based = IntegerVector(n_tree_subset);
  subset.offsets = IntegerVector(n_tree_subset + 1, 0);

  for (int j = 0; j < n_tree_subset; ++j) {
    const int t = tree_idx[j] - 1;
    if (t < 0 || t >= total_trees) {
      stop("[leaf scoring backend] tree_idx contains an out-of-range tree.");
    }
    subset.tree_zero_based[j] = t;
    subset.offsets[j + 1] = subset.offsets[j] + n_leaves[t];
  }

  return subset;
}

// [[Rcpp::export(name = ".stack_leaf_counts_rcpp")]]
NumericVector stack_leaf_counts_rcpp(List dense_leaf_ids,
                                     IntegerVector n_leaves,
                                     IntegerVector idx,
                                     IntegerVector tree_idx) {
  TreeSubset subset = make_tree_subset(dense_leaf_ids, n_leaves, tree_idx);
  NumericVector counts(subset.offsets[tree_idx.size()], 0.0);

  for (int j = 0; j < tree_idx.size(); ++j) {
    IntegerVector dense_ids = dense_leaf_ids[subset.tree_zero_based[j]];
    const int n_train = dense_ids.size();
    const int offset = subset.offsets[j];

    for (int k = 0; k < idx.size(); ++k) {
      const int row = idx[k] - 1;
      if (row < 0 || row >= n_train) {
        stop("[.stack_leaf_counts_rcpp] training index out of range.");
      }

      const int leaf = dense_ids[row];
      if (leaf > 0) {
        counts[offset + leaf - 1] += 1.0;
      }
    }
  }

  return counts;
}

// [[Rcpp::export(name = ".stack_leaf_masses_rcpp")]]
NumericVector stack_leaf_masses_rcpp(List dense_leaf_ids,
                                     IntegerVector n_leaves,
                                     IntegerVector idx,
                                     NumericVector weights,
                                     IntegerVector tree_idx) {
  if (idx.size() != weights.size()) {
    stop("[.stack_leaf_masses_rcpp] idx and weights must have the same length.");
  }

  TreeSubset subset = make_tree_subset(dense_leaf_ids, n_leaves, tree_idx);
  NumericVector masses(subset.offsets[tree_idx.size()], 0.0);

  for (int j = 0; j < tree_idx.size(); ++j) {
    IntegerVector dense_ids = dense_leaf_ids[subset.tree_zero_based[j]];
    const int n_train = dense_ids.size();
    const int offset = subset.offsets[j];

    for (int k = 0; k < idx.size(); ++k) {
      const int row = idx[k] - 1;
      if (row < 0 || row >= n_train) {
        stop("[.stack_leaf_masses_rcpp] training index out of range.");
      }

      const double w = weights[k];
      if (w == 0.0) {
        continue;
      }

      const int leaf = dense_ids[row];
      if (leaf > 0) {
        masses[offset + leaf - 1] += w;
      }
    }
  }

  return masses;
}

// [[Rcpp::export(name = ".leaf_llrs_rcpp")]]
List leaf_llrs_rcpp(NumericVector ce_all,
                                NumericVector cb_all,
                                IntegerVector n_leaves,
                                IntegerVector tree_idx,
                                double N_E,
                                double N_B,
                                double alpha = 0.5) {
  const int n_tree_subset = tree_idx.size();
  IntegerVector offsets(n_tree_subset + 1, 0);

  for (int j = 0; j < n_tree_subset; ++j) {
    const int t = tree_idx[j] - 1;
    if (t < 0 || t >= n_leaves.size()) {
      stop("[.leaf_llrs_rcpp] tree_idx contains an out-of-range tree.");
    }
    offsets[j + 1] = offsets[j] + n_leaves[t];
  }

  const int total_leaves = offsets[n_tree_subset];
  if (ce_all.size() != total_leaves || cb_all.size() != total_leaves) {
    stop("[.leaf_llrs_rcpp] leaf evidence vectors have the wrong length for tree_idx.");
  }

  const double eps = 1e-12;
  List leaf_llrs_by_tree(n_tree_subset);

  for (int j = 0; j < n_tree_subset; ++j) {
    const int L = n_leaves[tree_idx[j] - 1];
    const int offset = offsets[j];
    NumericVector llrs(L, NA_REAL);

    if (alpha == 0.0) {
      for (int leaf = 0; leaf < L; ++leaf) {
        const double ce = ce_all[offset + leaf];
        const double cb = cb_all[offset + leaf];

        if ((ce + cb) > 0.0) {
          double pE = ce / N_E;
          double pB = cb / N_B;
          if (pE < eps) pE = eps;
          if (pB < eps) pB = eps;
          llrs[leaf] = std::log(pE / pB);
        }
      }
    } else {
      const double denom_E = N_E + alpha * static_cast<double>(L);
      const double denom_B = N_B + alpha * static_cast<double>(L);

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

  return List::create(
    _["leaf_llrs_by_tree"] = leaf_llrs_by_tree,
    _["tree_idx"] = tree_idx
  );
}
