#include <Rcpp.h>
using namespace Rcpp;

/*
 .score_snps_rcpp

 Inputs:
 test_leaf_map: list with element "dense_leaf_ids", a list<IntegerVector>
 of length Tm. Each IntegerVector has length n and contains:
 0      = sentinel ("no TRAIN leaf")
 1..L_t = dense TRAIN leaf ID

 leaf_llrs_by_tree:
 $leaf_llrs_by_tree : list<NumericVector> length Tm
 (log(pE/pB) per leaf)

 Tm : number of trees
 n  : number of SNPs

 Output:
 List with numeric vector of length n:
 $scores    # mean leaf log likelihood ratio (NA if no leaf evidence)
 */

// [[Rcpp::export(name = ".score_snps_rcpp")]]
SEXP score_snps_rcpp(SEXP test_leaf_map,
                     SEXP leaf_llrs_by_tree,
                     int Tm,
                     int n) {

  // ----- Extract test dense leaf IDs -----
  List tlm(test_leaf_map);
  List pos_list = tlm["dense_leaf_ids"];   // list of IntegerVector, length Tm

  // ----- Extract LLR inputs  -----
  List le(leaf_llrs_by_tree);
  List leaf_llrs_list = le["leaf_llrs_by_tree"];  // list<NumericVector>

  // ----- Accumulators -----
  NumericVector llrs_sum(n, 0.0);
  IntegerVector used_llrs(n, 0);

  // ----- Main loop: iterate trees -----
  for (int t = 0; t < Tm; ++t) {

    IntegerVector pos_t = pos_list[t];
    if (pos_t.size() == 0) continue;   // empty tree for TEST

    NumericVector con_t = leaf_llrs_list[t];

    // Assume pos_t.size() == n
    for (int i = 0; i < n; ++i) {
      int p = pos_t[i];               // 0 = sentinel; 1..L = dense ID
      if (p == 0) continue;

      int idx = p - 1;                // convert to 0-based index

      double cv = con_t[idx];
      if (!NumericVector::is_na(cv)) {
        llrs_sum[i]   += cv;
        used_llrs[i]  += 1;
      }
    }
  }

  // ----- Finalize -----
  NumericVector scores(n);

  for (int i = 0; i < n; ++i) {
    // Mean leaf log likelihood ratio; SNPs with no leaf evidence are NA
    if (used_llrs[i] > 0)
      scores[i] = llrs_sum[i] / static_cast<double>(used_llrs[i]);
    else
      scores[i] = NA_REAL;
  }

  return List::create(
    _["scores"] = scores
  );
}
