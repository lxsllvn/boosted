#' Score test SNPs by averaging per-tree leaf LLRs
#'
#' For each test SNP, looks up the LLR of the leaf it fell into in each tree
#' (via \code{test_leaf_map}), accumulates the LLRs across trees, and returns
#' the mean over trees that contributed a non-\code{NA} LLR. Test SNPs that
#' land in leaves unseen during training (position 0 in the dense map) or
#' in leaves whose LLR is \code{NA} (no labelled training SNPs) do not
#' contribute to the mean for that SNP. SNPs with no contributing trees
#' receive a score of \code{NA}.
#'
#' This low-level R backend expects an already aligned leaf map: element
#' \code{j} of \code{test_leaf_map$dense_leaf_ids} must correspond to element
#' \code{j} of \code{leaf_llrs_by_tree$leaf_llrs_by_tree}. Use
#' \code{.score_snps()} in normal code; it receives the full test leaf map
#' from \code{make_boosted()} and slices it by the \code{tree_idx} metadata
#' returned by \code{.leaf_llrs()}.
#'
#' @param test_leaf_map Named list containing \code{dense_leaf_ids}: an
#'   already aligned list where element \code{j} is an integer vector of
#'   length \eqn{n_\text{test}}{} giving each test SNP's dense leaf position
#'   for the corresponding scored tree (0 for leaves unseen in training).
#' @param leaf_llrs_by_tree Named list returned by \code{.leaf_llrs()},
#'   containing \code{leaf_llrs_by_tree}: a list where element \code{j} is a
#'   numeric vector of LLRs indexed by dense leaf position for the same scored
#'   tree.
#' @param Tm Integer scalar: number of scored trees.
#' @param n Integer scalar: total number of test SNPs.
#'
#' @return A named list with a single element:
#' \describe{
#'   \item{\code{scores}}{Numeric vector of length \code{n}. Element \code{i}
#'     is the mean LLR of SNP \code{i} across all trees that provided a
#'     non-\code{NA} LLR for that SNP. SNPs with no contributing trees
#'     receive \code{NA}.}
#' }
#' @keywords internal
.score_snps_r <- function(test_leaf_map,
                          leaf_llrs_by_tree,
                          Tm,
                          n) {
  # Extract dense leaf IDs for test SNPs and corresponding leaf LLRs.
  pos_list <- test_leaf_map$dense_leaf_ids
  con_list <- leaf_llrs_by_tree$leaf_llrs_by_tree

  # Initialize accumulators.
  llrs_sum <- numeric(n)
  used_llrs <- integer(n)

  # Loop over scored trees.
  for (t in seq_len(Tm)) {
    pos_t <- pos_list[[t]]
    if (length(pos_t) == 0L)
      next

    # Only leaves seen in the labelled training set contribute.
    ok <- pos_t > 0L
    if (!any(ok))
      next

    idx <- which(ok)
    pidx <- pos_t[idx]
    v_con <- con_list[[t]][pidx]

    good <- !is.na(v_con)
    if (any(good)) {
      ii <- idx[good]
      llrs_sum[ii] <- llrs_sum[ii] + v_con[good]
      used_llrs[ii] <- used_llrs[ii] + 1L
    }
  }

  # Finalize: mean over contributing trees; SNPs with no leaf evidence are NA.
  scores <- rep(NA_real_, n)
  has_s <- used_llrs > 0L
  scores[has_s] <- llrs_sum[has_s] / used_llrs[has_s]
  list(scores = scores)
}


#' Internal dispatcher for SNP scoring
#'
#' Normal package callers pass the full \code{test_leaf_map} from
#' \code{make_boosted()}. The leaf-LLR object records which trees were used in
#' \code{tree_idx}; this dispatcher slices the full map to that tree order and
#' then hands the aligned blocks to the R or Rcpp backend.
#'
#' @keywords internal
.score_snps <- function(test_leaf_map,
                        leaf_llrs_by_tree,
                        Tm,
                        n,
                        use_rcpp = TRUE) {
  # The leaf-LLR object defines the scored tree subset and order.
  tree_idx <- leaf_llrs_by_tree$tree_idx
  aligned_test_leaf_map <- list(
    dense_leaf_ids = test_leaf_map$dense_leaf_ids[tree_idx]
  )
  aligned_Tm <- length(tree_idx)

  if (isTRUE(use_rcpp) &&
      exists(".score_snps_rcpp", mode = "function", inherits = TRUE)) {
    .score_snps_rcpp(
      aligned_test_leaf_map,
      leaf_llrs_by_tree,
      aligned_Tm,
      n
    )
  } else {
    .score_snps_r(
      aligned_test_leaf_map,
      leaf_llrs_by_tree,
      aligned_Tm,
      n
    )
  }
}
