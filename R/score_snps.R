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
#' @param test_leaf_map Named list returned by \code{.build_test_leaf_map},
#'   containing \code{dense_leaf_ids}: a list of length \code{Tm} where
#'   element \code{t} is an integer vector of length \eqn{n_\text{test}}{}
#'   giving each test SNP's position in the training leaf vocabulary for
#'   tree \code{t} (0 for leaves unseen in training).
#' @param leaf_llrs_by_tree Named list returned by \code{.leaf_llrs} or
#'   \code{.leaf_llrs_fast}, containing \code{leaf_llrs_by_tree}: a list of
#'   length \code{Tm} where element \code{t} is a numeric vector of LLRs
#'   indexed by dense leaf position.
#' @param Tm Integer scalar: number of trees.
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
  # Initialize accumulators
  llrs_sum  <- numeric(n)
  used_llrs <- integer(n)

  # Extract dense leaf IDs for test SNPs
  pos_list <- test_leaf_map$dense_leaf_ids
  # Extract leaf LLRs
  con_list <- leaf_llrs_by_tree$leaf_llrs_by_tree

  # Loop over trees
  for (t in seq_len(Tm)) {
    pos_t <- pos_list[[t]]
    if (length(pos_t) == 0L)
      next
    # only leaves seen in the labeled training set
    ok  <- pos_t > 0L
    if (!any(ok))
      next

    idx   <- which(ok)
    pidx  <- pos_t[idx]           # dense leaf indices 1..L_t
    v_con <- con_list[[t]][pidx]  # may contain NA

    good <- !is.na(v_con)
    if (any(good)) {
      ii            <- idx[good]
      llrs_sum[ii]  <- llrs_sum[ii]  + v_con[good]
      used_llrs[ii] <- used_llrs[ii] + 1L
    }
  }

  # Finalize: mean over contributing trees; SNPs with no leaf evidence are NA
  scores        <- rep(NA_real_, n)
  has_s         <- used_llrs > 0L
  scores[has_s] <- llrs_sum[has_s] / used_llrs[has_s]

  list(scores = scores)
}

#' Internal dispatcher for score_snps
#' @keywords internal
.score_snps <- function(test_leaf_map,
                        leaf_llrs_by_tree,
                        Tm,
                        n, use_rcpp = TRUE) {
  if (use_rcpp && exists(".score_snps_rcpp", mode = "function", inherits = TRUE)) {
    .score_snps_rcpp(test_leaf_map,
                     leaf_llrs_by_tree,
                     Tm,
                     n)
  } else {
    .score_snps_r(test_leaf_map,
                  leaf_llrs_by_tree,
                  Tm,
                  n)
  }
}
