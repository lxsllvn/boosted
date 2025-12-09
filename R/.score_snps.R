#' Title
#'
#' @param test_leaf_map
#' @param leaf_llrs_by_tree
#' @param Tm
#' @param n
#'
#' @return
#' @keywords internal
#'
#' @examples

.score_snps <- function(test_leaf_map,
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

  # Finalize: mean over contributing trees
  scores        <- numeric(n)
  has_s         <- used_llrs > 0L
  scores[has_s] <- llrs_sum[has_s] / used_llrs[has_s]

  list(scores = scores)
}
