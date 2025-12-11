#' Title
#'
#' @param extr_idx
#' @param bg_idx
#' @param train_leaf_map
#' @param N_extr
#' @param N_bg
#' @param tree_idx
#' @param alpha
#' @param return_ids
#'
#' @return
#' @keywords internal
#'
#' @examples

.leaf_llrs <- function(extr_idx,
                       bg_idx,
                       train_leaf_map,
                       N_extr,
                       N_bg,
                       tree_idx,
                       alpha = 0.5,
                       return_ids = FALSE) {

  J <- length(tree_idx)
  # dense leaf IDs for each tree and SNP
  invs <- train_leaf_map$dense_leaf_ids # invs[[1:Tm]][1:length(n_train)]
  # number of dense leaves per tree
  Lvec <- train_leaf_map$n_leaves       # vector length(1:Tm)

  # Initialize container for results
  leaf_llrs_by_tree <- vector("list", J)

  # Optionally return native leaf IDs
  native_leaf_ids <-
    if (isTRUE(return_ids))
      vector("list", J)
  else
    NULL

  # Near-empirical LLRs (alpha == 0) with tiny epsilon to avoid log(0)
  if (alpha == 0) {
    eps <- 1e-12
    # For each tree j in 1..J,
    for (j in seq_len(J)) {
      t   <- tree_idx[j] # 1-based tree index
      L   <- Lvec[t]     # number of leaves in tree t
      inv <- invs[[t]]   # length = n_train

      # Count number of extremes/background SNPs per dense leaf ID
      ce <- tabulate(inv[extr_idx], nbins = L)
      cb <- tabulate(inv[bg_idx],   nbins = L)

      llrs <- rep(NA_real_, L)

      # Only compute LLRs for leaves that have at least one labeled SNP
      has_any <- (ce + cb) > 0L
      if (any(has_any)) {
        # Empirical probabilities under each leaf for extreme/background SNPs
        pE_raw <- ce / N_extr
        pB_raw <- cb / N_bg

        # Epsilon floor to avoid zeros
        pE <- pmax(pE_raw, eps)
        pB <- pmax(pB_raw, eps)

        llrs[has_any] <- log(pE[has_any] / pB[has_any])
      }

      leaf_llrs_by_tree[[j]] <- llrs

      if (isTRUE(return_ids))
        native_leaf_ids[[j]] <- train_leaf_map$native_leaf_ids[[t]]
    }
  } else {
    # LLRs with Jeffreys prior
    for (j in seq_len(J)) {
      t   <- tree_idx[j] # 1-based tree index
      L   <- Lvec[t]     # number of leaves in tree t
      inv <- invs[[t]]   # length = n_train

      ce <- tabulate(inv[extr_idx], nbins = L)
      cb <- tabulate(inv[bg_idx],   nbins = L)

      llrs <- rep(NA_real_, L)

      # Compute only for leaves with ≥ 1 labeled SNP
      has_any <- (ce + cb) > 0L
      if (any(has_any)) {
        # Jeffreys-prior smoothing (applied per leaf)
        # pE = (ce + α) / (N_extr + α * L)
        # pB = (cb + α) / (N_bg + α * L)
        # α * L distributes smoothing across leaves
        denom_E <- N_extr + alpha * L
        denom_B <- N_bg   + alpha * L

        pE <- (ce + alpha) / denom_E
        pB <- (cb + alpha) / denom_B

        llrs[has_any] <- log(pE[has_any] / pB[has_any])
      }

      leaf_llrs_by_tree[[j]] <- llrs

      if (isTRUE(return_ids))
        native_leaf_ids[[j]] <- train_leaf_map$native_leaf_ids[[t]]
    }
  }

  # Return results
  if (isTRUE(return_ids)) {
    list(leaf_llrs_by_tree = leaf_llrs_by_tree,
         native_leaf_ids   = native_leaf_ids)
  } else {
    list(leaf_llrs_by_tree = leaf_llrs_by_tree)
  }
}
