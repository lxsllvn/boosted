#' Compute per-leaf log-likelihood ratios for a subset of trees
#'
#' For each tree in \code{tree_idx}, tabulates how many extreme and background
#' training SNPs fall in each leaf, then computes a log-likelihood ratio (LLR)
#' comparing the empirical probability of that leaf assignment under the extreme
#' distribution to the background distribution. Jeffreys-prior smoothing
#' (\code{alpha > 0}) shrinks the raw counts toward a uniform prior over
#' leaves, preventing extreme LLRs in sparsely occupied leaves. Leaves that
#' contain no labelled SNPs at all receive \code{NA} rather than a spurious
#' LLR.
#'
#' @param extr_idx Integer vector of 1-based SNP indices identifying the
#'   extreme training set.
#' @param bg_idx Integer vector of 1-based SNP indices identifying the
#'   background training set.
#' @param train_leaf_map Named list returned by \code{.build_train_leaf_map},
#'   containing \code{dense_leaf_ids} (per-tree compact leaf assignments for
#'   all training SNPs) and \code{n_leaves} (leaf count per tree).
#' @param N_extr Integer scalar: total number of extreme training SNPs
#'   (denominator for the extreme probability).
#' @param N_bg Integer scalar: total number of background training SNPs
#'   (denominator for the background probability).
#' @param tree_idx Integer vector of 1-based tree indices to evaluate. Allows
#'   the caller to score a subset of trees without reordering data.
#' @param alpha Numeric scalar \eqn{\ge 0}. Jeffreys-prior concentration; see
#'   \code{\link{.boosted_params}}.
#' @param return_ids Logical. If \code{TRUE}, the native `xgboost` leaf ID for
#'   each dense leaf index is included in the return value. Defaults to
#'   \code{FALSE}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{leaf_llrs_by_tree}}{A list of length \code{length(tree_idx)}.
#'     Element \code{j} is a numeric vector of length \eqn{L_t}{} (the number
#'     of dense leaves in tree \code{tree_idx[j]}) containing the LLR for each
#'     leaf. Leaves with no labelled SNPs carry \code{NA}.}
#'   \item{\code{native_leaf_ids}}{Only present when \code{return_ids = TRUE}.
#'     A list of length \code{length(tree_idx)} where element \code{j} is an
#'     integer vector of the native `xgboost` leaf IDs corresponding to each
#'     position in \code{leaf_llrs_by_tree[[j]]}.}
#' }
#' @keywords internal

.leaf_llrs_from_counts <- function(ce,
                                   cb,
                                   N_extr,
                                   N_bg,
                                   alpha = 0.5) {
  L <- length(ce)
  llrs <- rep(NA_real_, L)
  has_any <- (ce + cb) > 0L

  if (!any(has_any)) {
    return(llrs)
  }

  if (alpha == 0) {
    eps <- 1e-12
    pE <- pmax(ce / N_extr, eps)
    pB <- pmax(cb / N_bg, eps)
  } else {
    denom_E <- N_extr + alpha * L
    denom_B <- N_bg   + alpha * L
    pE <- (ce + alpha) / denom_E
    pB <- (cb + alpha) / denom_B
  }

  llrs[has_any] <- log(pE[has_any] / pB[has_any])
  llrs
}

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
#' @export
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
  invs <- train_leaf_map$dense_leaf_ids
  Lvec <- train_leaf_map$n_leaves

  leaf_llrs_by_tree <- vector("list", J)
  native_leaf_ids <- if (isTRUE(return_ids))
    vector("list", J)
  else
    NULL

  for (j in seq_len(J)) {
    t <- tree_idx[j]
    L <- Lvec[t]
    inv <- invs[[t]]

    ce <- tabulate(inv[extr_idx], nbins = L)
    cb <- tabulate(inv[bg_idx],   nbins = L)

    leaf_llrs_by_tree[[j]] <- .leaf_llrs_from_counts(
      ce     = ce,
      cb     = cb,
      N_extr = N_extr,
      N_bg   = N_bg,
      alpha  = alpha
    )

    if (isTRUE(return_ids)) {
      native_leaf_ids[[j]] <- train_leaf_map$native_leaf_ids[[t]]
    }
  }

  if (isTRUE(return_ids)) {
    list(
      leaf_llrs_by_tree = leaf_llrs_by_tree,
      native_leaf_ids   = native_leaf_ids
    )
  } else {
    list(leaf_llrs_by_tree = leaf_llrs_by_tree)
  }
}


#' Selected-index count + LLR backend for resampling paths
#'
#' A thin wrapper around the compiled resampling backend. The caller supplies
#' the training leaf assignments once via \code{train_leaf_map}, then chooses
#' one of four count modes:
#' \itemize{
#'   \item explicit \code{extr_idx + bg_idx}
#'   \item explicit \code{bg_idx} plus fixed extreme counts
#'   \item explicit \code{extr_idx} plus fixed background counts
#'   \item a fixed labelled pool plus exactly one varying side
#' }
#' The mathematical definitions match \code{.leaf_llrs()} exactly; only the
#' counting strategy changes.
#'
#' @param extr_idx Optional integer vector of 1-based training SNP indices for
#'   the extreme set.
#' @param bg_idx Optional integer vector of 1-based training SNP indices for
#'   the background set.
#' @param train_leaf_map Named list returned by \code{.build_train_leaf_map}.
#' @param N_extr Integer scalar: size of the extreme set.
#' @param N_bg Integer scalar: size of the background set.
#' @param tree_idx Integer vector of 1-based tree indices to evaluate.
#' @param alpha Numeric scalar \eqn{\ge 0}. Jeffreys-prior concentration.
#' @param fixed_cb_all Optional numeric vector of stacked background counts in
#'   \code{tree_idx} order.
#' @param fixed_ce_all Optional numeric vector of stacked extreme counts in
#'   \code{tree_idx} order.
#' @param pool_counts_all Optional numeric vector of stacked counts for a fixed
#'   relabelling pool in \code{tree_idx} order.
#' @param return_counts Logical. If \code{TRUE}, returns the stacked count
#'   vector for exactly one supplied index set instead of LLRs.
#'
#' @return Either a numeric vector of stacked counts (\code{return_counts =
#'   TRUE}) or a named list with element \code{leaf_llrs_by_tree}.
#' @keywords internal
.leaf_llrs_fast <- function(extr_idx = NULL,
                            bg_idx = NULL,
                            train_leaf_map,
                            N_extr,
                            N_bg,
                            tree_idx = seq_along(train_leaf_map$dense_leaf_ids),
                            alpha = 0.5,
                            fixed_cb_all = NULL,
                            fixed_ce_all = NULL,
                            pool_counts_all = NULL,
                            return_counts = FALSE) {

  dense_leaf_ids <- train_leaf_map$dense_leaf_ids
  n_leaves <- as.integer(train_leaf_map$n_leaves)
  tree_idx <- as.integer(tree_idx)

  has_extr <- !is.null(extr_idx)
  has_bg <- !is.null(bg_idx)
  has_fixed_ce <- !is.null(fixed_ce_all)
  has_fixed_cb <- !is.null(fixed_cb_all)
  has_pool <- !is.null(pool_counts_all)

  expected_counts <- sum(n_leaves[tree_idx])

  .leaf_llrs_backend_rcpp(
    dense_leaf_ids     = dense_leaf_ids,
    n_leaves           = n_leaves,
    N_extr             = as.integer(N_extr),
    N_bg               = as.integer(N_bg),
    alpha              = as.numeric(alpha),
    extr_idx_sexp      = if (has_extr) as.integer(extr_idx) else NULL,
    bg_idx_sexp        = if (has_bg) as.integer(bg_idx) else NULL,
    fixed_ce_all_sexp  = if (has_fixed_ce) as.numeric(fixed_ce_all) else NULL,
    fixed_cb_all_sexp  = if (has_fixed_cb) as.numeric(fixed_cb_all) else NULL,
    pool_counts_all_sexp = if (has_pool) as.numeric(pool_counts_all) else NULL,
    tree_idx_sexp      = tree_idx,
    return_counts      = return_counts
  )
}
