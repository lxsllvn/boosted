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


#' Sparse-matrix accelerated leaf LLR computation for permutation loops
#'
#' A drop-in replacement for \code{.leaf_llrs} designed for use inside
#' permutation hot loops. Instead of calling \code{tabulate()} once per tree,
#' a single sparse matrix-vector product (\code{A \%*\% y}) replaces all
#' \eqn{T_m}{} tabulation calls simultaneously, where \code{A} is a
#' pre-built leaf-by-SNP incidence matrix. This is substantially faster when
#' \code{.leaf_llrs} is called thousands of times with different permuted index
#' vectors but the same leaf structure. The pre-computation cost is paid once
#' by \code{.build_leaf_matrix}.
#'
#' @param perm_extr Integer vector of 1-based SNP indices for the permuted
#'   extreme set (drawn from the full labelled index pool).
#' @param leaf_mat Named list returned by \code{.build_leaf_matrix}, containing
#'   the sparse incidence matrix \code{A}, per-tree \code{offsets},
#'   \code{Lvec} (leaves per tree), and \code{labeled_row_sums} (labelled SNP
#'   count per leaf row).
#' @param N_extr Integer scalar: size of the permuted extreme set.
#' @param N_bg Integer scalar: size of the permuted background set
#'   (\code{N_index_train - N_extr}).
#' @param Tm Integer scalar: total number of trees.
#' @param alpha Numeric scalar \eqn{\ge 0}. Jeffreys-prior concentration; see
#'   \code{\link{.boosted_params}}.
#'
#' @return A named list with a single element \code{leaf_llrs_by_tree}: a list
#'   of length \code{Tm} where element \code{t} is a numeric vector of LLRs
#'   for the dense leaves of tree \code{t}. The format is identical to the
#'   \code{leaf_llrs_by_tree} element returned by \code{.leaf_llrs}.
#' @keywords internal
#'
#' @examples
.leaf_llrs_fast <- function(perm_extr,
                            leaf_mat,
                            N_extr,
                            N_bg,
                            Tm,
                            alpha = 0.5) {
  A        <- leaf_mat$A
  offsets  <- leaf_mat$offsets
  Lvec     <- leaf_mat$Lvec

  # Binary label vector: 1 = extreme, 0 = background
  y <- integer(ncol(A))
  y[perm_extr] <- 1L

  # Sparse matvec replaces all Tm tabulate() calls
  ce_all <- as.vector(A %*% y)
  cb_all <- leaf_mat$labeled_row_sums - ce_all

  # Unpack by tree and compute LLRs
  leaf_llrs_by_tree <- vector("list", Tm)

  for (t in seq_len(Tm)) {
    rows <- seq(offsets[t] + 1L, offsets[t + 1])
    ce   <- ce_all[rows]
    cb   <- cb_all[rows]
    L    <- Lvec[t]

    has_any <- (ce + cb) > 0L
    llrs    <- rep(NA_real_, L)

    if (any(has_any)) {
      denom_E <- N_extr + alpha * L
      denom_B <- N_bg   + alpha * L
      pE <- (ce + alpha) / denom_E
      pB <- (cb + alpha) / denom_B
      llrs[has_any] <- log(pE[has_any] / pB[has_any])
    }
    leaf_llrs_by_tree[[t]] <- llrs
  }

  list(leaf_llrs_by_tree = leaf_llrs_by_tree)
}
