#' Stack hard-label leaf counts across selected trees
#'
#' Counts how many selected training SNPs fall in each dense leaf of each
#' requested tree. The result is one numeric vector formed by concatenating
#' tree-specific count blocks in \code{tree_idx} order:
#' \deqn{(c_{1,1}, \ldots, c_{1,L_1}, c_{2,1}, \ldots, c_{2,L_2}, \ldots).}
#' Each selected SNP contributes one unit of mass to exactly one leaf per tree.
#'
#' This is the hard-label evidence primitive used before \code{\link{.leaf_llrs}}:
#' call it once for the E set and once for the B set, then pass the two stacked
#' vectors to \code{\link{.leaf_llrs}}. The compiled count backend is used when
#' available.
#'
#' @param idx Integer vector of selected training SNP indices, typically one
#'   hard-labelled set such as E or B.
#' @param train_leaf_map Named list returned by
#'   \code{\link{.build_train_leaf_map}}, containing compact 1-based
#'   \code{dense_leaf_ids} and per-tree \code{n_leaves}.
#' @param tree_idx Integer vector of 1-based tree indices to stack.
#' @param use_fast_counts Logical. Use the compiled count backend when
#'   available.
#'
#' @return Numeric vector of length
#'   \code{sum(train_leaf_map$n_leaves[tree_idx])}.
#' @keywords internal
.stack_leaf_counts <- function(
    idx,
    train_leaf_map,
    tree_idx = seq_along(train_leaf_map$dense_leaf_ids),
    use_fast_counts = TRUE) {
  # Downstream sensitivity functions pass validated boosted objects, so this
  # primitive keeps validation light and focuses on the tabulation operation.
  idx <- as.integer(idx)
  tree_idx <- as.integer(tree_idx)
  n_leaves <- as.integer(train_leaf_map$n_leaves)

  if (isTRUE(use_fast_counts) &&
      exists(".stack_leaf_counts_rcpp", mode = "function", inherits = TRUE)) {
    return(.stack_leaf_counts_rcpp(
      dense_leaf_ids = train_leaf_map$dense_leaf_ids,
      n_leaves = n_leaves,
      idx = idx,
      tree_idx = tree_idx
    ))
  }

  # Pure R fallback
  unlist(
    lapply(tree_idx, function(t) {
      # dense_leaf_ids[[t]][idx] are compact 1..L_t IDs, so tabulate() gives
      # this tree's per-leaf hard-label counts directly.
      tabulate(
        train_leaf_map$dense_leaf_ids[[t]][idx],
        nbins = n_leaves[[t]]
      )
    }),
    use.names = FALSE
  )
}


#' Stack soft-label weighted leaf masses across selected trees
#'
#' Weighted analogue of \code{\link{.stack_leaf_counts}}. Instead of adding one
#' count per selected SNP, SNP \code{idx[k]} contributes \code{weights[k]} to
#' its dense leaf in each requested tree. In the soft-label scoring path, these
#' weights come from \code{\link{.soft_label_support}}: \code{q} supplies
#' fractional E mass and \code{1 - q} supplies fractional B mass. The returned
#' vector uses the same stacked tree-block layout as
#' \code{\link{.stack_leaf_counts}}.
#'
#' This is the soft-label evidence primitive used before \code{\link{.leaf_llrs}}.
#' For example, if \code{q} is a vector of fractional E support for labelled
#' SNPs, use \code{weights = q} to form E masses and \code{weights = 1 - q} to
#' form B masses. The resulting fractional masses enter the same LLR formula as
#' hard counts.
#'
#' @param idx Integer vector of selected training SNP indices.
#' @param weights Numeric non-negative mass for each selected SNP, aligned
#'   position-by-position with \code{idx}.
#' @inheritParams .stack_leaf_counts
#' @param use_rcpp Logical. Use the compiled mass backend when available.
#'
#' @return Numeric vector of length
#'   \code{sum(train_leaf_map$n_leaves[tree_idx])}.
#' @keywords internal
.stack_leaf_masses <- function(
    idx,
    weights,
    train_leaf_map,
    tree_idx = seq_along(train_leaf_map$dense_leaf_ids),
    use_rcpp = TRUE) {
  idx <- as.integer(idx)
  weights <- as.numeric(weights)
  tree_idx <- as.integer(tree_idx)
  n_leaves <- as.integer(train_leaf_map$n_leaves)

  if (isTRUE(use_rcpp) &&
      exists(".stack_leaf_masses_rcpp", mode = "function", inherits = TRUE)) {
    return(.stack_leaf_masses_rcpp(
      dense_leaf_ids = train_leaf_map$dense_leaf_ids,
      n_leaves = n_leaves,
      idx = idx,
      weights = weights,
      tree_idx = tree_idx
    ))
  }

  #Pure R fallback
  unlist(
    lapply(tree_idx, function(t) {
      leaf_id <- train_leaf_map$dense_leaf_ids[[t]][idx]
      out <- numeric(n_leaves[[t]])
      # Zero weights are mathematically harmless but can be skipped in the R
      # fallback to avoid unnecessary rowsum() work.
      keep <- weights != 0
      if (any(keep)) {
        sums <- rowsum(weights[keep], group = leaf_id[keep], reorder = FALSE)
        out[as.integer(rownames(sums))] <- as.numeric(sums[, 1L])
      }
      out
    }),
    use.names = FALSE
  )
}


#' Calculate empirical-Bayes-adjusted fractional E support
#'
#' Fits E and B kernel density estimates to oriented, scaled OOF training
#' predictions for the observed E/B reference sets, then converts their density
#' log-ratio into fractional E support. The orientation is chosen so larger
#' margins mean more E-like evidence; when \code{lower_tail = NULL}, that
#' direction is inferred from the observed E and B response means.
#'
#' The returned \code{q} values are posterior-style soft labels:
#' \code{q[i]} is the fractional E mass for labelled SNP \code{labelled_idx[i]},
#' while \code{1 - q[i]} is its fractional B mass. These masses are stacked by
#' \code{\link{.stack_leaf_masses}} and then passed to \code{\link{.leaf_llrs}}.
#'
#' The model density log-ratio is clipped before applying the logistic transform
#' so a small number of KDE tail points cannot dominate the leaf evidence. A
#' zero \code{model_llr_weight} yields neutral support, \code{q = 0.5}, for all
#' labelled SNPs.
#'
#' @param extr_idx Integer vector of 1-based training SNP indices for the
#'   observed E reference set.
#' @param bg_idx Integer vector of 1-based training SNP indices for the observed
#'   B reference set.
#' @param oof_prediction_train Numeric vector of honest out-of-fold training
#'   predictions, one per training SNP.
#' @param yvar_train Numeric training response vector, used to infer orientation
#'   when \code{lower_tail = NULL} and to scale the OOF margin.
#' @param model_llr_weight Non-negative scalar multiplying the KDE log-density
#'   ratio before logistic conversion to \code{q}. Larger values make the soft
#'   support more model-driven; zero gives neutral \code{q = 0.5}.
#' @param lower_tail Optional logical. If \code{TRUE}, smaller response values
#'   are treated as more E-like; if \code{FALSE}, larger values are treated as
#'   more E-like; if \code{NULL}, orientation is inferred from the observed
#'   E/B means.
#'
#' @return Named list containing \code{labelled_idx}, observed labels, oriented
#'   margins, E/B density estimates at labelled SNPs, model log-density ratios,
#'   fractional E support \code{q}, and diagnostic fields including effective
#'   E and B masses.
#'
#' @keywords internal
.soft_label_support <- function(
    extr_idx,
    bg_idx,
    oof_prediction_train,
    yvar_train,
    model_llr_weight = 0.5,
    lower_tail = NULL) {
  extr_idx <- as.integer(extr_idx)
  bg_idx <- as.integer(bg_idx)
  oof_prediction_train <- as.numeric(oof_prediction_train)
  scale_y <- stats::sd(yvar_train, na.rm = TRUE)
  if (is.null(lower_tail)) {
    lower_tail <- mean(yvar_train[extr_idx], na.rm = TRUE) <
      mean(yvar_train[bg_idx], na.rm = TRUE)
  }

  # Orient and scale the OOF prediction margin so the KDE comparison is
  # invariant to tail direction and response scale.
  direction <- if (isTRUE(lower_tail)) -1 else 1
  margin_all <- direction * oof_prediction_train / scale_y
  labelled_idx <- c(extr_idx, bg_idx)
  observed_label <- c(rep("E", length(extr_idx)), rep("B", length(bg_idx)))
  is_observed_e <- observed_label == "E"
  margin_z <- margin_all[labelled_idx]
  e_margin <- margin_all[extr_idx]
  b_margin <- margin_all[bg_idx]
  pooled_margin <- c(e_margin, b_margin)

  margin_range <- range(pooled_margin)
  margin_width <- diff(margin_range)
  if (!is.finite(margin_width) || margin_width <= 0) {
    margin_width <- stats::sd(pooled_margin)
  }
  if (!is.finite(margin_width) || margin_width <= 0) {
    margin_width <- 1
  }
  density_from <- margin_range[[1L]] - 0.10 * margin_width
  density_to <- margin_range[[2L]] + 0.10 * margin_width

  density_bw <- stats::bw.nrd0(pooled_margin)
  if (!is.finite(density_bw) || density_bw <= 0) {
    density_bw <- stats::sd(pooled_margin) / 5
  }
  if (!is.finite(density_bw) || density_bw <= 0) {
    density_bw <- 0.1
  }

  # Use a common bandwidth/range for E and B so their log-density ratio is a
  # direct model-evidence comparison on the labelled SNP margins.
  e_kde <- stats::density(
    e_margin,
    bw = density_bw,
    from = density_from,
    to = density_to,
    n = 2048L,
    adjust = 1
  )
  b_kde <- stats::density(
    b_margin,
    bw = density_bw,
    from = density_from,
    to = density_to,
    n = 2048L,
    adjust = 1
  )
  e_density <- stats::approx(
    e_kde$x,
    e_kde$y,
    xout = margin_z,
    rule = 2
  )$y
  b_density <- stats::approx(
    b_kde$x,
    b_kde$y,
    xout = margin_z,
    rule = 2
  )$y

  density_floor <- 1e-4 / max(density_to - density_from, .Machine$double.eps)
  model_llr <- log(pmax(e_density, density_floor)) -
    log(pmax(b_density, density_floor))
  # Keep the soft labels empirical-Bayes-like: model evidence can adjust the
  # observed labels, but extreme KDE ratios are not allowed to overwhelm them.
  model_llr <- pmax(-6, pmin(6, model_llr))

  q <- stats::plogis(model_llr_weight * model_llr)
  q[!is.finite(q)] <- 0.5
  q <- pmin(1, pmax(0, q))

  n_extreme_effective <- sum(q)
  n_background_effective <- sum(1 - q)

  list(
    labelled_idx = labelled_idx,
    observed_label = observed_label,
    margin = margin_z,
    density_E = e_density,
    density_B = b_density,
    model_llr = model_llr,
    q = q,
    diagnostics = list(
      method = "soft_label",
      model_llr_weight = model_llr_weight,
      n_extreme_train = length(extr_idx),
      n_background_train = length(bg_idx),
      n_labelled_train = length(labelled_idx),
      n_extreme_effective = n_extreme_effective,
      n_background_effective = n_background_effective,
      mean_q_extreme_observed = mean(q[is_observed_e]),
      mean_q_background_observed = mean(q[!is_observed_e]),
      lower_tail = lower_tail,
      scale_y = scale_y,
      density_bw = density_bw,
      density_from = density_from,
      density_to = density_to,
      density_floor = density_floor
    )
  )
}

#' Compute per-leaf log-likelihood ratios from stacked leaf evidence
#'
#' Converts stacked E and B leaf evidence into per-tree leaf log-likelihood
#' ratios (LLRs). For each requested tree and each occupied dense leaf, the LLR
#' compares the estimated probability of landing in that leaf under the E
#' distribution against the corresponding probability under the B distribution:
#' \deqn{\log\left(p_E(\ell) / p_B(\ell)\right).}
#'
#' The inputs \code{ce_all} and \code{cb_all} may be ordinary hard-label counts
#' from \code{\link{.stack_leaf_counts}} or fractional soft-label masses from
#' \code{\link{.stack_leaf_masses}}. Both use the same stacked tree-block layout
#' defined by \code{tree_idx}.
#'
#' With \code{alpha > 0}, Jeffreys/prior smoothing shrinks each tree's leaf
#' probabilities toward a uniform distribution over that tree's \eqn{L} dense
#' leaves:
#' \deqn{p_E(\ell) = (c_E(\ell) + \alpha) / (N_E + \alpha L),}
#' \deqn{p_B(\ell) = (c_B(\ell) + \alpha) / (N_B + \alpha L).}
#' With \code{alpha = 0}, empirical probabilities are used with a tiny epsilon
#' floor to avoid \code{log(0)}. Leaves with no E or B evidence
#' (\code{ce + cb == 0}) are left as \code{NA} rather than being assigned a
#' spurious LLR; downstream scoring ignores those missing leaf contributions.
#'
#' @param ce_all Numeric vector of stacked extreme counts or masses.
#' @param cb_all Numeric vector of stacked background counts or masses.
#' @param n_leaves Integer vector giving the leaf count for every tree.
#' @param tree_idx Integer vector defining the tree order in the stacked
#'   vectors.
#' @param N_E Extreme count or effective mass.
#' @param N_B Background count or effective mass.
#' @param alpha Non-negative finite smoothing concentration.
#' @param use_rcpp Logical. Use the compiled LLR backend when available.
#'
#' @return Named list with \code{leaf_llrs_by_tree} and \code{tree_idx}.
#'   \code{leaf_llrs_by_tree[[j]]} is a numeric vector of length
#'   \code{n_leaves[tree_idx[j]]}; \code{NA} entries mark leaves with no
#'   labelled evidence.
#' @keywords internal
.leaf_llrs <- function(
    ce_all,
    cb_all,
    n_leaves,
    tree_idx,
    N_E,
    N_B,
    alpha = 0.5,
    use_rcpp = TRUE) {
  tree_idx <- as.integer(tree_idx)
  n_leaves <- as.integer(n_leaves)
  ce_all <- as.numeric(ce_all)
  cb_all <- as.numeric(cb_all)

  if (isTRUE(use_rcpp) &&
      exists(".leaf_llrs_rcpp", mode = "function", inherits = TRUE)) {
    return(.leaf_llrs_rcpp(
      ce_all = ce_all,
      cb_all = cb_all,
      n_leaves = n_leaves,
      tree_idx = tree_idx,
      N_E = as.numeric(N_E),
      N_B = as.numeric(N_B),
      alpha = as.numeric(alpha)
    ))
  }

  # Pure R implementation
  offsets <- c(0L, cumsum(n_leaves[tree_idx]))
  eps <- 1e-12
  leaf_llrs_by_tree <- vector("list", length(tree_idx))

  for (j in seq_along(tree_idx)) {
    L <- n_leaves[[tree_idx[[j]]]]
    # Slice this tree's block out of the stacked count/mass vectors.
    at <- seq.int(offsets[[j]] + 1L, offsets[[j + 1L]])
    ce <- ce_all[at]
    cb <- cb_all[at]
    has_any <- (ce + cb) > 0
    llrs <- rep(NA_real_, L)

    if (any(has_any)) {
      if (alpha == 0) {
        # Near-empirical LLRs. The epsilon floor avoids infinite LLRs when one
        # class has zero observed evidence in an occupied leaf.
        p_e <- pmax(ce / N_E, eps)
        p_b <- pmax(cb / N_B, eps)
      } else {
        # Per-tree smoothing: alpha is applied to each dense leaf, so alpha * L
        # is the total prior mass added to each class distribution in this tree.
        p_e <- (ce + alpha) / (N_E + alpha * L)
        p_b <- (cb + alpha) / (N_B + alpha * L)
      }
      llrs[has_any] <- log(p_e[has_any] / p_b[has_any])
    }

    leaf_llrs_by_tree[[j]] <- llrs
  }

  list(
    leaf_llrs_by_tree = leaf_llrs_by_tree,
    tree_idx = tree_idx
  )
}
