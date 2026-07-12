#' Sensitivity of boosted scores to the background set bandwidth
#'
#' Varies the bandwidth parameter \code{k} used to define the background SNP
#' set (the symmetric band within \code{k} standard units of the training
#' center), holding the extreme set fixed at \code{extreme_k}. For each
#' \code{k} in \code{k_grid}, leaf LLRs are recomputed from the resulting
#' training background set and gain curves are evaluated on the test partition.
#' When \code{R > 1}, background sets larger than the smallest in \code{k_grid}
#' are subsampled \code{R} times to their common minimum size, ensuring that
#' performance differences across bandwidths are attributable to which SNPs are
#' selected rather than how many.
#'
#' @param extreme_k Positive numeric scalar. Fixed bandwidth for the extreme
#'   set, in units of the training spread. Background \code{k} values above
#'   \code{extreme_k} are dropped with a warning. Equality is allowed because
#'   the background interval is strict and the extreme boundary is inclusive.
#' @param k_grid Numeric vector of positive background bandwidth values to
#'   evaluate. Values \eqn{>} \code{extreme_k} are dropped.
#' @inheritParams .boosted_params
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{results}}{A \code{data.table} with one row per
#'     \code{(k, rep, frac_screened)} combination, containing columns
#'     \code{mode} (\code{"background"}), \code{k}, \code{rep},
#'     \code{bg_size_used},
#'     \code{frac_screened}, \code{n_screened}, \code{recall},
#'     \code{lift_curve}, and \code{score_threshold}.}
#'   \item{\code{summary}}{A \code{data.table} summarizing mean and SD of
#'     \code{lift_curve} and \code{recall} across replicates, grouped by
#'     \code{k}, \code{score}, and \code{frac_screened}.}
#' }
#' @export
background_sensitivity <- function(boosted,
                                   extreme_k      = 1.51,
                                   lower_tail     = FALSE,
                                   center         = c("mean", "median"),
                                   scale          = c("sd", "mad"),
                                   k_grid         = seq(0.5, 1.5, 0.25),
                                   R              = 1000L,
                                   alpha          = 0.5,
                                   method         = c("hard_label", "soft_label"),
                                   model_llr_weight = 0.5,
                                   gain_grid      = seq(0.05, 0.50, 0.05),
                                   verbose        = FALSE,
                                   progress_every = NULL) {
  # Signature & basic checks
  FUN <- "background_sensitivity"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  if (!inherits(boosted, "boosted")) {
    stop(
      sprintf(
        "[%s] Input must be an object of class 'boosted' (from make_boosted()).",
        FUN
      )
    )
  }

  # Pull data from boosted
  yvar_train     <- boosted$yvar_train
  n_yvar_test    <- boosted$n_yvar_test
  extr_idx_test  <- boosted$extr_idx_test
  bg_idx_test    <- boosted$bg_idx_test
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  test_leaf_map  <- boosted$test_leaf_map
  oof_prediction_train <- boosted$oof_prediction_train

  # Parse input arguments
  center <- match.arg(center)
  scale  <- match.arg(scale)
  method <- match.arg(method)
  use_soft_labels <- method == "soft_label"
  use_hard_labels <- method == "hard_label"
  if (!is.numeric(model_llr_weight) || length(model_llr_weight) != 1L ||
      !is.finite(model_llr_weight) || model_llr_weight < 0) {
    stop(sprintf("[%s] model_llr_weight must be one non-negative finite number.", FUN))
  }
  if (!is.logical(lower_tail) || length(lower_tail) != 1L ||
      is.na(lower_tail)) {
    stop(sprintf("[%s] lower_tail must be TRUE or FALSE.", FUN))
  }
  if (use_soft_labels && is.null(oof_prediction_train)) {
    stop(sprintf("[%s] method = 'soft_label' requires boosted$oof_prediction_train.", FUN))
  }
  resample_on <- (R > 1L)
  fixed_bg_n  <- NA_integer_

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  # Ensure background k grid does not intersect with fixed extreme_k
  k_ok <- k_grid[k_grid <= extreme_k]
  # If problematic ks are present, try dropping them
  if (length(k_ok) < length(k_grid)) {
    message(
      sprintf(
        "[%s] dropped background ks to avoid overlap with extreme set. new k_grid is: %s",
        FUN,
        paste(format(k_ok, trim = TRUE), collapse = ", ")
      )
    )
  }
  # Fail if grid is empty
  if (!length(k_ok)) {
    stop(sprintf("[%s] all ks overlapped the extreme set definition.", FUN))
  }

  # Pre-allocate containers for per-k results and background set sizes
  per_k_sets <- vector("list", length(k_ok))
  names(per_k_sets) <- as.character(k_ok)
  bg_sizes <- integer(length(k_ok))

  # For each k, define sets and store set indices
  for (i in seq_along(k_ok)) {
    k <- k_ok[i]
    # Define set thresholds using train scale and center
    cls <- set_snp_class(
      yvar_train = yvar_train,
      # extreme band fixed to extreme_k
      extreme_k  = extreme_k,
      lower_tail = lower_tail,
      bg_band_k  = k,
      center     = center,
      scale      = scale
    )

    # Fetch train set indices. Test labels are fixed to boosted$extr_idx_test
    # throughout — the gain curve denominator does not shift with k.
    train_sets <- cls$apply(yvar_train)

    # Store set indices and background set size
    per_k_sets[[i]] <- list(
      # training extreme set (constant)
      E_tr = train_sets$extr_idx,
      # training background for this k
      B_tr = train_sets$bg_idx
    )
    bg_sizes[i] <- length(train_sets$bg_idx)

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%s] background k=%.3f → |extreme_train|=%d, |background_train|=%d",
          FUN,
          k,
          length(train_sets$extr_idx),
          length(train_sets$bg_idx)
        )
      )
    }
  }

  # If resampling is enabled, use the smallest set as the fixed background size
  if (resample_on) {
    min_bg     <- min(bg_sizes)
    fixed_bg_n <- min_bg
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%s] R=%d -> Resampling enabled; fixed background size per iteration: %d (smallest set size)",
          FUN,
          R,
          fixed_bg_n
        )
      )
    }
  }

  # Initialize container for results
  out <- vector("list", length(k_ok) * max(1L, R))
  rr  <- 1L

  # Main loop over background k
  for (i in seq_along(k_ok)) {
    k        <- k_ok[i]
    # training extreme set (constant)
    E_tr     <- per_k_sets[[i]]$E_tr
    # train background set for this k
    B_tr_all <- per_k_sets[[i]]$B_tr

    # Determine if this k is the smallest background set
    is_smallest_k <- (resample_on && length(B_tr_all) == fixed_bg_n)
    # Number of iterations to run for this k:
    R_eff <- if (is_smallest_k)
      1L
    else
      R

    ce_fixed_all <- NULL
    if (use_hard_labels) {
      # Hard-label path: E_tr is fixed for this k, while B_tr may vary by
      # replicate. Count the fixed E side once and reuse it in the hot loop.
      ce_fixed_all <- .stack_leaf_counts(
        idx = E_tr,
        train_leaf_map = train_leaf_map,
        tree_idx = seq_len(Tm),
        use_fast_counts = TRUE
      )
    }

    # Per-k iterations
    for (r in seq_len(R_eff)) {
      # If resample_on and |B_tr_all| > fixed_bg_n: sample fixed_bg_n indices
      if (resample_on && length(B_tr_all) > fixed_bg_n) {
        B_tr <- sort(sample(B_tr_all, fixed_bg_n, replace = FALSE))
      } else {
        if (resample_on && length(B_tr_all) < fixed_bg_n) {
          stop(
            sprintf(
              "[%s] internal: background set for k=%.3f smaller than fixed_bg_n=%d.",
              FUN,
              k,
              fixed_bg_n
            )
          )
        }
        # Otherwise (i.e. R == 1 or this is the smallest set),
        # use the full background set for this k
        B_tr <- B_tr_all
      }

      if (use_hard_labels) {
        # Hard-label path: fixed E counts, varying B counts.
        cb_all <- .stack_leaf_counts(
          idx = B_tr,
          train_leaf_map = train_leaf_map,
          tree_idx = seq_len(Tm),
          use_fast_counts = TRUE
        )
        enr <- .leaf_llrs(
          ce_all = ce_fixed_all,
          cb_all = cb_all,
          n_leaves = train_leaf_map$n_leaves,
          tree_idx = seq_len(Tm),
          N_E = length(E_tr),
          N_B = length(B_tr),
          alpha = alpha
        )
      } else {
        # Soft-label path: q depends on the current E/B sets, so support and
        # fractional leaf masses are recomputed for each replicate.
        support <- .soft_label_support(
          extr_idx       = E_tr,
          bg_idx         = B_tr,
          oof_prediction_train = oof_prediction_train,
          yvar_train     = yvar_train,
          model_llr_weight = model_llr_weight,
          lower_tail     = lower_tail
        )
        ce_all <- .stack_leaf_masses(
          idx = support$labelled_idx,
          weights = support$q,
          train_leaf_map = train_leaf_map,
          tree_idx = seq_len(Tm)
        )
        cb_all <- .stack_leaf_masses(
          idx = support$labelled_idx,
          weights = 1 - support$q,
          train_leaf_map = train_leaf_map,
          tree_idx = seq_len(Tm)
        )
        enr <- .leaf_llrs(
          ce_all = ce_all,
          cb_all = cb_all,
          n_leaves = train_leaf_map$n_leaves,
          tree_idx = seq_len(Tm),
          N_E = support$diagnostics$n_extreme_effective,
          N_B = support$diagnostics$n_background_effective,
          alpha = alpha
        )
        enr$diagnostics <- support$diagnostics
      }

      # Score test SNPs
      s_te <- .score_snps(
        test_leaf_map     = test_leaf_map,
        leaf_llrs_by_tree = enr,
        Tm                = Tm,
        n                 = n_yvar_test
      )

      # Optional progress update
      if (!is.null(progress_every) && progress_every > 0L &&
          (r %% progress_every == 0L)) {
        message(sprintf("[%s] k=%.3f → iteration %d/%d", FUN, k, r, R))
      }

      # Gain curve
      s  <- s_te$scores
      gc <- .gain_curve(
        scores   = s,
        n        = n_yvar_test,
        extr_idx = extr_idx_test,
        bg_idx   = bg_idx_test,
        grid     = grid_eval
      )

      gc[, `:=`(
        mode         = "background",
        k            = k,
        rep          = r,
        bg_size_used = as.integer(length(B_tr))
      )]

      data.table::setcolorder(
        gc,
        c(
          "mode",
          "k",
          "rep",
          "bg_size_used",
          "frac_screened",
          "n_screened",
          "tp",
          "fp_bg",
          "recall",
          "precision_vs_bg",
          "lift_curve",
          "score_threshold"
        )
      )

      # Stash results
      out[[rr]] <- gc
      rr <- rr + 1L
    }
  }

  # Bind all iterations
  results <-
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!nrow(results)) {
    warning(sprintf("[%s] no rows produced; check k_grid and inputs.", FUN))
    return(invisible(NULL))
  }

  # Summarize across iterations per k
  summary <- results[, .(
    lift_curve_mean      = mean(lift_curve),
    lift_curve_sd        = sd(lift_curve),
    recall_mean          = mean(recall),
    recall_sd            = sd(recall),
    precision_vs_bg_mean = mean(precision_vs_bg),
    precision_vs_bg_sd   = sd(precision_vs_bg)
  ), by = .(k, frac_screened)]

  list(results = results[],
       summary = summary[])
}
