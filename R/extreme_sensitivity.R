#' Sensitivity of boosted scores to the extreme set boundary
#'
#' Varies the boundary parameter \code{k} used to define the extreme SNP set
#' (the tail beyond \code{k} standard units from the training center), holding
#' the background set fixed at \code{background_k}. For each \code{k} in
#' \code{k_grid}, leaf LLRs are recomputed from the resulting training extreme
#' set and gain curves are evaluated on the test partition. When \code{R > 1},
#' extreme sets larger than the smallest in \code{k_grid} are subsampled
#' \code{R} times to their common minimum size, ensuring that performance
#' differences across boundaries are attributable to which SNPs are labelled
#' extreme rather than how many.
#'
#' @param background_k Positive numeric scalar. Fixed bandwidth for the
#'   background set, in units of the training spread. Extreme \code{k} values
#'   below \code{background_k} are dropped with a warning. Equality is allowed:
#'   the background interval is strict while the extreme boundary is inclusive,
#'   so the sets remain disjoint.
#' @param k_grid Numeric vector of positive extreme boundary values to
#'   evaluate. Values \eqn{<} \code{background_k} are dropped.
#' @inheritParams .boosted_params
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{results}}{A \code{data.table} with one row per
#'     \code{(k, rep, frac_screened)} combination, containing columns
#'     \code{mode} (\code{"extreme"}), \code{k}, \code{rep},
#'     \code{extr_size_used},
#'     \code{frac_screened}, \code{n_screened}, \code{recall},
#'     \code{lift_curve}, and \code{score_threshold}.}
#'   \item{\code{summary}}{A \code{data.table} summarizing mean and SD of
#'     \code{lift_curve} and \code{recall} across replicates, grouped by
#'     \code{k}, and \code{frac_screened}.}
#' }
#' @export
extreme_sensitivity <- function(boosted,
                                background_k   = 1.0,
                                lower_tail     = FALSE,
                                center         = c("mean", "median"),
                                scale          = c("sd", "mad"),
                                k_grid         = seq(1.25, 1.5, 0.25),
                                R              = 1000L,
                                alpha          = 0.5,
                                method         = c("hard_label", "soft_label"),
                                model_llr_weight = 0.5,
                                gain_grid      = seq(0.05, 0.50, 0.05),
                                verbose        = FALSE,
                                progress_every = NULL) {
  # Signature & basic checks
  FUN <- "extreme_sensitivity"
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
  fixed_ext_n <- NA_integer_

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  # Ensure extreme k grid does not intersect with fixed background_k
  k_ok <- k_grid[k_grid >= background_k]
  # If problematic ks are present, try dropping them
  if (length(k_ok) < length(k_grid)) {
    message(
      sprintf(
        "[%s] dropped extreme ks to avoid overlap with background set (background_k=%.3f). new k_grid is: %s",
        FUN,
        background_k,
        paste(format(k_ok, trim = TRUE), collapse = ", ")
      )
    )
  }
  # Fail if grid is empty
  if (!length(k_ok)) {
    stop(
      sprintf(
        "[%s] all extreme ks overlapped the background definition (background_k=%.3f).",
        FUN,
        background_k
      )
    )
  }

  # Pre-allocate containers for per-k results and extreme set sizes
  per_k_sets  <- vector("list", length(k_ok))
  names(per_k_sets) <- as.character(k_ok)
  extr_sizes  <- integer(length(k_ok))

  # For each k, define sets and store set indices
  for (i in seq_along(k_ok)) {
    k <- k_ok[i]
    # Define set thresholds using train scale and center
    cls <- set_snp_class(
      yvar_train = yvar_train,
      extreme_k  = k,
      lower_tail = lower_tail,
      # background band fixed to background_k
      bg_band_k  = background_k,
      center     = center,
      scale      = scale
    )

    # Fetch train set indices. Test labels are fixed to boosted$extr_idx_test
    # throughout — the gain curve denominator does not shift with k.
    train_sets <- cls$apply(yvar_train)

    # Store set indices and extreme set size
    per_k_sets[[i]] <- list(
      # training extreme set for this k
      E_tr = train_sets$extr_idx,
      # training background set (constant)
      B_tr = train_sets$bg_idx
    )
    extr_sizes[i] <- length(train_sets$extr_idx)

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%s] extreme k=%.3f → |extreme_train|=%d, |background_train|=%d",
          FUN,
          k,
          length(train_sets$extr_idx),
          length(train_sets$bg_idx)
        )
      )
    }
  }

  # If resampling is enabled, use the smallest extreme set as the fixed extreme size
  if (resample_on) {
    min_ext     <- min(extr_sizes)
    fixed_ext_n <- min_ext
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%s] R=%d -> Resampling enabled; fixed extreme size per iteration: %d (smallest set size)",
          FUN,
          R,
          fixed_ext_n
        )
      )
    }
  }

  # Initialize container for results
  out <- vector("list", length(k_ok) * max(1L, R))
  rr  <- 1L

  # Main loop over extreme k
  for (i in seq_along(k_ok)) {
    k        <- k_ok[i]
    # training extreme set for this k
    E_tr_all <- per_k_sets[[i]]$E_tr
    # training background set (constant)
    B_tr     <- per_k_sets[[i]]$B_tr

    # Determine if this k is the smallest extreme set
    is_smallest_k <-
      (resample_on && length(E_tr_all) == fixed_ext_n)
    # Number of iterations to run for this k
    R_eff <- if (is_smallest_k)
      1L
    else
      R

    cb_fixed_all <- NULL
    if (use_hard_labels) {
      # Hard-label path: B_tr is fixed for this k, while E_tr may vary by
      # replicate. Count the fixed B side once and reuse it in the hot loop.
      cb_fixed_all <- .stack_leaf_counts(
        idx = B_tr,
        train_leaf_map = train_leaf_map,
        tree_idx = seq_len(Tm),
        use_fast_counts = TRUE
      )
    }

    # Per-k iterations
    for (r in seq_len(R_eff)) {
      # If resample_on and |E_tr_all| > fixed_ext_n: sample fixed_ext_n indices
      if (resample_on && length(E_tr_all) > fixed_ext_n) {
        E_tr <- sort(sample(E_tr_all, fixed_ext_n, replace = FALSE))
      } else {
        # Otherwise (i.e., R == 1 or this is the smallest set),
        # use the full extreme set for this k
        E_tr <- E_tr_all
      }

      if (use_hard_labels) {
        # Hard-label path: varying E counts, fixed B counts.
        ce_all <- .stack_leaf_counts(
          idx = E_tr,
          train_leaf_map = train_leaf_map,
          tree_idx = seq_len(Tm),
          use_fast_counts = TRUE
        )
        enr <- .leaf_llrs(
          ce_all = ce_all,
          cb_all = cb_fixed_all,
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
        mode           = "extreme",
        k              = k,
        rep            = r,
        extr_size_used = as.integer(length(E_tr))
      )]

      data.table::setcolorder(
        gc,
        c(
          "mode",
          "k",
          "rep",
          "extr_size_used",
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
