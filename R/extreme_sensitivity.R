#' Title
#'
#' @param boosted
#' @param background_k
#' @param lower_tail
#' @param center
#' @param scale
#' @param k_grid
#' @param R
#' @param alpha
#' @param gain_grid
#' @param verbose
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples
extreme_sensitivity <- function(boosted,
                                background_k   = 1.0,
                                lower_tail     = FALSE,
                                center         = c("mean", "median"),
                                scale          = c("sd", "mad"),
                                k_grid         = seq(1.25, 1.5, 0.25),
                                R              = 1000L,
                                alpha          = 0.5,
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
  yvar_test      <- boosted$yvar_test
  n_yvar_test    <- boosted$n_yvar_test
  extr_idx_test  <- boosted$extr_idx_test
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  test_leaf_map  <- boosted$test_leaf_map
  tree_idx       <- seq_len(Tm)

  # Parse input arguments
  center <- match.arg(center)
  scale  <- match.arg(scale)
  resample_on <- (R > 1L)
  fixed_ext_n <- NA_integer_

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  # Ensure extreme k grid does not intersect with fixed background_k
  k_ok <- k_grid[k_grid > background_k]
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

    # Fetch train and test set indices
    train_sets <- cls$apply(yvar_train)
    test_sets  <- cls$apply(yvar_test)

    # Store set indices and extreme set size
    per_k_sets[[i]] <- list(
      # training extreme set for this k
      E_tr = train_sets$extr_idx,
      # training background set (constant)
      B_tr = train_sets$bg_idx,
      E_te = test_sets$extr_idx
    )
    extr_sizes[i] <- length(train_sets$extr_idx)

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%s] extreme k=%.3f → |extreme_train|=%d, |background_train|=%d, |extreme_test|=%d",
          FUN,
          k,
          length(train_sets$extr_idx),
          length(train_sets$bg_idx),
          length(test_sets$extr_idx)
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
    # test extreme set (constant)
    E_te     <- per_k_sets[[i]]$E_te

    # Determine if this k is the smallest extreme set
    is_smallest_k <-
      (resample_on && length(E_tr_all) == fixed_ext_n)
    # Number of iterations to run for this k
    R_eff <- if (is_smallest_k)
      1L
    else
      R

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

      # Compute leaf LLRs from training SNPs
      enr <- .leaf_llrs(
        extr_idx       = E_tr,
        bg_idx         = B_tr,
        train_leaf_map = train_leaf_map,
        N_extr         = length(E_tr),
        N_bg           = length(B_tr),
        tree_idx       = tree_idx,
        alpha          = alpha
      )

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
        grid     = grid_eval
      )

      gc[, `:=`(
        mode           = "extreme",
        k              = k,
        rep            = r,
        score          = "contrast",
        extr_size_used = as.integer(length(E_tr))
      )]

      data.table::setcolorder(
        gc,
        c(
          "mode",
          "k",
          "rep",
          "score",
          "extr_size_used",
          "frac_screened",
          "n_screened",
          "recall",
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
    lift_curve_mean = mean(lift_curve),
    lift_curve_sd   = sd(lift_curve),
    recall_mean     = mean(recall),
    recall_sd       = sd(recall)
  ), by = .(k, score, frac_screened)]

  list(results = results[],
       summary = summary[])
}
