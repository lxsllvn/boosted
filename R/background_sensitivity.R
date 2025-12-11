#' Title
#'
#' @param boosted
#' @param extreme_k
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
background_sensitivity <- function(boosted,
                                   extreme_k      = 1.51,
                                   lower_tail     = FALSE,
                                   center         = c("mean", "median"),
                                   scale          = c("sd", "mad"),
                                   k_grid         = seq(0.5, 1.5, 0.25),
                                   R              = 1000L,
                                   alpha          = 0.5,
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
  fixed_bg_n  <- NA_integer_

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  # Ensure background k grid does not intersect with fixed extreme_k
  k_ok <- k_grid[k_grid < extreme_k]
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

    # Fetch train and test set indices
    train_sets <- cls$apply(yvar_train)
    test_sets  <- cls$apply(yvar_test)

    # Store set indices and background set size
    per_k_sets[[i]] <- list(
      E_tr = train_sets$extr_idx,
      B_tr = train_sets$bg_idx,
      E_te = test_sets$extr_idx
    )
    bg_sizes[i] <- length(train_sets$bg_idx)

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%s] background k=%.3f → |extreme_train|=%d, |background_train|=%d, |extreme_test|=%d",
          FUN,
          k,
          length(train_sets$extr_idx),
          length(train_sets$bg_idx),
          length(test_sets$extr_idx)
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
    # test extreme set (constant)
    E_te     <- per_k_sets[[i]]$E_te

    # Determine if this k is the smallest background set
    is_smallest_k <- (resample_on && length(B_tr_all) == fixed_bg_n)

    # Number of iterations to run for this k:
    R_eff <- if (is_smallest_k)
      1L
    else
      R

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
        mode         = "background",
        k            = k,
        rep          = r,
        score        = "contrast",
        bg_size_used = as.integer(length(B_tr))
      )]

      data.table::setcolorder(
        gc,
        c(
          "mode",
          "k",
          "rep",
          "score",
          "bg_size_used",
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
