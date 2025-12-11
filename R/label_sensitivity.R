#' Title
#'
#' @param boosted
#' @param noise_grid
#' @param R
#' @param alpha
#' @param gain_grid
#' @param target
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples
label_sensitivity <- function(boosted,
                              noise_grid    = seq(0.0, 0.50, 0.1),
                              R             = 1000L,
                              alpha         = 0.5,
                              gain_grid     = seq(0.05, 0.50, 0.05),
                              target        = c("extreme", "background"),
                              progress_every = NULL) {
  # Signature & basic checks
  FUN <- "label_noise_sensitivity"
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
  extr_idx_train <- boosted$extr_idx_train
  bg_idx_train   <- boosted$bg_idx_train
  extr_idx_test  <- boosted$extr_idx_test
  N_extr_train   <- boosted$N_extr_train
  N_bg_train     <- boosted$N_bg_train
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  test_leaf_map  <- boosted$test_leaf_map
  tree_idx       <- seq_len(Tm)
  n_yvar_test    <- boosted$n_yvar_test

  # Parse input arguments
  target    <- match.arg(target)

  # Guard against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  # Initialize container for results
  out <- vector("list", length(noise_grid) * R * length(grid_eval))
  rr <- 1L

  # Precompute the number of training labels to flip to obtain the qith
  # noise level in noise_grid
  for (qi in seq_along(noise_grid)) {
    q <- noise_grid[qi]

    # Reset counters
    n_flip_E <- 0L
    n_flip_B <- 0L

    if (target == "extreme") {
      # Add n_flip_B background impostors to extreme set
      n_flip_B <- floor(q * N_extr_train)
    } else {
      # Add n_flip_E extreme imposters to background set
      n_flip_E <- floor(q * N_bg_train)
    }

    # In each iteration, apply the class flip to the training SNPs
    for (r in seq_len(R)) {
      if (target == "extreme") {
        # Randomly select n_flip_B background labels to add to extreme set
        flip_B <-
          if (n_flip_B > 0L)
            sort(sample(bg_idx_train, n_flip_B, replace = FALSE))
        else
          integer(0L)
        # Add background contamination to extreme set
        E_cor  <- sort(c(extr_idx_train, flip_B))
        # Remove these labels from background set
        B_cor  <- setdiff(bg_idx_train, flip_B)
      } else {
        # Randomly select n_flip_E extreme labels to add to background set
        flip_E <-
          if (n_flip_E > 0L)
            sort(sample(extr_idx_train, n_flip_E, replace = FALSE))
        else
          integer(0L)
        # Add extreme contamination to background set
        B_cor  <- sort(c(bg_idx_train, flip_E))
        # Remove these labels from extreme set
        E_cor  <- setdiff(extr_idx_train, flip_E)
      }

      # Compute the leaf enrichment values using the corrupted index sets
      enr <- .leaf_llrs(
        extr_idx       = E_cor,
        bg_idx         = B_cor,
        N_extr         = length(E_cor),
        N_bg           = length(B_cor),
        alpha          = alpha,
        train_leaf_map = train_leaf_map,
        tree_idx       = tree_idx
      )

      # Score the test SNPs based on the noisy leaf enrichment values
      scores <- .score_snps(
        test_leaf_map     = test_leaf_map,
        leaf_llrs_by_tree = enr,
        Tm                = Tm,
        n                 = n_yvar_test
      )

      # Gain curve
      s  <- scores$scores
      gc <- .gain_curve(
        scores   = s,
        n        = n_yvar_test,
        extr_idx = extr_idx_test,
        grid     = grid_eval
      )

      gc[, `:=`(
        mode  = target,
        q     = q,
        rep   = r,
        score = "contrast"
      )]
      data.table::setcolorder(
        gc,
        c(
          "mode",
          "q",
          "rep",
          "score",
          "frac_screened",
          "n_screened",
          "recall",
          "lift_curve",
          "score_threshold"
        )
      )
      out[[rr]] <- gc
      rr <- rr + 1L

      # Optional progress update
      if (!is.null(progress_every) &&
          is.numeric(progress_every) && progress_every > 0L) {
        if (r %% progress_every == 0L) {
          message(sprintf(
            "[%s] q=%.3f (%d/%d), iteration %d/%d",
            FUN,
            q,
            qi,
            length(noise_grid),
            r,
            R
          ))
        }
      }
    }
  }

  # Collate results
  results <-
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)

  # Summarize across iterations
  summary <- results[, .(
    lift_curve_mean = mean(lift_curve),
    lift_curve_sd   = sd(lift_curve),
    recall_mean     = mean(recall),
    recall_sd       = sd(recall)
  ), by = .(mode, score, q, frac_screened)]

  list(results = results[], summary = summary[])
}
