#' Sensitivity of boosted scores to training label noise
#'
#' Injects controlled label noise into the training extreme or background set
#' by randomly adding a fraction \code{q} of impostor labels (background SNPs
#' added to the extreme set, or vice versa), then recomputing leaf LLRs and
#' gain curves across \code{R} resampling iterations per noise level. This
#' quantifies how much the gain curve degrades as a function of label
#' contamination and helps assess robustness to imperfect class definitions.
#'
#' @param noise_grid Numeric vector of contamination fractions in
#'   \code{[0, 1]}. For \code{target = "extreme"}, a fraction \code{q} of the
#'   extreme set size is drawn from the background set and added to the extreme
#'   labels; for \code{target = "background"}, a fraction \code{q} of the
#'   background set size is drawn from the extreme set.
#' @param target Character string, one of \code{"extreme"} or
#'   \code{"background"}. Which training label set receives the impostor
#'   contamination.
#' @inheritParams .boosted_params
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{results}}{A \code{data.table} with one row per
#'     \code{(q, rep, frac_screened)} combination, containing columns
#'     \code{mode} (the value of \code{target}), \code{q}, \code{rep},
#'     \code{score} (\code{"contrast"}), \code{frac_screened},
#'     \code{n_screened}, \code{recall}, \code{lift_curve}, and
#'     \code{score_threshold}.}
#'   \item{\code{summary}}{A \code{data.table} summarising mean and SD of
#'     \code{lift_curve} and \code{recall} across replicates, grouped by
#'     \code{mode}, \code{score}, \code{q}, and \code{frac_screened}.}
#' }
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
