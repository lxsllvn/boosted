#' Title
#'
#' @param boosted
#' @param gain_grid
#' @param R
#' @param alpha
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples
permute_null <- function(boosted,
                         gain_grid       = seq(0.05, 0.50, 0.05),
                         R               = 1000L,
                         alpha           = 0.5,
                         method          = c("hard_label", "soft_label"),
                         model_llr_weight = 0.5,
                         progress_every  = NULL) {
  # Signature & basic checks
  FUN <- "permute_null"
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
  extr_idx_test  <- boosted$extr_idx_test
  bg_idx_test    <- boosted$bg_idx_test
  N_extr_train   <- boosted$N_extr_train
  N_bg_train     <- boosted$N_bg_train
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  yvar_train     <- boosted$yvar_train
  n_yvar_train   <- boosted$n_yvar_train
  n_yvar_test    <- boosted$n_yvar_test
  test_leaf_map  <- boosted$test_leaf_map
  oof_prediction_train <- boosted$oof_prediction_train

  # Parse input arguments
  method <- match.arg(method)
  use_soft_labels <- method == "soft_label"
  use_hard_labels <- method == "hard_label"
  if (!is.numeric(model_llr_weight) || length(model_llr_weight) != 1L ||
      !is.finite(model_llr_weight) || model_llr_weight < 0) {
    stop(sprintf("[%s] model_llr_weight must be one non-negative finite number.", FUN))
  }
  if (use_soft_labels && is.null(oof_prediction_train)) {
    stop(sprintf("[%s] method = 'soft_label' requires boosted$oof_prediction_train.", FUN))
  }
  tree_idx <- seq_len(Tm)

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  # Initialize container for results
  out <- vector("list", R)
  rr   <- 1L

  # Hot loop
  for (r in seq_len(R)) {
    # Shuffle indices across the full training partition; preserve set sizes
    perm      <- sample.int(
      n_yvar_train,
      N_extr_train + N_bg_train,
      replace = FALSE
    )
    perm_extr <- perm[seq_len(N_extr_train)]
    perm_bg   <- perm[seq_len(N_bg_train) + N_extr_train]

    # Compute leaf LLRs using the randomized indices
    if (use_hard_labels) {
      ce_all <- .stack_leaf_counts(
        idx = perm_extr,
        train_leaf_map = train_leaf_map,
        tree_idx = tree_idx,
        use_fast_counts = TRUE
      )
      cb_all <- .stack_leaf_counts(
        idx = perm_bg,
        train_leaf_map = train_leaf_map,
        tree_idx = tree_idx,
        use_fast_counts = TRUE
      )
      enr_perm <- .leaf_llrs(
        ce_all = ce_all,
        cb_all = cb_all,
        n_leaves = train_leaf_map$n_leaves,
        tree_idx = tree_idx,
        N_E = N_extr_train,
        N_B = N_bg_train,
        alpha = alpha
      )
    } else {
      support <- .soft_label_support(
        extr_idx = perm_extr,
        bg_idx = perm_bg,
        oof_prediction_train = oof_prediction_train,
        yvar_train = yvar_train,
        model_llr_weight = model_llr_weight,
        lower_tail = NULL
      )
      ce_all <- .stack_leaf_masses(
        idx = support$labelled_idx,
        weights = support$q,
        train_leaf_map = train_leaf_map,
        tree_idx = tree_idx
      )
      cb_all <- .stack_leaf_masses(
        idx = support$labelled_idx,
        weights = 1 - support$q,
        train_leaf_map = train_leaf_map,
        tree_idx = tree_idx
      )
      enr_perm <- .leaf_llrs(
        ce_all = ce_all,
        cb_all = cb_all,
        n_leaves = train_leaf_map$n_leaves,
        tree_idx = tree_idx,
        N_E = support$diagnostics$n_extreme_effective,
        N_B = support$diagnostics$n_background_effective,
        alpha = alpha
      )
      enr_perm$diagnostics <- support$diagnostics
    }

    # Score test SNPs
    scores <- .score_snps(
      test_leaf_map     = test_leaf_map,
      leaf_llrs_by_tree = enr_perm,
      Tm                = Tm,
      n                 = n_yvar_test
    )

    # Gain curve
    s  <- scores$scores
    gc <- .gain_curve(
      scores   = s,
      n        = n_yvar_test,
      extr_idx = extr_idx_test,
      bg_idx   = bg_idx_test,
      grid     = grid_eval
    )

    gc[, `:=`(rep = r)]
    out[[rr]] <- gc
    rr <- rr + 1L

    # Optional progress update
    if (!is.null(progress_every) &&
        is.numeric(progress_every) && progress_every > 0L) {
      if (r %% progress_every == 0L) {
        message(sprintf("[%s] iteration %d/%d", FUN, r, R))
      }
    }
  }

  # Return results
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}
