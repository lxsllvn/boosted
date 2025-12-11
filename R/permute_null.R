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
  extr_idx_train <- boosted$extr_idx_train
  bg_idx_train   <- boosted$bg_idx_train
  extr_idx_test  <- boosted$extr_idx_test
  N_index_train  <- boosted$N_index_train
  N_extr_train   <- boosted$N_extr_train
  N_bg_train     <- boosted$N_bg_train
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  n_yvar_test    <- boosted$n_yvar_test
  test_leaf_map  <- boosted$test_leaf_map
  tree_idx       <- seq_len(Tm)

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
    # Shuffle indices; preserve set sizes
    perm      <- sample.int(N_index_train)
    perm_extr <- perm[seq_len(N_extr_train)]
    perm_bg   <- perm[seq_len(N_bg_train) + N_extr_train]

    # Compute leaf LLRs using the randomized indices
    enr_perm <- .leaf_llrs(
      extr_idx       = perm_extr,
      bg_idx         = perm_bg,
      N_extr         = N_extr_train,
      N_bg           = N_bg_train,
      alpha          = alpha,
      train_leaf_map = train_leaf_map,
      tree_idx       = tree_idx
    )

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
      grid     = grid_eval
    )

    gc[, `:=`(rep = r, score = "contrast")]
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

