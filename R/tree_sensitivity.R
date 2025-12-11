#' Title
#'
#' @param boosted
#' @param subsample_grid
#' @param R
#' @param alpha
#' @param gain_grid
#' @param topk_frac
#' @param verbose
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples
tree_sensitivity <- function(boosted,
                             subsample_grid = seq(0.10, 0.50, 0.10),
                             R              = 50L,
                             alpha          = 0.5,
                             gain_grid      = seq(0.05, 0.50, 0.05),
                             topk_frac      = 0.05,
                             verbose        = FALSE,
                             progress_every = NULL) {
  # Signature & basic checks
  FUN <- "tree_sensitivity"
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
  N_extr_train   <- boosted$N_extr_train
  N_bg_train     <- boosted$N_bg_train
  extr_idx_test  <- boosted$extr_idx_test
  n_yvar_test    <- boosted$n_yvar_test
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  test_leaf_map  <- boosted$test_leaf_map

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  k_top <- max(1L, as.integer(round(topk_frac * n_yvar_test)))

  # Helper: compute enrichment and score using a subset of trees
  .fit_subset <- function(trees) {
    trees <- sort(unique(as.integer(trees)))
    if (!length(trees))
      stop(sprintf("[%s] internal: empty tree subset.", FUN))

    enr <- .leaf_llrs(
      extr_idx       = extr_idx_train,
      bg_idx         = bg_idx_train,
      train_leaf_map = train_leaf_map,
      N_extr         = N_extr_train,
      N_bg           = N_bg_train,
      tree_idx       = trees,
      alpha          = alpha
    )
    sc <- .score_snps(
      test_leaf_map     = test_leaf_map,
      leaf_llrs_by_tree = enr,
      Tm                = length(trees),
      n                 = n_yvar_test
    )$scores
    list(scores = sc, trees = trees)
  }

  # Baseline: full-model scores for diagnostics
  enr_full <- .leaf_llrs(
    extr_idx       = extr_idx_train,
    bg_idx         = bg_idx_train,
    train_leaf_map = train_leaf_map,
    N_extr         = N_extr_train,
    N_bg           = N_bg_train,
    tree_idx       = seq_len(Tm),
    alpha          = alpha
  )
  s_full <- .score_snps(
    test_leaf_map     = test_leaf_map,
    leaf_llrs_by_tree = enr_full,
    Tm                = Tm,
    n                 = n_yvar_test
  )$scores

  top_full <- utils::head(order(s_full, decreasing = TRUE), k_top)

  # Subsampling runs
  n_runs   <- length(subsample_grid) * R
  out      <- vector("list", n_runs)   # one block per run
  rr       <- 1L

  # Initialize container for per-iteration diagnostics
  diag_rows <- vector("list", n_runs)
  dr <- 1L

  # Initialize container for within-q Jaccard across iterations
  top_sets <- vector("list", length(subsample_grid))
  for (qi in seq_along(subsample_grid)) {
    top_sets[[qi]] <- vector("list", R)
  }

  # Main loop over tree subsets
  for (qi in seq_along(subsample_grid)) {
    q <- subsample_grid[qi]
    n_trees <- as.integer(min(Tm, max(1L, round(q * Tm))))
    if (isTRUE(verbose)) {
      message(sprintf("[%s] q=%.2f → n_trees=%d/%d", FUN, q, n_trees, Tm))
    }

    # Per-qi iterations
    for (r in seq_len(R)) {
      trees <- sort(sample.int(Tm, n_trees, replace = FALSE))
      fit   <- .fit_subset(trees)
      s_vec <- fit$scores

      # Gain curve
      gc <- .gain_curve(
        scores   = s_vec,
        n        = n_yvar_test,
        extr_idx = extr_idx_test,
        grid     = grid_eval
      )

      gc[, `:=`(
        q_trees = as.numeric(q),
        n_trees = as.integer(n_trees),
        rep     = as.integer(r)
      )]

      # Stash results
      out[[rr]] <- gc
      rr <- rr + 1L

      # Diagnostics vs full
      cor_vs_full  <-
        stats::cor(s_vec, s_full, method = "spearman", use = "pairwise.complete.obs")
      top_sub      <-
        utils::head(order(s_vec, decreasing = TRUE), k_top)
      jacc_vs_full <- .jacc(top_full, top_sub)

      diag_rows[[dr]] <- data.table::data.table(
        q_trees         = as.numeric(q),
        n_trees         = as.integer(n_trees),
        rep             = as.integer(r),
        cor_vs_full     = as.numeric(cor_vs_full),
        jaccard_vs_full = as.numeric(jacc_vs_full),
        trees_sampled   = paste0(fit$trees, collapse = ",")
      )
      dr <- dr + 1L

      # Save for within-q Jaccard across iterations
      top_sets[[qi]][[r]] <- top_sub

      # Optional progress update
      if (!is.null(progress_every) &&
          progress_every > 0L && (r %% progress_every == 0L)) {
        message(sprintf("[%s] q=%.2f → iteration %d/%d", FUN, q, r, R))
      }
    }
  }

  # Assemble outputs
  results   <-
    data.table::rbindlist(out,      use.names = TRUE, fill = TRUE)
  diag_runs <-
    data.table::rbindlist(diag_rows, use.names = TRUE, fill = TRUE)

  # Mean within-q Jaccard across iterations
  jacc_within_q <- {
    vals <- lapply(seq_along(subsample_grid), function(qi) {
      L  <- top_sets[[qi]]
      Rq <- length(L)
      if (Rq < 2L)
        return(NA_real_)
      idx <- utils::combn(Rq, 2L)
      mean(vapply(seq_len(ncol(idx)), function(j) {
        .jacc(L[[idx[1L, j]]], L[[idx[2L, j]]])
      }, numeric(1)), na.rm = TRUE)
    })
    data.table::data.table(
      q_trees = as.numeric(subsample_grid),
      n_trees = as.integer(pmax(1L, round(
        subsample_grid * Tm
      ))),
      jaccard = as.numeric(vals)
    )
  }

  # Summarize across iterations per q
  summary_core <- if (nrow(results)) {
    results[, .(
      lift_curve_mean = mean(lift_curve),
      lift_curve_sd   = sd(lift_curve),
      recall_mean     = mean(recall),
      recall_sd       = sd(recall)
    ), by = .(q_trees, n_trees)]
  } else
    data.table::data.table()

  cor_tbl <- if (nrow(diag_runs)) {
    diag_runs[, .(
      mean_cor_vs_full  = mean(cor_vs_full,     na.rm = TRUE),
      mean_jacc_vs_full = mean(jaccard_vs_full, na.rm = TRUE)
    ), by = .(q_trees, n_trees)]
  } else
    data.table::data.table()

  summary <- merge(
    summary_core,
    jacc_within_q,
    by = c("q_trees", "n_trees"),
    all.x = TRUE
  )
  summary <- merge(summary,
                   cor_tbl,
                   by = c("q_trees", "n_trees"),
                   all.x = TRUE)

  # gain curves for all qs and iterations
  list(results   = results[],
       # correlation and Jaccard vs full
       diag_runs = diag_runs[],
       # summarized metrics per q
       summary   = summary[])
}
