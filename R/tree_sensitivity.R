#' Sensitivity of boosted scores to tree subsampling
#'
#' Assesses how much the gain curve and SNP score rankings depend on which
#' trees contribute to scoring. For each subsample fraction \code{q} in
#' \code{subsample_grid}, a random subset of \code{round(q * Tm)} trees is
#' drawn \code{R} times; leaf LLRs and SNP scores are recomputed from each
#' subset and compared against the full-model baseline via Spearman correlation
#' and Jaccard similarity of the top-\code{k} ranked SNPs. This helps diagnose
#' whether a small number of influential trees dominate the score or whether
#' signal is distributed across the ensemble.
#'
#' @param subsample_grid Numeric vector of tree subsample fractions in
#'   \code{(0, 1]}. Each value \code{q} determines the proportion of the
#'   \code{Tm} trees to retain in each replicate.
#' @param topk_frac Numeric scalar in \code{(0, 1]}. Fraction of test SNPs
#'   used to define the top-k set for Jaccard comparisons.
#' @inheritParams .boosted_params
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{results}}{A \code{data.table} with one row per
#'     \code{(q_trees, n_trees, rep, frac_screened)} combination, containing
#'     gain curve columns \code{recall}, \code{lift_curve}, and
#'     \code{score_threshold}.}
#'   \item{\code{diag_runs}}{A \code{data.table} with one row per replicate
#'     containing per-run diagnostics: \code{q_trees}, \code{n_trees},
#'     \code{rep}, \code{cor_vs_full} (Spearman correlation of subset scores
#'     vs full-model scores), \code{jaccard_vs_full} (Jaccard similarity of
#'     top-k sets), and \code{trees_sampled} (comma-separated tree indices).}
#'   \item{\code{summary}}{A \code{data.table} with one row per
#'     \code{(q_trees, n_trees)} summarising mean and SD of
#'     \code{lift_curve} and \code{recall} across replicates, plus
#'     \code{jaccard} (mean within-\code{q} Jaccard across replicate pairs)
#'     and \code{mean_cor_vs_full} and \code{mean_jacc_vs_full}.}
#' }
#' @export
#'
#' @examples
tree_sensitivity <- function(boosted,
                             subsample_grid = seq(0.10, 0.50, 0.10),
                             R              = 50L,
                             alpha          = 0.5,
                             method         = c("hard_label", "soft_label"),
                             model_llr_weight = 0.5,
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
  yvar_train     <- boosted$yvar_train
  extr_idx_test  <- boosted$extr_idx_test
  bg_idx_test    <- boosted$bg_idx_test
  n_yvar_test    <- boosted$n_yvar_test
  Tm             <- boosted$Tm
  train_leaf_map <- boosted$train_leaf_map
  test_leaf_map  <- boosted$test_leaf_map
  oof_prediction_train <- boosted$oof_prediction_train
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

  # Guard requested grid against duplicates and values == 0 or > 1
  grid_eval <- sort(unique(pmin(0.999, pmax(1e-3, gain_grid))))
  if (!length(grid_eval)) {
    stop(sprintf("[%s] gain_grid is empty after clipping/deduplication.", FUN))
  }

  k_top <- max(1L, as.integer(round(topk_frac * n_yvar_test)))

  soft_support <- NULL
  if (use_soft_labels) {
    # Tree subsampling changes only the tree set; soft support for the labelled
    # training SNPs is fixed and can be computed once.
    soft_support <- .soft_label_support(
      extr_idx = extr_idx_train,
      bg_idx = bg_idx_train,
      oof_prediction_train = oof_prediction_train,
      yvar_train = yvar_train,
      model_llr_weight = model_llr_weight,
      lower_tail = NULL
    )
  }

  # Helper: compute enrichment and score using a subset of trees
  .fit_subset <- function(trees) {
    trees <- sort(unique(as.integer(trees)))
    if (!length(trees))
      stop(sprintf("[%s] internal: empty tree subset.", FUN))

    if (use_hard_labels) {
      # Hard-label path: stack E/B counts on the sampled tree subset.
      ce_all <- .stack_leaf_counts(
        idx            = extr_idx_train,
        train_leaf_map = train_leaf_map,
        tree_idx       = trees,
        use_fast_counts = TRUE
      )
      cb_all <- .stack_leaf_counts(
        idx            = bg_idx_train,
        train_leaf_map = train_leaf_map,
        tree_idx       = trees,
        use_fast_counts = TRUE
      )
      enr <- .leaf_llrs(
        ce_all = ce_all,
        cb_all = cb_all,
        n_leaves = train_leaf_map$n_leaves,
        tree_idx = trees,
        N_E = N_extr_train,
        N_B = N_bg_train,
        alpha = alpha
      )
    } else {
      # Soft-label path: reuse fixed q values, but restack their masses on the
      # sampled tree subset.
      ce_all <- .stack_leaf_masses(
        idx = soft_support$labelled_idx,
        weights = soft_support$q,
        train_leaf_map = train_leaf_map,
        tree_idx = trees
      )
      cb_all <- .stack_leaf_masses(
        idx = soft_support$labelled_idx,
        weights = 1 - soft_support$q,
        train_leaf_map = train_leaf_map,
        tree_idx = trees
      )
      enr <- .leaf_llrs(
        ce_all = ce_all,
        cb_all = cb_all,
        n_leaves = train_leaf_map$n_leaves,
        tree_idx = trees,
        N_E = soft_support$diagnostics$n_extreme_effective,
        N_B = soft_support$diagnostics$n_background_effective,
        alpha = alpha
      )
      enr$diagnostics <- soft_support$diagnostics
    }
    sc <- .score_snps(
      test_leaf_map     = test_leaf_map,
      leaf_llrs_by_tree = enr,
      Tm                = length(trees),
      n                 = n_yvar_test
    )$scores
    list(scores = sc, trees = trees)
  }

  # Baseline: full-model scores for diagnostics
  full_fit <- .fit_subset(seq_len(Tm))
  s_full <- full_fit$scores

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
        bg_idx   = bg_idx_test,
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
      cor_vs_full <- suppressWarnings(stats::cor(
        s_vec,
        s_full,
        method = "spearman",
        use = "pairwise.complete.obs"
      ))
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
