#' Evaluate candidate rules and all their prefixes across depths
#'
#' For each rule string in \code{candidate_rules}, retrieves all of its
#' prefixes (i.e. shorter rules formed by the first \eqn{1, 2, \ldots,
#' \text{rule\_len}}{} splits of each leaf path that produced the candidate)
#' from \code{boosted$leaf_rule_cache}. Each prefix is then scored in the same
#' way as \code{\link{validate_rules}}: SNPs covered by matching
#' \code{(Tree, leaf_id)} pairs are pooled and enrichment statistics are
#' computed. The result makes it straightforward to see how performance changes
#' as splits are added along each candidate's path, which is useful for
#' pruning overly deep rules or identifying the split that contributes most to
#' enrichment.
#'
#' @inheritParams .boosted_params
#'
#' @return A \code{data.table} with one row per unique rule string (candidate
#'   or prefix) that covers at least one SNP, ordered by \code{n_clauses}
#'   ascending then \code{enrichment} and \code{support_labeled} descending.
#'   Columns:
#' \describe{
#'   \item{\code{rule_str}}{Character. The rule or prefix string.}
#'   \item{\code{n_clauses}}{Integer. Number of split conditions after
#'     monotone tightening.}
#'   \item{\code{n_extreme}, \code{n_bg}}{Integer. Extreme and background SNP
#'     counts within the rule's SNP bucket.}
#'   \item{\code{support_labeled}}{Integer. Total labelled SNPs in the bucket.}
#'   \item{\code{support_all}}{Integer. Total SNPs in the bucket (labelled and
#'     unlabelled).}
#'   \item{\code{enrichment}}{Numeric. Rule-level LLR with Jeffreys smoothing.}
#'   \item{\code{precision}}{Numeric. \code{n_extreme / support_labeled}.}
#'   \item{\code{recall}}{Numeric. \code{n_extreme / N_extr}.}
#'   \item{\code{lift}}{Numeric. \code{precision / base_rate}.}
#'   \item{\code{med_y_extreme}, \code{med_y_bg}, \code{med_y_overall}}{
#'     Numeric. Median response values for extreme, background, and all SNPs
#'     in the bucket.}
#'   \item{\code{n_leaf_occurrences}}{Integer. Number of \code{(Tree, leaf_id)}
#'     pairs in the ensemble that match this rule.}
#'   \item{\code{n_trees}}{Integer. Number of distinct trees in which the rule
#'     occurs.}
#'   \item{\code{frac_trees}}{Numeric. \code{n_trees / Tm}.}
#' }
#' @export
#'
#' @examples

analyze_rule_depth <- function(boosted,
                               candidate_rules,
                               which = c("train", "test"),
                               alpha          = 0.5,
                               fold_indices   = NULL,
                               progress_every = NULL
) {
  # Signature & basic checks
  FUN <- "analyze_rule_depth"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  if (!inherits(boosted, "boosted_harvest")) {
    stop(sprintf(
      "[%s] boosted is not ready; run prepare_harvest() first.",
      FUN
    ))
  }

  # Determine if we're working with test or training data
  which <- match.arg(which)

  # Pull test/training indices and yvar data from `boosted` as requested
  extr_idx  <- boosted[[sprintf("extr_idx_%s",  which)]]
  bg_idx    <- boosted[[sprintf("bg_idx_%s",    which)]]
  yvar      <- boosted[[sprintf("yvar_%s",      which)]]
  n_all     <- boosted[[sprintf("n_yvar_%s",    which)]]
  N_extr    <- boosted[[sprintf("N_extr_%s",    which)]]
  N_bg      <- boosted[[sprintf("N_bg_%s",      which)]]
  base_rate <- boosted[[sprintf("base_rate_%s", which)]]

  # Pull (Tree, leaf_id) → SNP map
  snps_all_by_leaf <- boosted[[sprintf("snps_all_by_leaf_%s",   which)]]
  # Pull rule → (Tree, leaf_id) map
  leaf_rule_cache <- boosted$leaf_rule_cache

  # If fold_indices != null, ensure that they are a valid vector of integers
  if (is.null(fold_indices)) {
    all_idx <- seq_len(n_all)
  } else {
    all_idx <- .check_idx(fold_indices, n_all, FUN, "fold_indices")
    # Restrict background/extreme labels to fold_indices
    extr_idx  <- intersect(extr_idx, all_idx)
    bg_idx    <- intersect(bg_idx, all_idx)
    N_extr    <- length(extr_idx)
    N_bg      <- length(bg_idx)
    n_all     <- length(all_idx)
    base_rate <- if ((N_extr + N_bg) > 0L) N_extr / (N_extr + N_bg) else NA_real_
  }

  # Membership vectors for fast counting within buckets
  in_fold    <- rep(FALSE, n_all)
  is_extreme <- rep(FALSE, n_all)
  is_bg      <- rep(FALSE, n_all)

  in_fold[all_idx]     <- TRUE
  is_extreme[extr_idx] <- TRUE
  is_bg[bg_idx]        <- TRUE

  # Check if the candidate rule strings are found in `boosted$leaf_rule_cache`
  # and return the (Tree, leaf_id) rows under each rule.
  rules <- .validate_candidate_rules(
    candidate_rules = candidate_rules,
    leaf_rule_cache = leaf_rule_cache,
    caller          = FUN
  )
  pairs_all <- rules$pairs_all

  # Get (Tree, leaf_id) pairs for candidate rules and their prefixes
  pairs0 <- unique(pairs_all[, list(Tree = Tree,
                                    leaf_id = leaf_id,
                                    rule_len_full = rule_len)],
                   by = c("Tree", "leaf_id", "rule_len_full"))

  # Exclude rules longer than the candidate rule (cap by rule_len_full per (Tree, leaf_id))
  rules_to_score <- unique(leaf_rule_cache[pairs0,
                                           on = c("Tree", "leaf_id"),
                                           nomatch = 0L][rule_len <= rule_len_full,
                                                         rule_str])

  # Ensemble leaf membership for each rule_str (union over depths) and ensemble
  # metadata, computed once
  # membership: one row per (rule_str, Tree, leaf_id)
  ens_pairs <- unique(leaf_rule_cache[rule_str %in% rules_to_score,
                                      list(
                                        rule_str = rule_str,
                                        Tree = Tree,
                                        leaf_id = leaf_id,
                                        n_clauses = n_clauses
                                      )],
                      by = c("rule_str", "Tree", "leaf_id"))
  data.table::setindex(ens_pairs, rule_str)

  ens_meta <- ens_pairs[, list(
    n_leaf_occurrences = .N,
    n_trees            = data.table::uniqueN(Tree),
    n_clauses          = n_clauses[[1L]]
  ), by = list(rule_str)]
  data.table::setindex(ens_meta, rule_str)

  n_trees_total <- data.table::uniqueN(leaf_rule_cache$Tree)

  # Container for results
  pooled <- vector("list", length(rules_to_score))

  # For each rule string, pool all (Tree, leaf_ids) in the ensemble where it
  # occurs, then use the look-up to count the labeled (extreme, background) and
  # unlabeled SNPs in these leaves.
  for (i in seq_along(rules_to_score)) {
    rs <- rules_to_score[[i]]

    meta <- ens_meta[rule_str == rs]
    if (!nrow(meta)) next

    # All (Tree, leaf_id) pairs under this rule in the ensemble
    pairs <- ens_pairs[rule_str == rs, list(Tree = Tree, leaf_id = leaf_id)]
    if (!nrow(pairs)) next

    # Look up SNP indices covered by these (Tree, leaf_id) pairs.
    # .snp_lookup() handles deduplication (no double-counting SNPs
    # that land in the same rule via multiple trees).
    pref_buckets <- .snp_lookup(
      pairs            = pairs,
      snps_all_by_leaf = snps_all_by_leaf
    )
    bucket_all <- pref_buckets$bucket_all

    # Restrict bucket to fold_indices
    if (length(bucket_all)) {
      bucket_all <- bucket_all[in_fold[bucket_all]]
    }
    n_all_rule <- length(bucket_all)
    if (!n_all_rule) next

    # Counts: how many extreme vs background SNPs the prefix captures
    # (fold-restricted)
    n_e <- sum(is_extreme[bucket_all])
    n_b <- sum(is_bg[bucket_all])

    support     <- n_e + n_b  # labeled support
    support_all <- n_all_rule # total support (labeled + unlabeled)

    # Performance metrics
    precision   <- n_e / pmax(1e-12, support)
    recall <-
      if (N_extr > 0L)
        n_e / N_extr
    else
      NA_real_
    lift <- if (base_rate > 0)
      precision / base_rate
    else
      NA_real_

    # Medians of yvars under this prefix (fold-restricted)
    bucket_ext <- if (n_e) bucket_all[is_extreme[bucket_all]] else integer()
    bucket_bg  <- if (n_b) bucket_all[is_bg[bucket_all]] else integer()

    med_e <- if (length(bucket_ext)) median(yvar[bucket_ext], na.rm = TRUE) else NA_real_
    med_b <- if (length(bucket_bg)) median(yvar[bucket_bg], na.rm = TRUE) else NA_real_
    med_o <- if (length(bucket_all)) median(yvar[bucket_all], na.rm = TRUE) else NA_real_

    # Rule-level LLRs with Jeffreys smoothing on contingency counts
    enrichment <- .rule_llr(
      n_extreme    = n_e,
      n_bg         = n_b,
      N_extr_total = N_extr,
      N_bg_total   = N_bg,
      alpha        = alpha
    )

    # Add to results
    pooled[[i]] <- data.table::data.table(
      rule_str        = rs,
      n_clauses       = as.integer(meta$n_clauses[1L]),

      n_extreme       = as.integer(n_e),
      n_bg            = as.integer(n_b),
      support_labeled = as.integer(support),
      support_all     = as.integer(support_all),

      enrichment = as.numeric(enrichment),
      precision  = as.numeric(precision),
      recall     = as.numeric(recall),
      lift       = as.numeric(lift),

      med_y_extreme = as.numeric(med_e),
      med_y_bg      = as.numeric(med_b),
      med_y_overall = as.numeric(med_o),

      # Ensemble metadata
      n_leaf_occurrences = as.integer(meta$n_leaf_occurrences[[1L]]),
      n_trees            = as.integer(meta$n_trees[[1L]]),
      frac_trees         = as.numeric(meta$n_trees[[1L]] / n_trees_total)
    )

    # Optional progress over rule space
    if (!is.null(progress_every) && progress_every > 0L &&
        (i %% progress_every == 0L)) {
      message(
        sprintf(
          "[%s] processed %d / %d rule strings (%.1f%%)",
          FUN,
          i,
          length(rules_to_score),
          100 * i / length(rules_to_score)
        )
      )
    }
  }
  R_tbl <- data.table::rbindlist(pooled, use.names = TRUE, fill = TRUE)

  if (!nrow(R_tbl)) {
    warning(sprintf(
      "[%s] no prefix diagnostics generated (check inputs).", FUN))
    return(invisible(NULL))
  }

  data.table::setorderv(
    R_tbl,
    cols  = c("n_clauses", "enrichment", "support_labeled"),
    order = c(1L, -1L, -1L)
  )

  R_tbl <- .drop_all_na_cols(R_tbl)
  R_tbl[]
}
