#' Score and rank candidate rules by enrichment for extreme SNPs
#'
#' For each candidate rule, pools all \code{(Tree, leaf_id)} pairs in the
#' ensemble where the rule occurs (via \code{boosted$leaf_rule_cache}), looks
#' up the SNPs covered by those leaves, and computes enrichment statistics
#' comparing the rule's extreme-to-background composition against the overall
#' base rate. Rules with labeled support below \code{min_support} are
#' silently dropped. By default, the terminal (deepest) rule of every leaf is
#' evaluated; a curated subset can be supplied via \code{candidate_rules}.
#'
#' @param min_support Positive integer. Minimum number of labelled SNPs
#'   (extreme + background) a rule must cover to be included in the output.
#' @param compute_pvals Logical or \code{NULL}. If \code{TRUE}, a one-sided
#'   hypergeometric p-value is computed for each rule under the null hypothesis
#'   that it captures extreme SNPs at the background base rate. \code{NULL}
#'   (default) computes p-values only when \code{which = "test"}.
#' @param return_ledger Logical. If \code{TRUE}, returns a named list
#'   containing the rule table (\code{R}), the full \code{pairs_all}
#'   \code{data.table} mapping rule strings to \code{(Tree, leaf_id)} pairs,
#'   and a \code{meta} list with fold-level counts. If \code{FALSE} (default),
#'   returns only the rule table.
#' @inheritParams .boosted_params
#'
#' @return When \code{return_ledger = FALSE} (default): a \code{data.table}
#'   with one row per rule string surviving the \code{min_support} filter,
#'   ordered by p-value (if computed) then by lift, recall, and precision.
#'   Columns:
#' \describe{
#'   \item{\code{rule_str}}{Character. The rule string.}
#'   \item{\code{n_clauses}}{Integer. Number of distinct split conditions in
#'     the rule after monotone tightening.}
#'   \item{\code{n_extreme}, \code{n_bg}}{Integer. Extreme and background SNP
#'     counts within the rule's SNP bucket.}
#'   \item{\code{support_labeled}}{Integer. Total labelled SNPs in the bucket
#'     (\code{n_extreme + n_bg}).}
#'   \item{\code{support_all}}{Integer. Total SNPs in the bucket (labelled and
#'     unlabelled).}
#'   \item{\code{enrichment}}{Numeric. Rule-level LLR with Jeffreys smoothing
#'     (see \code{.rule_llr}).}
#'   \item{\code{precision}}{Numeric. \code{n_extreme / support_labeled}.}
#'   \item{\code{recall}}{Numeric. \code{n_extreme / N_extr}.}
#'   \item{\code{lift}}{Numeric. \code{precision / base_rate}.}
#'   \item{\code{pval}}{Numeric. One-sided hypergeometric p-value; \code{NA}
#'     if \code{compute_pvals = FALSE}.}
#'   \item{\code{med_y_extreme}, \code{med_y_bg}, \code{med_y_overall}}{
#'     Numeric. Median response values for extreme, background, and all SNPs
#'     in the bucket.}
#' }
#'   When \code{return_ledger = TRUE}: a named list with elements \code{R}
#'   (the table above), \code{pairs_all}, and \code{meta}.
#' @import data.table
#' @export
validate_rules <- function(boosted,
                           which = c("train", "test"),
                           min_support     = 1L,
                           alpha           = 0.5,
                           candidate_rules = NULL,
                           fold_indices    = NULL,
                           compute_pvals   = NULL,
                           progress_every  = NULL,
                           return_ledger   = FALSE) {
  # Signature & basic checks
  FUN <- "validate_rules"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  if (!inherits(boosted, "boosted_harvest")) {
    stop(sprintf(
      "[%s] boosted is not ready; run prepare_harvest() first.",
      FUN
    ))
  }

  # Determine if we're working with test or training data
  which <- match.arg(which)

  # Compute p values only if validating on test data
  if (is.null(compute_pvals)) {
    compute_pvals <- identical(which, "test")
  }

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

  # If candidate_rules !=null, we check if the rule strings are found in
  # `boosted$leaf_rule_cache` and return the (Tree, leaf_id) rows under each
  # rule.
  # Otherwise, we pull all terminal rules from `boosted$leaf_rule_cache`, i.e.,
  # the longest rule for a (Tree, leaf_id), and their (Tree, leaf_id) pairs. .
  rules <- .validate_candidate_rules(
    candidate_rules = candidate_rules,
    leaf_rule_cache = leaf_rule_cache,
    caller          = FUN
  )

  uniq_rules  <- rules$candidate_rules
  pairs_all   <- rules$pairs_all
  pairs_by_rule <- split(pairs_all, by = "rule_str")

  # Pre-allocate output vectors — one slot per candidate rule.
  # Filled by index inside the loop; rows where keep[i] = FALSE are skipped
  # rules (empty bucket or below min_support) and are excluded when building
  # the final data.table.
  n_rules         <- length(uniq_rules)
  keep            <- logical(n_rules)

  out_rule_str    <- character(n_rules)
  out_n_clauses   <- integer(n_rules)
  out_n_e         <- integer(n_rules)
  out_n_b         <- integer(n_rules)
  out_support_lab <- integer(n_rules)
  out_support_all <- integer(n_rules)
  out_enrichment  <- numeric(n_rules)
  out_precision   <- numeric(n_rules)
  out_recall      <- numeric(n_rules)
  out_lift        <- numeric(n_rules)
  out_pval        <- rep(NA_real_, n_rules)
  out_med_e       <- numeric(n_rules)
  out_med_b       <- numeric(n_rules)
  out_med_o       <- numeric(n_rules)

  # For each rule string, pool all (Tree, leaf_ids) in the ensemble where it
  # occurs, then use the look-up to count the labeled (extreme, background) and
  # unlabeled SNPs in these leaves.
  for (i in seq_along(uniq_rules)) {
    rs        <- uniq_rules[i]
    pr        <- pairs_by_rule[[rs]]
    n_clauses <- pr$n_clauses[1L]

    # All (Tree, leaf_id) pairs under this rule in the ensemble
    if (!nrow(pr))
      next

    # Look up SNP indices covered by these (Tree, leaf_id) pairs.
    # .snp_lookup() handles deduplication (no double-counting SNPs
    # that land in the same rule via multiple trees).
    buckets <- .snp_lookup(
      pairs            = pr,
      snps_all_by_leaf = snps_all_by_leaf
    )
    bucket_all <- buckets$bucket_all

    # Restrict bucket to fold_indices
    if (length(bucket_all)) {
      bucket_all <- bucket_all[in_fold[bucket_all]]
    }
    n_all_rule <- length(bucket_all)
    if (!n_all_rule) next

    # Counts: how many extreme vs background SNPs the rule captures
    # (fold-restricted)
    n_e <- sum(is_extreme[bucket_all])
    n_b <- sum(is_bg[bucket_all])

    support     <- n_e + n_b  # labeled support
    support_all <- n_all_rule # total support (labeled + unlabeled)

    # Drop rules with labeled support < min_support
    if (support < as.integer(min_support)) {
      next
    }

    # Performance metrics
    precision <- n_e / pmax(1e-12, support)
    recall    <- if (N_extr > 0L)
      n_e / N_extr
    else
      NA_real_
    lift <-
      if (base_rate > 0)
        precision / base_rate
    else
      NA_real_

    # Hypergeometric p-value under null (rule draws extremes at baseline rate)
    # Default: only compute on test
    pval <- NA_real_
    if (isTRUE(compute_pvals)) {
      pval <- if (support > 0L) {
        stats::phyper(
          q = n_e - 1L,
          m = N_extr,
          n = N_bg,
          k = support,
          lower.tail = FALSE
        )
      } else {
        NA_real_
      }
    }

    # Medians of yvars under this rule (fold-restricted)
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

    # Store results in pre-allocated vectors
    keep[i]            <- TRUE
    out_rule_str[i]    <- rs
    out_n_clauses[i]   <- as.integer(n_clauses)
    out_n_e[i]         <- as.integer(n_e)
    out_n_b[i]         <- as.integer(n_b)
    out_support_lab[i] <- as.integer(support)
    out_support_all[i] <- as.integer(support_all)
    out_enrichment[i]  <- as.numeric(enrichment)
    out_precision[i]   <- as.numeric(precision)
    out_recall[i]      <- as.numeric(recall)
    out_lift[i]        <- as.numeric(lift)
    out_pval[i]        <- as.numeric(pval)
    out_med_e[i]       <- as.numeric(med_e)
    out_med_b[i]       <- as.numeric(med_b)
    out_med_o[i]       <- as.numeric(med_o)

    # Optional progress over rule space
    if (!is.null(progress_every) && progress_every > 0L &&
        (i %% progress_every == 0L)) {
      message(
        sprintf(
          "[%s] processed %d / %d rule strings (%.1f%%)",
          FUN,
          i,
          length(uniq_rules),
          100 * i / length(uniq_rules)
        )
      )
    }
  }

  # Build output table in one shot from pre-allocated vectors
  R_tbl <- data.table::data.table(
    rule_str        = out_rule_str[keep],
    n_clauses       = out_n_clauses[keep],
    n_extreme       = out_n_e[keep],
    n_bg            = out_n_b[keep],
    support_labeled = out_support_lab[keep],
    support_all     = out_support_all[keep],
    enrichment      = out_enrichment[keep],
    precision       = out_precision[keep],
    recall          = out_recall[keep],
    lift            = out_lift[keep],
    pval            = out_pval[keep],
    med_y_extreme   = out_med_e[keep],
    med_y_bg        = out_med_b[keep],
    med_y_overall   = out_med_o[keep]
  )
  if (!nrow(R_tbl)) {
    warning(sprintf(
      "[%s] no pooled rules on %s (after min_support filtering).",
      FUN,
      which
    ))
    return(invisible(NULL))
  }

  # Ordering: if pval exists use it; otherwise use lift
  if (isTRUE(compute_pvals)) {
    data.table::setorderv(
      R_tbl,
      cols = c("pval", "lift", "recall", "precision"),
      order = c(1L,-1L,-1L,-1L,-1L)
    )
  } else {
    ord_cols <-
      c("lift", "recall", "precision")
    ord_cols <- ord_cols[ord_cols %in% names(R_tbl)]
    data.table::setorderv(R_tbl,
                          cols = ord_cols,
                          order = rep(-1L, length(ord_cols)))
  }

  R_tbl <- .drop_all_na_cols(R_tbl)

  if (isTRUE(return_ledger)) {
    out <- list(
      R             = R_tbl[],
      pairs_all     = pairs_all[],
      meta          = list(
        which       = which,
        fold_n_all  = length(all_idx),
        fold_n_extr = N_extr,
        fold_n_bg   = N_bg,
        base_rate   = base_rate
      )
    )

  } else {
    R_tbl[]
  }
}
