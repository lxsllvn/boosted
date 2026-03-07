#' Title
#'
#' @description
#'
#' @param boosted
#' @param which use the "train" or "test" data partition; "train" is used by default.
#' @param min_support
#' @param alpha
#' @param candidate_rules
#' @param fold_indices (optional) integer vector of SNP indices; if provided, the analysis is restricted to these SNPs. By default, all SNPs in "train" or "test" are used.
#' @param compute_pvals (optional) if TRUE, return hypergeometric p-values under the null hypothesis that a rule enriches for extreme SNPs at the baseline rate. By default, p-values are only returned if which = "test".
#' @param progress_every (optional) print an update message after every N evaluated rules.
#' @param return_ledger
#'
#' @return
#' @import data.table
#' @export
#'
#' @examples
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
  data.table::setindex(pairs_all, rule_str)

  # Container for results
  pooled     <- vector("list", length(uniq_rules))

  # For each rule string, pool all (Tree, leaf_ids) in the ensemble where it
  # occurs, then use the look-up to count the labeled (extreme, background) and
  # unlabeled SNPs in these leaves.
  for (i in seq_along(uniq_rules)) {
    rs    <- uniq_rules[i]
    n_clauses <- pairs_all[rule_str == rs, n_clauses][1L]

    # All (Tree, leaf_id) pairs under this rule in the ensemble
    pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
    if (!nrow(pairs))
      next

    # Look up SNP indices covered by these (Tree, leaf_id) pairs.
    # .snp_lookup() handles deduplication (no double-counting SNPs
    # that land in the same rule via multiple trees).
    buckets <- .snp_lookup(
      pairs            = pairs,
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

    pooled[[i]] <- data.table::data.table(
      rule_str        = rs,
      n_clauses       = as.integer(n_clauses),

      n_extreme       = as.integer(n_e),
      n_bg            = as.integer(n_b),
      support_labeled = as.integer(support),
      support_all     = as.integer(support_all),

      enrichment = as.numeric(enrichment),
      precision  = as.numeric(precision),
      recall     = as.numeric(recall),
      lift       = as.numeric(lift),
      pval       = as.numeric(pval),

      med_y_extreme = as.numeric(med_e),
      med_y_bg      = as.numeric(med_b),
      med_y_overall = as.numeric(med_o)
    )

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

  R_tbl <- data.table::rbindlist(pooled, use.names = TRUE, fill = TRUE)
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
