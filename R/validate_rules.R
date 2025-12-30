#' Title
#'
#' @param boosted
#' @param harvest
#' @param candidate_rules
#' @param min_support
#' @param progress_every
#' @param return_ledger
#' @param which
#' @param fold_indices
#' @param shrink_m
#' @param compute_pq
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

validate_rules <- function(boosted,
                           harvest,
                           which = c("train", "test"),
                           min_support     = 1L,
                           shrink_m        = 0,
                           candidate_rules = NULL,
                           fold_indices    = NULL,
                           compute_pq      = NULL,
                           progress_every  = NULL,
                           return_ledger   = FALSE) {
  # Signature & basic checks
  FUN <- "validate_rules"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  if (!inherits(boosted, "boosted")) {
    stop(sprintf("[%s]  Input must be an object of class 'boosted' (from make_boosted()).",
                 FUN))
  }
  if (!inherits(boosted, "boosted_binned")) {
    stop(sprintf(
      "[%s] boosted is not ready; run prepare_harvest() first.",
      FUN
    ))
  }
  if (!inherits(harvest, "boosted_harvest")) {
    stop(
      sprintf(
        "[%s] harvest must be the object returned by harvest/validate_rules(..., return_ledger = TRUE).",
        FUN
      )
    )
  }

  # Determine if we're working with test or training data
  which <- match.arg(which)
  # Compute p/q values only if validating on test data
  if (is.null(compute_pq)) {
    compute_pq <- identical(which, "test")
  }

  # Pull test/training indices and yvar data from `boosted` as requested
  extr_idx <- boosted[[sprintf("extr_idx_%s", which)]]
  bg_idx   <- boosted[[sprintf("bg_idx_%s",   which)]]
  yvar     <- boosted[[sprintf("yvar_%s",     which)]]
  n_all    <- boosted[[sprintf("n_yvar_%s",   which)]]

  # Pull per-tree leaf → SNP maps; reused for all rules and prefixes.
  snps_all_by_leaf <- boosted[[sprintf("snps_all_by_leaf_%s",   which)]]

  # If fold_indices != null, ensure that they are a valid vector of integers
  if (is.null(fold_indices)) {
    all_idx <- seq_len(n_all)
  } else {
    all_idx <- .check_idx(fold_indices, n_all, FUN, "fold_indices")
  }

  # Restrict background/extreme labels to fold_indices
  extr_idx  <- intersect(extr_idx, all_idx)
  bg_idx    <- intersect(bg_idx, all_idx)
  N_extr    <- length(extr_idx)
  N_bg      <- length(bg_idx)
  base_rate <- if ((N_extr + N_bg) > 0L) N_extr / (N_extr + N_bg) else NA_real_

  # Membership vectors for fast counting within buckets
  in_fold    <- rep(FALSE, n_all)
  is_extreme <- rep(FALSE, n_all)
  is_bg      <- rep(FALSE, n_all)

  in_fold[all_idx]     <- TRUE
  is_extreme[extr_idx] <- TRUE
  is_bg[bg_idx]        <- TRUE

  # Harvested rules
  pairs_all_ledger <- data.table::as.data.table(harvest$pairs_all)
  rule_meta_tbl    <- unique(pairs_all_ledger[, .(rule_str, rule_len)], by = "rule_str")
  data.table::setkey(rule_meta_tbl, rule_str)

  # Harvest metadata
  max_depth        <- harvest$meta$max_depth
  tighten_monotone <- harvest$meta$tighten_monotone
  n_train          <- harvest$meta$n_train
  Tm               <- harvest$meta$Tm

  # Decide which rules to validate
  if (is.null(candidate_rules)) {
    candidate_vec <- unique(rule_meta_tbl$rule_str)
  } else if (is.character(candidate_rules)) {
    candidate_vec <- unique(candidate_rules)
  } else {
    cd <- data.table::as.data.table(candidate_rules)
    if ("rule_str" %in% names(cd)) {
      candidate_vec <- unique(cd$rule_str)
    } else if ("rule_prefix" %in% names(cd)) {
      # allow this, but note: validate_rules can only score rules that exist in pairs_all ledger
      candidate_vec <- unique(cd$rule_prefix)
    } else {
      stop(sprintf("[%s] candidate_rules must be a character vector or have 'rule_str' (or 'rule_prefix').", FUN))
    }
  }
  candidate_vec <- candidate_vec[!is.na(candidate_vec)]
  candidate_vec <- candidate_vec[nzchar(trimws(candidate_vec))]
  if (!length(candidate_vec)) {
    warning(sprintf("[%s] candidate_rules is empty; nothing to validate.", FUN))
    return(invisible(NULL))
  }

  # Filter (Tree, leaf_id, rule_str) map down to candidate rules
  CR <- data.table::data.table(rule_str = candidate_vec)
  data.table::setkey(CR, rule_str)
  data.table::setkey(pairs_all_ledger, rule_str)
  data.table::setkey(rule_meta_tbl, rule_str)

  pairs_all <- pairs_all_ledger[CR, nomatch = 0L]
  if (!nrow(pairs_all)) {
    warning(
      sprintf(
        "[%s] none of the candidate rules were found in harvest$pairs_all.",
        FUN
      )
    )
    return(invisible(NULL))
  }

  # Identify all unique rule strings to evaluate
  uniq_rules <- sort(unique(pairs_all$rule_str))
  pooled     <- vector("list", length(uniq_rules))

  # Main pooling over rules
  for (i in seq_along(uniq_rules)) {
    rs    <- uniq_rules[i]
    # All (Tree, leaf_id) pairs under this rule
    pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
    if (!nrow(pairs))
      next

    # Look up SNP indices under this rule via the lookup
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

    # Optional beta–binomial smoothing (selection stability)
    precision_shrink <- NA_real_
    lift_shrink <- NA_real_
    if (is.numeric(shrink_m) && shrink_m > 0 && !is.na(base_rate) && base_rate > 0) {
      m <- as.numeric(shrink_m)
      precision_shrink <- (n_e + m * base_rate) / (support + m)
      lift_shrink <- precision_shrink / base_rate
    }

    # Hypergeometric p-value under null (rule draws extremes at baseline rate)
    # Default: only compute on test
    pval <- NA_real_
    qval <- NA_real_
    if (isTRUE(compute_pq)) {
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

    # Haldane–Anscombe corrected odds ratio (finite even with zeros)
    or_ha <- ((n_e + 0.5) * (N_bg - n_b + 0.5)) /
      ((N_extr - n_e + 0.5) * (n_b + 0.5))

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
      N_bg_total   = N_bg
    )

    # Rule length from the harvest catalog
    rl <- rule_meta_tbl[rs, rule_len]

    pooled[[i]] <- data.table::data.table(
      rule_str      = rs,
      rule_len      = as.integer(rl),

      n_extreme     = as.integer(n_e),
      n_bg          = as.integer(n_b),
      support       = as.integer(support),
      support_all   = as.integer(support_all),

      enrichment = as.numeric(enrichment),
      precision  = as.numeric(precision),
      recall     = as.numeric(recall),
      lift       = as.numeric(lift),
      odds_ratio = as.numeric(or_ha),
      pval       = as.numeric(pval),

      precision_shrink = as.numeric(precision_shrink),
      lift_shrink      = as.numeric(lift_shrink),

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

  # q-values
  if (isTRUE(compute_pq)) {
    R_tbl[, qval := stats::p.adjust(pval, method = "BH")]
  } else {
    R_tbl[, qval := NA_real_]
  }

  # Ordering: if qval exists use it; otherwise use lift (or lift_shrink if present)
  if (isTRUE(compute_pq)) {
    data.table::setorderv(
      R_tbl,
      cols = c("qval", "lift", "odds_ratio", "recall", "precision"),
      order = c(1L,-1L,-1L,-1L,-1L)
    )
  } else {
    ord_cols <-
      c("lift_shrink", "lift", "odds_ratio", "recall", "precision")
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
        base_rate   = base_rate,
        shrink_m    = shrink_m,
        max_depth   = max_depth,
        Tm          = Tm,
        n_train     = n_train,
        tighten_monotone = tighten_monotone
      )
    )

    # Tag as a rule-harvest ledger without interfering with other list-based tools
    class(out) <- c("boosted_harvest", class(out))
    return(out)
  } else {
    R_tbl[]
  }
}
