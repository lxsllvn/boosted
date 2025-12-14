#' Title
#'
#' @param boosted
#' @param harvest
#' @param candidate_rules
#' @param min_support
#' @param progress_every
#' @param return_ledger
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

validate_rules <- function(boosted,
                           harvest,
                           candidate_rules = NULL,
                           min_support     = 1L,
                           progress_every  = NULL,
                           return_ledger   = TRUE) {
  # Signature & basic checks
  FUN <- "validate_rules"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  if (!inherits(boosted, "boosted")) {
    stop(sprintf("[%s]  Input must be an object of class 'boosted' (from make_boosted()).",
                 FUN))
  }
  if (!inherits(harvest, "boosted_harvest")) {
    stop(
      sprintf(
        "[%s] harvest must be the object returned by harvest/validate_rules(..., return_ledger = TRUE).",
        FUN
      )
    )
  }

  # Pull test data from boosted
  extr_idx_test  <- boosted$extr_idx_test
  bg_idx_test    <- boosted$bg_idx_test
  N_extr_test    <- boosted$N_extr_test
  N_bg_test      <- boosted$N_bg_test
  Tm             <- boosted$Tm
  base_rate_test <- boosted$base_rate_test
  n_yvar_test    <- boosted$n_yvar_test
  test_leaf_map  <- boosted$test_leaf_map
  train_leaf_map <- boosted$train_leaf_map

  # Extract native and dense leaf IDs
  dense_leaf_ids_test <- test_leaf_map$dense_leaf_ids
  native_ids_train    <- train_leaf_map$native_leaf_ids

  # Optional yvar (for medians)
  yvar_test <- boosted$yvar_test

  # Harvested rules from training data
  pairs_all_train <- data.table::as.data.table(harvest$pairs_all)
  rule_meta_tbl   <- data.table::as.data.table(harvest$R[, .(rule_str, rule_len)])

  # Harvest metadata
  max_depth        <- harvest$meta$max_depth
  tighten_monotone <- harvest$meta$tighten_monotone
  n_train          <- harvest$meta$n_train

  # Decide which rules to validate
  if (is.null(candidate_rules)) {
    # Default: all rules from harvest
    candidate_vec <- unique(rule_meta_tbl$rule_str)
  } else if (is.character(candidate_rules)) {
    candidate_vec <- unique(candidate_rules)
  } else {
    cd <- data.table::as.data.table(candidate_rules)
    if (!"rule_str" %in% names(cd)) {
      stop(
        sprintf(
          "[%s] candidate_rules must have a column 'rule_str' (or be a character vector).",
          FUN
        )
      )
    }
    candidate_vec <- unique(cd$rule_str)
  }

  if (!length(candidate_vec)) {
    warning(sprintf("[%s] candidate_rules is empty; nothing to validate.", FUN))
    return(invisible(NULL))
  }

  # Filter training (Tree, leaf_id, rule_str) map down to candidate rules
  CR <- data.table::data.table(rule_str = candidate_vec)
  data.table::setkey(CR, rule_str)
  data.table::setkey(pairs_all_train, rule_str)
  data.table::setkey(rule_meta_tbl, rule_str)

  pairs_all <- pairs_all_train[CR, nomatch = 0L]
  if (!nrow(pairs_all)) {
    warning(
      sprintf(
        "[%s] none of the candidate rules were found in harvest_ledger$pairs_all.",
        FUN
      )
    )
    return(invisible(NULL))
  }


  # Build per-tree leaf → SNP maps, reused for all rules and prefixes.
  # For each tree tt and native leaf_id:
  # snps_ext_by_leaf[[tt]][[leaf_id]] = extreme SNP indices
  # snps_bg_by_leaf [[tt]][[leaf_id]] = background SNP indices
  # snps_all_by_leaf[[tt]][[leaf_id]] = all SNP indices (labeled + unlabeled)
  lookup <- .build_lookup(
    dense_leaf_ids  = dense_leaf_ids_test,
    native_leaf_ids = native_ids_train,
    extr_idx        = extr_idx_test,
    bg_idx          = bg_idx_test,
    all             = TRUE
  )
  snps_ext_by_leaf <- lookup$snps_ext_by_leaf
  snps_bg_by_leaf  <- lookup$snps_bg_by_leaf
  snps_all_by_leaf <- lookup$snps_all_by_leaf

  # Identify all unique rule strings to validate
  uniq_rules <- sort(unique(pairs_all$rule_str))
  pooled     <- vector("list", length(uniq_rules))

  # Main pooling over rules: evaluate each rule on the held-out test set
  for (i in seq_along(uniq_rules)) {
    rs <- uniq_rules[i]
    # All (Tree, leaf_id) pairs where this rule held in the training set
    pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
    if (!nrow(pairs))
      next

    # Look up SNP indices under this rule via the lookup
    buckets <- .snp_lookup(
      pairs            = pairs,
      snps_ext_by_leaf = snps_ext_by_leaf,
      snps_bg_by_leaf  = snps_bg_by_leaf,
      snps_all_by_leaf = snps_all_by_leaf
    )
    bucket_ext <- buckets$bucket_ext
    bucket_bg  <- buckets$bucket_bg
    bucket_all <- buckets$bucket_all

    # Counts: how many test extreme vs background SNPs the rule captures
    n_e         <- length(bucket_ext)
    n_b         <- length(bucket_bg)
    support     <- n_e + n_b          # labeled support
    support_all <- length(bucket_all) # total support (labeled + unlabeled)

    # Drop rules with labeled support < min_support
    if (support < as.integer(min_support)) {
      next
    }

    # Performance metrics on test set
    precision <- n_e / pmax(1e-12, support)
    recall    <- if (N_extr_test > 0L)
      n_e / N_extr_test
    else
      NA_real_
    lift <-
      if (base_rate_test > 0)
        precision / base_rate_test
    else
      NA_real_

    # Hypergeometric p-value under null (rule draws extremes at baseline rate)
    pval <- if (support > 0L) {
      stats::phyper(
        q = n_e - 1L,
        m = N_extr_test,
        n = N_bg_test,
        k = support,
        lower.tail = FALSE
      )
    } else {
      NA_real_
    }

    # Haldane–Anscombe corrected odds ratio (finite even with zeros)
    or_ha <- ((n_e + 0.5) * (N_bg_test - n_b + 0.5)) /
      ((N_extr_test - n_e + 0.5) * (n_b + 0.5))

    # Optional: medians of test yvars under this rule
    med_e <- if (!is.null(yvar_test) && n_e) {
      median(yvar_test[bucket_ext], na.rm = TRUE)
    } else {
      NA_real_
    }
    med_b <- if (!is.null(yvar_test) && n_b) {
      median(yvar_test[bucket_bg], na.rm = TRUE)
    } else {
      NA_real_
    }
    # Overall median: all SNPs under the rule (labeled + unlabeled)
    med_o <- if (!is.null(yvar_test) && length(bucket_all)) {
      median(yvar_test[bucket_all], na.rm = TRUE)
    } else {
      NA_real_
    }

    # Rule-level LLRs with Jeffreys smoothing on contingency counts
    enrichment <- .rule_llr(
      n_extreme    = n_e,
      n_bg         = n_b,
      N_extr_total = N_extr_test,
      N_bg_total   = N_bg_test
    )

    # Rule length from the training-time catalog
    rl <- rule_meta_tbl[rs, rule_len]

    pooled[[i]] <- data.table::data.table(
      rule_str      = rs,
      rule_len      = as.integer(rl),
      n_extreme     = n_e,
      n_bg          = n_b,
      support       = support,
      support_all   = support_all,
      precision     = precision,
      recall        = recall,
      lift          = lift,
      odds_ratio    = or_ha,
      enrichment    = enrichment,
      pval          = pval,
      med_y_extreme = med_e,
      med_y_bg      = med_b,
      med_y_overall = med_o
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
      "[%s] no pooled rules on the test set (after min_support filtering).",
      FUN
    ))
    return(invisible(NULL))
  }

  # q-values & ordering
  R_tbl[, qval := p.adjust(pval, method = "BH")]
  data.table::setorderv(
    R_tbl,
    cols  = c("qval", "lift", "odds_ratio", "recall", "precision"),
    order = c(1L,-1L,-1L,-1L,-1L)
  )

  if (isTRUE(return_ledger)) {
    out <- list(
      R            = R_tbl[],
      pairs_all    = pairs_all[],
      meta         = list(
        max_depth        = max_depth,
        tighten_monotone = tighten_monotone,
        n_train          = n_yvar_train,
        n_test           = n_yvar_test,
        Tm               = Tm
      )
    )
    # Tag as a rule-harvest ledger without interfering with other list-based tools
    class(out) <- c("boosted_harvest", class(out))
    return(out)
  } else {
    R_tbl[]
  }
}
