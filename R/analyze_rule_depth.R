#' Title
#'
#' @param boosted
#' @param harvest
#' @param candidate_rules
#' @param which
#' @param shrink_m
#' @param alpha
#' @param fold_indices
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples
analyze_rule_depth <- function(boosted,
                               harvest,
                               candidate_rules,
                               which = c("train", "test"),
                               shrink_m       = 0,
                               alpha          = 0.5,
                               fold_indices   = NULL,
                               progress_every = NULL
                               ) {
  # Signature & basic checks
  FUN <- "analyze_rule_depth"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted")) {
    stop(sprintf(
      "[%s] Input must be an object of class 'boosted' (from make_boosted()).",
      FUN
    ))
  }
  if (!inherits(harvest, "boosted_harvest")) {
    stop(
      sprintf(
        "[%s] harvest must have ledger (harvest$pairs_all) created by harvest_rules(...) or validate_rules(..., return_ledger = TRUE).",
        FUN
      )
    )
  }

  # Validate candidate rules
  if (!is.character(candidate_rules)) {
    if (is.list(candidate_rules) || is.data.frame(candidate_rules)) {
      if ("rule_str" %in% names(candidate_rules)) {
        candidate_rules <- candidate_rules[["rule_str"]]
      } else if ("rule_prefix" %in% names(candidate_rules)) {
        candidate_rules <- candidate_rules[["rule_prefix"]]
      }
    }
  }
  if (!is.character(candidate_rules)) {
    stop(
      sprintf(
        "[%s] candidate_rules must be a character vector, or have a 'rule_str' or 'rule_prefix' field.",
        FUN
      )
    )
  }
  candidate_rules <-
    unique(candidate_rules[!is.na(candidate_rules)])
  candidate_rules <-
    candidate_rules[nzchar(trimws(candidate_rules))]
  if (!length(candidate_rules)) {
    stop(sprintf("[%s] candidate_rules is empty.",
                 FUN))
  }

  # Determine if we're working with test or training data
  which <- match.arg(which)

  # Pull test/training indices and yvar data from `boosted` as requested
  extr_idx <- boosted[[sprintf("extr_idx_%s", which)]]
  bg_idx   <- boosted[[sprintf("bg_idx_%s",   which)]]
  yvar     <- boosted[[sprintf("yvar_%s",     which)]]
  n_all    <- boosted[[sprintf("n_yvar_%s",   which)]]

  # Pull per-tree leaf → SNP maps; reused for all rules and prefixes.
  snps_all_by_leaf <- boosted[[sprintf("snps_all_by_leaf_%s",   which)]]
  leaf_paths       <- boosted$leaf_paths
  harvest_bins     <- boosted$harvest_bins

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

  # Harvest ledger pieces: (Tree, leaf_id, rule_str) mapping and metadata
  pairs_all        <- data.table::as.data.table(harvest$pairs_all)
  max_depth        <- as.integer(harvest$meta$max_depth)
  tighten_monotone <- isTRUE(harvest$meta$tighten_monotone)

  # Make sure max_depth exists
  if (is.na(max_depth) || max_depth < 1L) {
    stop("[analyze_rule_depth] harvest$meta$max_depth must be a positive integer")
  }

  anchor_rules <- pairs_all[rule_str %in% candidate_rules]
  rm(pairs_all)

  if (!nrow(anchor_tbl)) {
    stop(sprintf("[%s] No candidate_rules matched harvest$pairs_all$rule_str.", FUN))
  }
  anchor_tbl <- unique(anchor_rules, by = "rule_str")

  prefix_universe <- anchor_tbl[, {
    rs_full <- rule_str[1L]
    conds   <- strsplit(rs_full, " | ", fixed = TRUE)[[1]]
    n_cond  <- length(conds)

    max_d <- min(as.integer(max_depth), n_cond - 1L)
    if (!is.finite(max_d) || max_d < 1L) NULL else
      data.table::data.table(
        rule_str_full   = rs_full,
        full_rule_len   = as.integer(n_cond),
        prefix_len      = as.integer(seq_len(max_d)),
        rule_str_prefix = vapply(seq_len(max_d),
                                 function(d) paste(conds[1:d], collapse = " | "),
                                 character(1))
      )
  }, by = rule_str][, rule_str := NULL]

  if (!nrow(prefix_universe)) {
    warning(sprintf("[%s] no prefixes generated from anchor rules (check inputs).", FUN))
    return(list(prefixes = NULL, pairs_all = pairs_all[], map = NULL))
  }

  prefix_list <- prefix_universe[, .(
    n_full_rules_with_prefix = data.table::uniqueN(rule_str_full),
    min_full_rule_len        = min(full_rule_len),
    max_full_rule_len        = max(full_rule_len)
  ), by = .(prefix_len, rule_str_prefix)]

  data.table::setorderv(prefix_list,
                        cols  = c("prefix_len", "n_full_rules_with_prefix"),
                        order = c(1L, -1L))

  want_prefix <- unique(prefix_list[, .(prefix_len, rule_str_prefix)])
  want_prefix[, prefix_len := as.integer(prefix_len)]

  pairs_by_prefix <- .map_prefixes_to_leaves(
    want_prefix      = want_prefix,
    leaf_paths       = leaf_paths,
    harvest_bins     = harvest_bins,
    max_depth        = max_depth,
    tighten_monotone = tighten_monotone,
    proper_only      = TRUE
  )

  # Add prefix (Tree, leaf_id, rule_str) to candidate rule ledger for reuse
  pairs_all <- data.table::rbindlist(
    list(
      anchor_rules[, .(rule_str, Tree, leaf_id, rule_len)],
      pairs_by_prefix[, .(rule_str = rule_str_prefix, Tree, leaf_id, rule_len = prefix_len)]
    ),
    use.names = TRUE
  )
  pairs_all <- unique(pairs_all, by = c("rule_str", "Tree", "leaf_id"))
  data.table::setkey(pairs_all, rule_str)


  # Pre-allocate container for results (one per unique prefix)
  out <- vector("list", nrow(prefix_list))
  rr  <- 1L

  # Compute metrics per prefix
  for (i in seq_len(nrow(prefix_list))) {
    prefix_str <- prefix_list$rule_str_prefix[i]

    # All (Tree, leaf_id) pairs whose rule_str starts with this prefix
    sub <-
      pairs_by_prefix[.(prefix_list$prefix_len[i], prefix_str), .(Tree, leaf_id)]

    if (!nrow(sub)) {
      next
    }

    pref_buckets <- .snp_lookup(
      pairs            = sub,
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

    # Optional beta–binomial smoothing (selection stability)
    precision_shrink <- NA_real_
    lift_shrink <- NA_real_
    if (is.numeric(shrink_m) && shrink_m > 0 && !is.na(base_rate) && base_rate > 0) {
      m <- as.numeric(shrink_m)
      precision_shrink <- (n_e + m * base_rate) / (support + m)
      lift_shrink <- precision_shrink / base_rate
    }

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

    out[[rr]] <- data.table::data.table(
      rule_str_prefix = prefix_str,
      prefix_len      = as.integer(prefix_list$prefix_len[i]),

      n_extreme   = as.integer(n_e),
      n_bg        = as.integer(n_b),
      support     = as.integer(support),
      support_all = as.integer(support_all),

      enrichment = as.numeric(enrichment),
      precision  = as.numeric(precision),
      recall     = as.numeric(recall),
      lift       = as.numeric(lift),

      precision_shrink = as.numeric(precision_shrink),
      lift_shrink      = as.numeric(lift_shrink),

      med_y_extreme = as.numeric(med_e),
      med_y_bg      = as.numeric(med_b),
      med_y_overall = as.numeric(med_o),

      # Context: how often this prefix occurs among anchor rules
      n_rules_with_prefix      = as.integer(prefix_list$n_full_rules_with_prefix[i]),
      min_full_rule_len        = as.integer(prefix_list$min_full_rule_len[i]),
      max_full_rule_len        = as.integer(prefix_list$max_full_rule_len[i])
    )
    rr <- rr + 1L

    if (!is.null(progress_every) &&
        progress_every > 0L &&
        (i %% progress_every == 0L)) {
      message(sprintf(
        "[%s] processed %d / %d unique prefixes",
        FUN,
        i,
        nrow(prefix_list)
      ))
    }
  }

  if (rr == 1L) {
    warning(sprintf("[%s] no results generated (check inputs).", FUN))
    return(list(prefixes = NULL, pairs_all = pairs_all[], map = prefix_universe[]))
  }
  out <- out[seq_len(rr - 1L)]

  R_tbl <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  R_tbl <- .drop_all_na_cols(R_tbl)

  if (!nrow(R_tbl)) {
    warning(sprintf("[%s] no prefix diagnostics generated (check inputs).", FUN))
    return(list(prefixes = NULL, pairs_all = pairs_all[], map = prefix_universe[]))
  }

  data.table::setorderv(
    R_tbl,
    cols  = c("prefix_len", "enrichment", "support"),
    order = c(1L, -1L, -1L)
  )

  # Return list with prefix results, ledger with appended prefixes, prefix-full rule map,
  # and harvest metadata
  out <- list(
    R         = R_tbl[],
    pairs_all = pairs_all[],
    meta      = list(max_depth        = max_depth,
                     tighten_monotone = tighten_monotone),
    map       = prefix_universe[]
  )
  # Tag as a rule-harvest ledger without interfering with other list-based tools
  class(out) <- c("boosted_prefixes", "boosted_harvest", class(out))
  return(out)
}
