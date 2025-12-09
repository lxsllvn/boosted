#' Title
#'
#' @param boosted
#' @param harvest
#' @param max_depth
#' @param top_n
#' @param progress_every
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples

analyze_rule_depth <- function(boosted,
                               harvest,
                               max_depth      = NULL,
                               top_n          = 1000L,
                               progress_every = NULL,
                               alpha          = 0.5) {
  # Signature & basic checks
  FUN <- "analyze_rule_depth"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted")) {
    stop(
      sprintf(
        "[%s] Input must be an object of class 'boosted' (from make_boosted()).",
        FUN
      )
    )
  }
  if (!inherits(harvest, "boosted_harvest")) {
    stop(
      sprintf(
        "[%s] harvest must be the object returned by harvest/validate_rules(..., return_ledger = TRUE).",
        FUN
      )
    )
  }

  # If max_depth is NULL, default to the depth used in harvest_rules()
  if (is.null(max_depth)) {
    max_depth <- as.integer(harvest$meta$max_depth)
  } else {
    max_depth <- as.integer(max_depth)
  }

  # Pull data from boosted
  extr_idx_train <- boosted$extr_idx_train
  bg_idx_train   <- boosted$bg_idx_train
  N_extr_train   <- boosted$N_extr_train
  N_bg_train     <- boosted$N_bg_train
  Tm             <- boosted$Tm

  train_leaf_map <- boosted$train_leaf_map
  dense_leaf_ids <- train_leaf_map$dense_leaf_ids
  native_ids_all <- train_leaf_map$native_leaf_ids

  base_rate <- boosted$base_rate_train
  y_num     <- boosted$yvar_train

  # Harvest ledger pieces: rule table + mapping from rules to (Tree, leaf_id)
  R_tbl        <- data.table::as.data.table(harvest$R)
  pairs_all    <- data.table::as.data.table(harvest$pairs_all)

  # Build per-tree leaf → SNP maps, reused for all rules and prefixes.
  # For each tree tt and native leaf_id:
  # snps_ext_by_leaf[[tt]][[leaf_id]] = extreme SNP indices
  # snps_bg_by_leaf [[tt]][[leaf_id]] = background SNP indices
  # snps_all_by_leaf[[tt]][[leaf_id]] = all SNP indices (labeled + unlabeled)
  lookup <- .build_lookup(
    dense_leaf_ids  = dense_leaf_ids,
    native_leaf_ids = native_ids_all,
    extr_idx        = extr_idx_train,
    bg_idx          = bg_idx_train,
    all             = TRUE
  )
  snps_ext_by_leaf <- lookup$snps_ext_by_leaf
  snps_bg_by_leaf  <- lookup$snps_bg_by_leaf
  snps_all_by_leaf <- lookup$snps_all_by_leaf

  # Choose anchor rules (top_n by lift, odds_ratio, recall, precision)
  data.table::setorderv(
    R_tbl,
    cols  = c("lift", "odds_ratio", "recall", "precision"),
    order = c(-1L,-1L,-1L,-1L)
  )
  if (is.finite(top_n) && top_n > 0L && top_n < nrow(R_tbl)) {
    R_anchor <- R_tbl[seq_len(top_n)]
  } else {
    R_anchor <- R_tbl
  }

  # Pre-allocate container for results
  out <- vector("list", nrow(R_anchor) * max_depth)
  # Initialize counter
  rr  <- 1L

  # For each full rule (anchor rule),
  for (i in seq_len(nrow(R_anchor))) {
    rs_full <- R_anchor$rule_str[i]
    # All (Tree, leaf_id) pairs where this full rule held in the ensemble
    pairs_full <- unique(pairs_all[rule_str == rs_full, .(Tree, leaf_id)])
    if (!nrow(pairs_full))
      next

    # Look up SNP indices under the full rule via the lookup
    full_buckets <- .snp_lookup(
      pairs            = pairs_full,
      snps_ext_by_leaf = snps_ext_by_leaf,
      snps_bg_by_leaf  = snps_bg_by_leaf,
      snps_all_by_leaf = snps_all_by_leaf
    )
    bucket_ext_full <- full_buckets$bucket_ext
    bucket_bg_full  <- full_buckets$bucket_bg
    bucket_all_full <- full_buckets$bucket_all

    n_e_full     <- length(bucket_ext_full)
    n_b_full     <- length(bucket_bg_full)
    support_full <- n_e_full + n_b_full # labeled support

    # LLR for the full rule with Jeffreys smoothing on contingency counts
    enrichment_full <- .rule_llr(
      n_extreme    = n_e_full,
      n_bg         = n_b_full,
      N_extr_total = N_extr_train,
      N_bg_total   = N_bg_train,
      alpha        = alpha
    )

    # Parse full rule into its ordered list of conditions
    conds  <- strsplit(rs_full, " \\| ", fixed = FALSE)[[1]]
    n_cond <- length(conds)

    # Only proper prefixes: up to min(max_depth, n_cond - 1)
    max_d <- min(as.integer(max_depth), n_cond - 1L)
    if (!is.finite(max_d) || max_d < 1L) {
      if (!is.null(progress_every) && progress_every > 0L &&
          (i %% progress_every == 0L)) {
        message(sprintf(
          "[%s] processed %d / %d anchor rules",
          FUN,
          i,
          nrow(R_anchor)
        ))
      }
      next
    }

    # For each prefix of the rule (1 condition, 2 conditions, ..., max_d)
    for (d in seq_len(max_d)) {
      prefix_str <- paste(conds[1:d], collapse = " | ")

      # All (Tree, leaf_id) pairs whose rule_str starts with this prefix.
      # This uses the harvested (Tree, leaf_id, rule_str) catalog
      sub <- pairs_all[substr(rule_str, 1L, nchar(prefix_str)) == prefix_str,
                       .(Tree, leaf_id)]
      if (!nrow(sub))
        next
      # Deduplicate (Tree, leaf_id) pairs: the same leaf can appear multiple
      # times in pairs_all, but we only want to look up each leaf once.
      sub <- unique(sub)
      # SNP row indices captured by this prefix across all trees
      pref_buckets <- .snp_lookup(
        pairs            = sub,
        snps_ext_by_leaf = snps_ext_by_leaf,
        snps_bg_by_leaf  = snps_bg_by_leaf,
        snps_all_by_leaf = snps_all_by_leaf
      )
      bucket_ext <- pref_buckets$bucket_ext
      bucket_bg  <- pref_buckets$bucket_bg
      bucket_all <- pref_buckets$bucket_all

      # Compute how many extreme vs background SNPs the prefix captures
      n_e <- length(bucket_ext)
      n_b <- length(bucket_bg)

      # Performance metrics: quantify how well the prefix enriches for extreme SNPs
      support <- n_e + n_b # labeled support
      precision <- n_e / pmax(1e-12, support)
      recall <-
        if (N_extr_train > 0L)
          n_e / N_extr_train
      else
        NA_real_
      lift <- if (base_rate > 0)
        precision / base_rate
      else
        NA_real_

      # LLR for this prefix with Jeffreys smoothing on contingency counts
      enrichment <- .rule_llr(
        n_extreme    = n_e,
        n_bg         = n_b,
        N_extr_total = N_extr_train,
        N_bg_total   = N_bg_train,
        alpha        = alpha
      )

      # Overlap between prefix and full rule (extreme/background/all SNPs)
      j_ext <- .jacc(bucket_ext_full, bucket_ext)
      j_bg  <- .jacc(bucket_bg_full,  bucket_bg)
      j_all <- .jacc(bucket_all_full, bucket_all)

      # yvar medians under this prefix (extreme/background/all SNPs)
      med_e <- if (!is.null(y_num) && n_e) {
        median(y_num[bucket_ext], na.rm = TRUE)
      } else {
        NA_real_
      }
      med_b <- if (!is.null(y_num) && n_b) {
        median(y_num[bucket_bg], na.rm = TRUE)
      } else {
        NA_real_
      }
      med_o <- if (!is.null(y_num) && length(bucket_all) > 0L) {
        median(y_num[bucket_all], na.rm = TRUE)
      } else {
        NA_real_
      }

      # Save results for this rule prefix
      out[[rr]] <- data.table::data.table(
        rule_str_full   = rs_full,
        rule_str_prefix = prefix_str,
        prefix_len      = as.integer(d),
        full_rule_len   = as.integer(n_cond),

        n_extreme_prefix = as.integer(n_e),
        n_bg_prefix      = as.integer(n_b),
        support_prefix   = as.integer(support),

        enrichment_prefix = as.numeric(enrichment),
        enrichment_full   = as.numeric(enrichment_full),

        precision_prefix = as.numeric(precision),
        recall_prefix    = as.numeric(recall),
        lift_prefix      = as.numeric(lift),

        jaccard_ext_vs_full = as.numeric(j_ext),
        jaccard_bg_vs_full  = as.numeric(j_bg),
        jaccard_all_vs_full = as.numeric(j_all),

        med_y_extreme = med_e,
        med_y_bg      = med_b,
        med_y_overall = med_o
      )
      rr <- rr + 1L
    }

    # Optional per-rule progress
    if (!is.null(progress_every) && progress_every > 0L &&
        (i %% progress_every == 0L)) {
      message(sprintf("[%s] processed %d / %d anchor rules",
                      FUN,
                      i,
                      nrow(R_anchor)))
    }
  }

  res <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!nrow(res)) {
    warning(sprintf("[%s] no prefix diagnostics generated (check inputs).", FUN))
    return(invisible(NULL))
  }

  # Order results
  data.table::setorderv(res,
                        cols  = c("rule_str_full", "prefix_len"),
                        order = c(1L, 1L))
  res[]
}
