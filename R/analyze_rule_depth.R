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
    order = c(-1L, -1L, -1L, -1L)
  )
  if (is.finite(top_n) && top_n > 0L && top_n < nrow(R_tbl)) {
    R_anchor <- R_tbl[seq_len(top_n)]
  } else {
    R_anchor <- R_tbl
  }

  # Derive prefixes from the selected anchor full rules, then compute
  # prefix statistics ONCE per unique (prefix_len, rule_str_prefix).
  # Make per-row grouping explicit (avoids .I confusion)
  R_anchor[, row_id := .I]

  prefix_universe <- R_anchor[, {
    rs_full <- rule_str
    conds   <- strsplit(rs_full, " \\| ", fixed = FALSE)[[1]]
    n_cond  <- length(conds)

    max_d <- min(as.integer(max_depth), n_cond - 1L)

    if (!is.finite(max_d) || max_d < 1L) {
      NULL
    } else {
      data.table::data.table(
        rule_str_full   = rs_full,
        full_rule_len   = as.integer(n_cond),
        prefix_len      = as.integer(seq_len(max_d)),
        rule_str_prefix = vapply(
          seq_len(max_d),
          function(d) paste(conds[1:d], collapse = " | "),
          character(1)
        )
      )
    }
  }, by = row_id]

  R_anchor[, row_id := NULL]


  if (!nrow(prefix_universe)) {
    warning(sprintf(
      "[%s] no prefixes generated from anchor rules (check inputs).",
      FUN
    ))
    return(invisible(NULL))
  }

  # One row per unique prefix; add a little context about how often it appears
  prefix_list <- prefix_universe[, .(
    n_full_rules_with_prefix = data.table::uniqueN(rule_str_full),
    min_full_rule_len        = min(full_rule_len),
    max_full_rule_len        = max(full_rule_len),
    example_full_rule        = rule_str_full[1]
  ), by = .(prefix_len, rule_str_prefix)]

  data.table::setorderv(
    prefix_list,
    cols  = c("prefix_len", "n_full_rules_with_prefix"),
    order = c(1L,-1L)
  )

  # Pre-allocate container for results (one per unique prefix)
  out <- vector("list", nrow(prefix_list))
  rr  <- 1L

  # Compute metrics per prefix
  for (i in seq_len(nrow(prefix_list))) {
    prefix_str <- prefix_list$rule_str_prefix[i]

    # All (Tree, leaf_id) pairs whose rule_str starts with this prefix
    sub <-
      pairs_all[substr(rule_str, 1L, nchar(prefix_str)) == prefix_str,
                .(Tree, leaf_id)]
    if (!nrow(sub)) {
      next
    }
    sub <- unique(sub)

    pref_buckets <- .snp_lookup(
      pairs            = sub,
      snps_ext_by_leaf = snps_ext_by_leaf,
      snps_bg_by_leaf  = snps_bg_by_leaf,
      snps_all_by_leaf = snps_all_by_leaf
    )
    bucket_ext <- pref_buckets$bucket_ext
    bucket_bg  <- pref_buckets$bucket_bg
    bucket_all <- pref_buckets$bucket_all

    n_e <- length(bucket_ext)
    n_b <- length(bucket_bg)

    support   <- n_e + n_b
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

    enrichment <- .rule_llr(
      n_extreme    = n_e,
      n_bg         = n_b,
      N_extr_total = N_extr_train,
      N_bg_total   = N_bg_train,
      alpha        = alpha
    )

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
    med_o <-
      if (!is.null(y_num) && length(bucket_all) > 0L) {
        median(y_num[bucket_all], na.rm = TRUE)
      } else {
        NA_real_
      }

    out[[rr]] <- data.table::data.table(
      rule_str_prefix = prefix_str,
      prefix_len      = as.integer(prefix_list$prefix_len[i]),

      n_extreme_prefix = as.integer(n_e),
      n_bg_prefix      = as.integer(n_b),
      support_prefix   = as.integer(support),

      enrichment_prefix = as.numeric(enrichment),
      precision_prefix  = as.numeric(precision),
      recall_prefix     = as.numeric(recall),
      lift_prefix       = as.numeric(lift),

      med_y_extreme = med_e,
      med_y_bg      = med_b,
      med_y_overall = med_o,

      # Context: how often this prefix occurs among anchor rules
      n_full_rules_with_prefix = as.integer(prefix_list$n_full_rules_with_prefix[i]),
      min_full_rule_len        = as.integer(prefix_list$min_full_rule_len[i]),
      max_full_rule_len        = as.integer(prefix_list$max_full_rule_len[i]),
      example_full_rule        = as.character(prefix_list$example_full_rule[i])
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

  res <-
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!nrow(res)) {
    warning(sprintf("[%s] no prefix diagnostics generated (check inputs).", FUN))
    return(invisible(NULL))
  }

  data.table::setorderv(
    res,
    cols  = c("prefix_len", "enrichment_prefix", "support_prefix"),
    order = c(1L,-1L,-1L)
  )
  res[]
}
