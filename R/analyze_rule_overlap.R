#' Title
#'
#' @param boosted
#' @param harvest
#' @param candidate_rules
#' @param which character; evaluate on "train" or "test"
#' @param fold_indices optional integer vector of SNP indices (1-based) to restrict evaluation universe
#' @param shrink_m non-negative numeric; beta–binomial smoothing strength for precision/lift (0 disables)
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples
analyze_rule_overlap <- function(boosted,
                                 harvest,
                                 candidate_rules,
                                 which = c("train", "test"),
                                 fold_indices = NULL,
                                 shrink_m = 0,
                                 progress_every = NULL) {
  # Signature & basic checks
  FUN <- "analyze_rule_overlap"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted")) {
    stop(
      sprintf(
        "[%s]  Input must be an object of class 'boosted' (from make_boosted()).",
        FUN
      )
    )
  }
  if (!inherits(boosted, "boosted_binned")) {
    stop(sprintf("[%s] boosted is not ready; run prepare_harvest() first.",
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
  candidate_rules <- unique(candidate_rules[!is.na(candidate_rules)])
  candidate_rules <- candidate_rules[nzchar(trimws(candidate_rules))]

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
  n_tot    <- boosted[[sprintf("n_yvar_%s",   which)]]

  # Pull per-tree leaf → SNP maps; reused for all rules and prefixes.
  snps_all_by_leaf <- boosted[[sprintf("snps_all_by_leaf_%s",   which)]]
  leaf_paths       <- boosted$leaf_paths
  harvest_bins     <- boosted$harvest_bins

  # If fold_indices != null, ensure that they are a valid vector of integers
  if (is.null(fold_indices)) {
    all_idx <- seq_len(n_tot)
  } else {
    all_idx <- .check_idx(fold_indices, n_tot, FUN, "fold_indices")
  }

  # Restrict background/extreme labels to fold_indices
  extr_idx  <- intersect(extr_idx, all_idx)
  bg_idx    <- intersect(bg_idx, all_idx)
  N_extr    <- length(extr_idx)
  N_bg      <- length(bg_idx)
  base_rate <-
    if ((N_extr + N_bg) > 0L)
      N_extr / (N_extr + N_bg)
  else
    NA_real_

  # Membership vectors for fast counting within buckets
  in_fold    <- rep(FALSE, n_tot)
  is_extreme <- rep(FALSE, n_tot)
  is_bg      <- rep(FALSE, n_tot)

  in_fold[all_idx]     <- TRUE
  is_extreme[extr_idx] <- TRUE
  is_bg[bg_idx]        <- TRUE

  # Harvest ledger pieces: rule table + mapping from rules to (Tree, leaf_id)
  R_tbl     	     <- data.table::as.data.table(harvest$R)
  pairs_tbl 	     <- data.table::as.data.table(harvest$pairs_all)
  max_depth        <- as.integer(harvest$meta$max_depth)
  tighten_monotone <- isTRUE(harvest$meta$tighten_monotone)

  # Full rules among candidate_rules
  R_tbl <- R_tbl[rule_str %chin% candidate_rules]

  # Keep what we already have mapped in the ledger for these candidate strings
  pairs_have <-
    pairs_tbl[rule_str %chin% candidate_rules, .(rule_str, Tree, leaf_id)]
  missing_rules <-
    setdiff(candidate_rules, unique(pairs_have$rule_str))

  # Default: what we have
  pairs_all <- pairs_have

  if (length(missing_rules)) {
    want_prefix <-
      data.table::data.table(rule_str_prefix = missing_rules)
    want_prefix[, prefix_len := lengths(strsplit(rule_str_prefix, " | ", fixed = TRUE))]
    want_prefix <- want_prefix[prefix_len <= max_depth,
                               .(prefix_len = as.integer(prefix_len), rule_str_prefix)]
    want_prefix <-
      unique(want_prefix, by = c("prefix_len", "rule_str_prefix"))

    if (nrow(want_prefix)) {
      pairs_need <- .map_prefixes_to_leaves(
        want_prefix      = want_prefix,
        leaf_paths       = leaf_paths,
        harvest_bins     = harvest_bins,
        max_depth        = max_depth,
        tighten_monotone = tighten_monotone,
        proper_only      = TRUE
      )[, .(rule_str = rule_str_prefix, Tree, leaf_id)]

      pairs_all <- data.table::rbindlist(list(pairs_have, pairs_need),
                                         use.names = TRUE,
                                         fill = TRUE)
      pairs_all <- unique(pairs_all, by = c("rule_str", "Tree", "leaf_id"))
    }
  }

  data.table::setkey(pairs_all, rule_str)

  # Rule universe for overlap
  rule_ids <- unique(pairs_all$rule_str)
  rule_ids <- sort(rule_ids)
  K <- length(rule_ids)

  if (K < 2L) {
    warning(sprintf("[%s] fewer than two rules selected; nothing to compare.", FUN))
    return(
      list(
        overlap  = data.table::data.table(),
        summary  = data.table::data.table(),
        rule_ids = rule_ids,
        jaccard_ext = NULL,
        jaccard_bg  = NULL,
        jaccard_all = NULL
      )
    )
  }

  # Map rule_str -> index (1..K) for compact matrix operations
  ridx <- data.table::data.table(rule_str = rule_ids, i_index = seq_len(K))
  data.table::setkey(ridx, rule_str)
  data.table::setkey(pairs_all, rule_str)

  # Precompute buckets per rule (fold-restricted) and counts
  bucket_ext <- vector("list", K)
  bucket_bg  <- vector("list", K)
  bucket_all <- vector("list", K)
  n_all <- integer(K)
  n_ext <- integer(K)
  n_bg  <- integer(K)

  # Compute metrics per prefix
  for (k in seq_len(K)) {
    rs <- rule_ids[k]

    # All (Tree, leaf_id) pairs under this rule
    pairs <-
      unique(pairs_all[J(rs), .(Tree, leaf_id)], by = c("Tree", "leaf_id"))

    # Look up SNP indices covered by these (Tree, leaf_id) pairs.
    # .snp_lookup() handles deduplication (no double-counting SNPs
    # that land in the same rule via multiple trees).
    b <- .snp_lookup(pairs = pairs,
                     snps_all_by_leaf = snps_all_by_leaf)$bucket_all

    # Restrict bucket to fold_indices
    if (length(b))
      b <- b[in_fold[b]]
    bucket_all[[k]] <- b

    # Split buckets into extreme/background subsets (fold-restricted)
    # Needed for tcrossprod-based intersections
    if (length(b)) {
      bucket_ext[[k]] <- b[is_extreme[b]]
      bucket_bg[[k]]  <- b[is_bg[b]]
    } else {
      bucket_ext[[k]] <- integer(0)
      bucket_bg[[k]]  <- integer(0)
    }

    # Optional progress over rule buckets
    if (!is.null(progress_every) && progress_every > 0L &&
        (k %% progress_every == 0L)) {
      message(sprintf("[%s] built SNP buckets for %d / %d rules",
                      FUN, k, K))
    }
  }

  # Compute pairwise intersections in compiled code via sparse incidence
  # matrices and tcrossprod(). This replaces the O(K^2) pairwise
  # set-intersection loops.
  #
  # Note: SNP indices in buckets_* are row positions in the train/test leaf
  # matrices, so they are not contiguous. We remap them to a compact 1..N_used
  # index first.

  # Build a compact SNP universe over the selected rules (all buckets)
  snps_universe <-
    sort.int(unique(unlist(bucket_all, use.names = FALSE)))
  N_used <- length(snps_universe)

  # Guard: if all selected rules are empty, return empty structures
  if (N_used == 0L) {
    labs <- paste0("r", seq_len(K))

    J_ext <-
      matrix(0,
             nrow = K,
             ncol = K,
             dimnames = list(labs, labs))
    J_bg  <-
      matrix(0,
             nrow = K,
             ncol = K,
             dimnames = list(labs, labs))
    J_all <-
      matrix(0,
             nrow = K,
             ncol = K,
             dimnames = list(labs, labs))

    summary_tbl <- data.table::data.table()
    overlap_tbl <- data.table::data.table()

    return(
      list(
        summary     = summary_tbl[],
        jaccard_ext = J_ext,
        jaccard_bg  = J_bg,
        jaccard_all = J_all,
        overlap     = overlap_tbl[],
        rule_ids    = rule_ids
      )
    )
  } else {
    # Remap each bucket to 1..N_used (based on snps_universe)
    buckets_all_m <-
    lapply(bucket_all, function(b)
        match(b, snps_universe))
    buckets_ext_m <-
    lapply(bucket_ext, function(b)
        match(b, snps_universe))
    buckets_bg_m  <-
    lapply(bucket_bg,  function(b)
        match(b, snps_universe))

    # Build sparse incidence matrices: rows = rules, cols = SNPs (remapped
    # 1..N_used)
    .build_incidence <- function(buckets_m) {
      lens <- lengths(buckets_m)
      if (!any(lens)) {
        return(Matrix::sparseMatrix(
          i = integer(0),
          j = integer(0),
          x = 1L,
          dims = c(K, N_used)
        ))
      }
      i <- rep.int(seq_len(K), lens)
      j <- unlist(buckets_m, use.names = FALSE)
      Matrix::sparseMatrix(
        i = i,
        j = j,
        x = 1L,
        dims = c(K, N_used)
      )
    }

    M_all <- .build_incidence(buckets_all_m)
    M_ext <- .build_incidence(buckets_ext_m)
    M_bg  <- .build_incidence(buckets_bg_m)

    # Intersection count matrices (K x K)
    I_all <- as.matrix(Matrix::tcrossprod(M_all))
    I_ext <- as.matrix(Matrix::tcrossprod(M_ext))
    I_bg  <- as.matrix(Matrix::tcrossprod(M_bg))

    # Sizes per rule
    n_all <- as.integer(diag(I_all))
    n_ext <- as.integer(diag(I_ext))
    n_bg  <- as.integer(diag(I_bg))

    # Jaccard matrices over extreme / background / all SNPs
    # Use the same convention as before: if union is 0, Jaccard = 0.
    .jaccard_from_intersections <- function(I, nA) {
      U <- outer(nA, nA, "+") - I
      J <- matrix(0, nrow = K, ncol = K)
      ok <- (U > 0)
      J[ok] <- I[ok] / U[ok]
      diag(J) <- 1
      J
    }

    labs <- paste0("r", seq_len(K))
    J_all <- .jaccard_from_intersections(I_all, n_all)
    J_ext <- .jaccard_from_intersections(I_ext, n_ext)
    J_bg  <- .jaccard_from_intersections(I_bg,  n_bg)

    dimnames(J_all) <- list(labs, labs)
    dimnames(J_ext) <- list(labs, labs)
    dimnames(J_bg)  <- list(labs, labs)

    # Also build a rich long-form overlap table:
    # for each (i,j) we store intersection sizes and directional proportions.
    n_pairs <- K * (K - 1L) / 2L

    i_index <- integer(n_pairs)
    j_index <- integer(n_pairs)

    # Fill (i,j) indices for i < j without allocating large row()/col() matrices
    pos <- 1L
    for (i in seq_len(K - 1L)) {
      nj <- K - i
      idx <- pos:(pos + nj - 1L)
      i_index[idx] <- i
      j_index[idx] <- (i + 1L):K
      pos <- pos + nj
    }

    # Extract intersections for i < j
    inter_all <- I_all[cbind(i_index, j_index)]
    keep <- (inter_all > 0)
    if (!any(keep)) {
      overlap_tbl <- data.table::data.table()
    } else {
      i_index <- i_index[keep]
      j_index <- j_index[keep]

      inter_all <- inter_all[keep]
      inter_ext <- I_ext[cbind(i_index, j_index)]
      inter_bg  <- I_bg [cbind(i_index, j_index)]

      n_all_i <- as.integer(n_all[i_index])
      n_all_j <- as.integer(n_all[j_index])

      n_ext_i <- as.integer(n_ext[i_index])
      n_ext_j <- as.integer(n_ext[j_index])

      n_bg_i <- as.integer(n_bg[i_index])
      n_bg_j <- as.integer(n_bg[j_index])

      n_all_intersect <- as.integer(inter_all)
      n_ext_intersect <- as.integer(inter_ext)
      n_bg_intersect  <- as.integer(inter_bg)

      n_all_unique_i  <- n_all_i - n_all_intersect
      n_all_unique_j  <- n_all_j - n_all_intersect

      prop_all_i_in_j <-
        ifelse(n_all_i > 0L, n_all_intersect / n_all_i, NA_real_)
      prop_all_j_in_i <-
        ifelse(n_all_j > 0L, n_all_intersect / n_all_j, NA_real_)

      # Jaccard vectors (duplicated from matrices, but handy here)
      jacc_all <- J_all[cbind(i_index, j_index)]
      jacc_ext <- J_ext[cbind(i_index, j_index)]
      jacc_bg  <- J_bg [cbind(i_index, j_index)]

      overlap_all <-
        n_all_intersect / pmax.int(1L, pmin.int(n_all_i, n_all_j))
      overlap_ext <- ifelse(
        pmin.int(n_ext_i, n_ext_j) > 0L,
        n_ext_intersect / pmin.int(n_ext_i, n_ext_j),
        NA_real_
      )
      overlap_bg  <- ifelse(pmin.int(n_bg_i,  n_bg_j)  > 0L,
                            n_bg_intersect / pmin.int(n_bg_i,  n_bg_j),
                            NA_real_)

      overlap_tbl <- data.table::data.table(
        i_index = i_index,
        j_index = j_index,

        # bucket sizes (all SNPs)
        n_all_i = n_all_i,
        n_all_j = n_all_j,
        n_all_intersect = n_all_intersect,
        n_all_unique_i  = n_all_unique_i,
        n_all_unique_j  = n_all_unique_j,

        # extremes & background intersection sizes
        n_ext_i = n_ext_i,
        n_ext_j = n_ext_j,
        n_bg_i  = n_bg_i,
        n_bg_j  = n_bg_j,

        n_bg_intersect  = n_bg_intersect,
        n_ext_intersect = n_ext_intersect,

        # directional overlap for all-SNP buckets
        prop_all_i_in_j = prop_all_i_in_j,
        prop_all_j_in_i = prop_all_j_in_i,

        # Jaccard (duplicated from matrices, but handy here)
        jacc_ext = jacc_ext,
        jacc_bg  = jacc_bg,
        jacc_all = jacc_all,

        overlap_all = overlap_all,
        overlap_ext = overlap_ext,
        overlap_bg  = overlap_bg
      )
    }
  }

  # Performance metrics
  precision <- n_ext / pmax(1e-12, (n_ext + n_bg))
  recall <- if (N_extr > 0L)
    n_ext / N_extr
  else
    rep(NA_real_, K)
  lift <-
    if (!is.na(base_rate) &&
        base_rate > 0)
      precision / base_rate
  else
    rep(NA_real_, K)

  # Optional beta–binomial smoothing (selection stability)
  precision_shrink <- rep(NA_real_, K)
  lift_shrink      <- rep(NA_real_, K)
  if (is.numeric(shrink_m) &&
      shrink_m > 0 && !is.na(base_rate) && base_rate > 0) {
    m <- as.numeric(shrink_m)
    precision_shrink <-
      (n_ext + m * base_rate) / (pmax(1e-12, (n_ext + n_bg)) + m)
    lift_shrink      <- precision_shrink / base_rate
  }

  # Build summary table: start with R_tbl (may exclude prefixes), then fill
  # missing rows
  summary_tbl <-
    R_tbl[, .(
      rule_str,
      rule_len,
      n_extreme,
      n_bg,
      support,
      support_all,
      precision,
      recall,
      lift,
      enrichment,
      med_y_extreme,
      med_y_bg,
      med_y_overall
    )]

  # Add rows for rules without summaries (e.g., prefixes) using fold-evaluated
  # counts
  have_sum <- summary_tbl$rule_str
  miss_sum <- setdiff(rule_ids, have_sum)

  if (length(miss_sum)) {
    add <- data.table::data.table(
      rule_str    = rule_ids,
      rule_len    = vapply(strsplit(rule_ids, " | ", fixed = TRUE), length, integer(1)),
      n_extreme   = as.integer(n_ext),
      n_bg        = as.integer(n_bg),
      support     = as.integer(n_ext + n_bg),
      support_all = as.integer(n_all),
      precision   = as.numeric(precision),
      recall      = as.numeric(recall),
      lift        = as.numeric(lift),
      enrichment  = as.numeric(NA_real_),
      med_y_extreme = as.numeric(NA_real_),
      med_y_bg      = as.numeric(NA_real_),
      med_y_overall = as.numeric(NA_real_)
    )[rule_str %chin% miss_sum]

    summary_tbl <-
      data.table::rbindlist(list(summary_tbl, add),
                            use.names = TRUE,
                            fill = TRUE)
  }

  # Attach smoothed columns (always present, may be NA)
  data.table::setkey(summary_tbl, rule_str)
  tmp_sm <- data.table::data.table(
    rule_str = rule_ids,
    precision_shrink = precision_shrink,
    lift_shrink      = lift_shrink
  )
  data.table::setkey(tmp_sm, rule_str)
  summary_tbl <- summary_tbl[tmp_sm, on = "rule_str"]
  data.table::setorderv(summary_tbl, c("lift", "precision"), c(-1L,-1L))

  summary_tbl <- .drop_all_na_cols(summary_tbl)
  overlap_tbl <- .drop_all_na_cols(overlap_tbl)

  list(
    summary     = summary_tbl[],
    overlap     = overlap_tbl[],
    rule_ids    = rule_ids,
    jaccard_ext = J_ext,
    jaccard_bg  = J_bg,
    jaccard_all = J_all,
    meta        = list(
      which       = which,
      fold_n_all  = length(n_tot),
      fold_n_extr = N_extr,
      fold_n_bg   = N_bg,
      base_rate   = base_rate,
      shrink_m    = shrink_m
    )
  )
}
