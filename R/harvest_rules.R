#' Title
#'
#' @param boosted
#' @param max_depth
#' @param min_support
#' @param trees_subset
#' @param trees_per_batch
#' @param progress_every
#' @param tighten_monotone
#' @param return_ledger
#'
#' @return
#' @export
#'
#' @examples

harvest_rules <- function(boosted,
                          max_depth        = 5L,
                          min_support      = 20L,
                          trees_subset     = NULL,
                          trees_per_batch  = 250L,
                          progress_every   = NULL,
                          tighten_monotone = TRUE,
                          return_ledger    = TRUE) {
  # Signature & basic checks
  FUN <- "harvest_rules"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted")) {
    stop(
      sprintf(
        "[%s] Input must be an object of class 'boosted' (from make_boosted()).",
        FUN
      )
    )
  }
  if (!inherits(boosted, "boosted_binned")) {
    stop(sprintf(
      "[%s] boosted is not harvest-ready; run prepare_harvest() first.",
      FUN
    ))
  }

  # Pull train data from boosted
  Tm             <- boosted$Tm
  extr_idx_train <- boosted$extr_idx_train
  bg_idx_train   <- boosted$bg_idx_train
  N_extr_train   <- boosted$N_extr_train
  N_bg_train     <- boosted$N_bg_train
  n_yvar_train   <- boosted$n_yvar_train

  train_leaf_map <- boosted$train_leaf_map
  dense_leaf_ids <- train_leaf_map$dense_leaf_ids
  native_ids_all <- train_leaf_map$native_leaf_ids
  n_leaves       <- as.integer(train_leaf_map$n_leaves)

  # Feature bins
  bin_spec <- boosted$harvest_bins

  # Base rate & optional yvar (for medians)
  base_rate <- boosted$base_rate_train
  y_num     <- boosted$yvar_train

  # -----------------------------------
  # Step 1: turn leaf-level paths + bins into canonical rule strings
  # and a rule-length map.
  # -----------------------------------

  # Build rule strings to max_depth using precomputed leaf paths
  leaf_paths <- boosted$leaf_paths
  PATHS_full <- .build_rule_strings(
    leaf_paths,
    bin_spec         = bin_spec,
    max_depth        = max_depth,
    tighten_monotone = tighten_monotone
  )

  # Extract the Tree, leaf_id, and rule_str columns and key by (Tree, leaf_id)
  # This is the map we later use to attach a rule string to each leaf.
  PATHS <- PATHS_full[, .(Tree, leaf_id, rule_str)]
  data.table::setkey(PATHS, Tree, leaf_id)
  rm(PATHS_full); gc()

  # Optional: restrict trees to a supplied subset (0-based indices)
  all_trees0 <- 0:(Tm - 1L)
  if (!is.null(trees_subset)) {
    trees_subset <- as.integer(sort(unique(trees_subset)))
    trees_subset <- trees_subset[trees_subset %in% all_trees0]
    if (!length(trees_subset)) {
      stop(sprintf("[%s] trees_subset is empty after filtering.", FUN))
    }
    use_tt <- trees_subset + 1L     # convert 0-based to 1-based indices
  } else {
    use_tt <- seq_len(Tm)
  }

  # -----------------------------------
  # Step 2: find leaves with labeled support >= min_support
  # and attach rule strings
  # -----------------------------------

  # Initialize global containers and counters
  # - pairs_by_batch: (Tree, leaf_id, rule_str) triples from each batch
  pairs_by_batch <- vector("list", ceiling(length(use_tt) / trees_per_batch))
  pbb <- 1L

  # For each batch of trees:
  # 1) compute labeled support per leaf (extreme/background)
  # 2) filter leaves by min_support
  # 3) attach rule_str via PATHS
  # 4) accumulate (Tree, leaf_id, rule_str) pairs
  for (batch_start in seq(1L, length(use_tt), by = trees_per_batch)) {
    # Find where to stop for this batch
    batch_end   <- min(batch_start + trees_per_batch - 1L, length(use_tt))

    # Initialize container and counter
    # - batch_pairs: (Tree, leaf_id, rule_str) records for this batch
    batch_pairs <- vector("list", batch_end - batch_start + 1L)
    bp <- 1L

    # Loop over trees in this batch
    for (tt in use_tt[batch_start:batch_end]) {
      tr    <- tt - 1L              # 0-based tree index
      inv_t <- dense_leaf_ids[[tt]] # leaf assignment per SNP
      nb    <- n_leaves[tt]         # number of leaves in this tree

      # Count extreme/background SNPs per leaf
      ce_dense <- tabulate(inv_t[extr_idx_train], nbins = nb)
      cb_dense <- tabulate(inv_t[bg_idx_train],   nbins = nb)

      # Filter out leaves with labeled support < min_support
      support_dense <- ce_dense + cb_dense
      keep_dense    <- which(support_dense >= as.integer(min_support))

      # If no leaves pass, skip this tree
      if (!length(keep_dense)) {
        bp <- bp + 1L
        next
      }

      # Map dense IDs back to native IDs for compatibility with rule strings
      native_ids_t <- native_ids_all[[tt]]
      if (length(native_ids_t) < nb) {
        stop(
          sprintf(
            "[%s] internal: native_leaf_ids[[%d]] has length %d < n_leaves = %d.",
            FUN,
            tt,
            length(native_ids_t),
            nb
          )
        )
      }
      leaf_ids_native <- native_ids_t[keep_dense]

      # Per-leaf labeled counts for this tree (extreme/background)
      DT0 <- data.table::data.table(
        Tree      = tr,
        leaf_id   = leaf_ids_native,
        n_extreme = ce_dense[keep_dense],
        n_bg      = cb_dense[keep_dense]
      )

      # Attach rule string for each (Tree, leaf_id); drop unmapped leaves
      J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
      if (!nrow(J)) {
        bp <- bp + 1L
        next
      }

      # Store which leaves fall under each rule
      batch_pairs [[bp]]    <- J[, .(Tree, leaf_id, rule_str)]
      bp <- bp + 1L
    }

    # Add results from this batch to global container
    pairs_by_batch[[pbb]] <- data.table::rbindlist(batch_pairs,
                                                   use.names = TRUE,
                                                   fill      = TRUE)
    # Increment counter
    pbb <- pbb + 1L

    # Tree-level progress if requested
    if (!is.null(progress_every) && progress_every > 0L &&
        ((batch_end %% progress_every) == 0L ||
         batch_end == length(use_tt))) {
      message(
        sprintf(
          "[%s] processed tree %d / %d (%.1f%%)",
          FUN,
          use_tt[batch_end],
          length(use_tt),
          100 * batch_end / length(use_tt)
        )
      )
    }
  }

  # Collate mapping from all batches:
  # each row says "in Tree t, leaf_id L corresponds to rule_str r".
  pairs_all <- data.table::rbindlist(pairs_by_batch,
                                     use.names = TRUE,
                                     fill = TRUE)
  if (!nrow(pairs_all)) {
    stop(sprintf(
      "[%s] no rule (Tree, leaf_id) pairs to pool after filtering.",
      FUN
    ))
  }

  # -----------------------------------
  # Step 3: pool SNPs under each unique rule and evaluate enrichment
  # -----------------------------------

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

  # Identify unique rule strings across trees
  uniq_rules <- sort(unique(pairs_all$rule_str))
  # Pre-allocate container for pooled rule-level summaries
  pooled     <- vector("list", length(uniq_rules))

  # Pool SNPs under each rule and compute enrichment metrics
  for (i in seq_along(uniq_rules)) {
    rs <- uniq_rules[i]
    # All (Tree, leaf_id) pairs where this rule held in the ensemble
    pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
    if (!nrow(pairs))
      next
    # Look up SNP indices under this rule via the per-tree lookup
    buckets <- .snp_lookup(
      pairs            = pairs,
      snps_ext_by_leaf = snps_ext_by_leaf,
      snps_bg_by_leaf  = snps_bg_by_leaf,
      snps_all_by_leaf = snps_all_by_leaf
    )
    bucket_ext <- buckets$bucket_ext
    bucket_bg  <- buckets$bucket_bg
    bucket_all <- buckets$bucket_all

    # Compute how many extreme vs background SNPs the rule captures
    n_e         <- length(bucket_ext)
    n_b         <- length(bucket_bg)
    support     <- n_e + n_b          # labeled support
    support_all <- length(bucket_all) # total support (labeled + unlabeled)

    # Performance metrics: how well the rule enriches for extreme SNPs
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

    # Haldane–Anscombe corrected odds ratio (finite with zeros)
    or_ha <- ((n_e + 0.5) * (N_bg_train - n_b + 0.5)) /
      ((N_extr_train - n_e + 0.5) * (n_b + 0.5))

    # Rule-level LLRs with Jeffreys smoothing on contingency counts
    enrichment <- .rule_llr(
      n_extreme    = n_e,
      n_bg         = n_b,
      N_extr_total = N_extr_train,
      N_bg_total   = N_bg_train
    )

    # Optional: include median of yvar in output tables
    med_e <- if (!is.null(y_num) && n_e) {
      median(y_num[bucket_ext], na.rm = TRUE)
    } else
      NA_real_

    med_b <- if (!is.null(y_num) && n_b) {
      median(y_num[bucket_bg], na.rm = TRUE)
    } else
      NA_real_

    med_o <- if (!is.null(y_num) && support_all > 0L) {
      median(y_num[bucket_all], na.rm = TRUE)
    } else
      NA_real_

    # Compute rule length from rule string
    rl <- length(strsplit(rs, " \\| ", fixed = FALSE)[[1]])

    # Save pooled result for this rule
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
      med_y_extreme = med_e,
      med_y_bg      = med_b,
      med_y_overall = med_o
    )

    # Optional per-rule progress
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

  # -----------------------------------
  # Step 4: prepare and return the results
  # -----------------------------------

  R_tbl <- data.table::rbindlist(pooled, use.names = TRUE, fill = TRUE)
  if (!nrow(R_tbl)) {
    warning(sprintf("[%s] no pooled rules.", FUN))
    return(invisible(NULL))
  }

  # Order results
  data.table::setorderv(
    R_tbl,
    cols  = c("lift", "odds_ratio", "precision", "recall"),
    order = c(-1L,-1L,-1L,-1L)
  )

  if (isTRUE(return_ledger)) {
    out <- list(
      R            = R_tbl[],
      pairs_all    = pairs_all[],
      meta         = list(max_depth        = max_depth,
                          tighten_monotone = tighten_monotone,
                          n_train          = n_yvar_train,
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
