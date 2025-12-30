#' Title
#'
#' @param boosted
#' @param max_depth
#' @param min_support
#' @param trees_subset
#' @param trees_per_batch
#' @param progress_every
#' @param tighten_monotone
#'
#' @return
#' @export
#'
#' @examples
harvest_rules <- function(boosted,
                          min_support      = 20L,
                          trees_per_batch  = 250L,
                          max_depth        = NULL,
                          progress_every   = NULL,
                          trees_subset     = NULL,
                          tighten_monotone = TRUE) {
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

  # If max_depth is provided, check if it's an integer <= boosted$max_depth
  model_md <- as.integer(boosted$max_depth)

  if (is.null(max_depth)) {
    max_depth <- model_md
  } else {
    bad <- length(max_depth) != 1L ||
      !is.numeric(max_depth) ||
      is.na(max_depth) ||
      !is.finite(max_depth) ||
      max_depth != as.integer(max_depth) ||
      max_depth < 1L ||
      max_depth > model_md

    if (bad) {
      stop(sprintf(
        "[%s] max_depth must be a single finite integer in [1, %d].",
        FUN, model_md
      ))
    }
    max_depth <- as.integer(max_depth)
  }

  # Pull train data from boosted
  Tm             <- boosted$Tm
  extr_idx_train <- boosted$extr_idx_train
  bg_idx_train   <- boosted$bg_idx_train
  n_yvar_train   <- boosted$n_yvar_train

  # Pull native and dense leaf IDs
  dense_leaf_ids <- boosted$train_leaf_map$dense_leaf_ids
  native_ids_all <- boosted$train_leaf_map$native_leaf_ids
  n_leaves       <- as.integer(boosted$train_leaf_map$n_leaves)

  # Feature bins
  bin_spec <- boosted$harvest_bins

  # Minimum labeled support
  min_support <- as.integer(min_support)

  # -----------------------------------
  # Step 1: turn leaf-level paths + bins into canonical rule strings
  # and a rule-length map.
  # -----------------------------------

  # Build per-leaf rule definitions (rule_str, rule_len) for labeling leaves
  leaf_paths <- boosted$leaf_paths

  PATHS <- .build_rule_strings(
    leaf_paths,
    bin_spec         = bin_spec,
    max_depth        = max_depth,
    tighten_monotone = tighten_monotone
  )

  data.table::setkey(PATHS, Tree, leaf_id)


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
    # - batch_pairs: (Tree, leaf_id, rule_str, rule_len) records for this batch
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
      keep_dense    <- which(support_dense >= min_support)

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
      J <-
        PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L][, .(Tree, leaf_id, rule_str, rule_len)]

      if (!nrow(J)) {
        bp <- bp + 1L
        next
      }

      # Store which leaves fall under each rule
      batch_pairs [[bp]] <- J[, .(Tree, leaf_id, rule_str, rule_len)]
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
        sprintf("[%s] processed tree %d / %d", FUN, batch_end, length(use_tt))
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
      "[%s] no rule (Tree, leaf_id) after filtering.",
      FUN
    ))
  }

  # -----------------------------------
  # Step 3: return results
  # -----------------------------------
  out <- list(
    pairs_all    = pairs_all[],
    meta         = list(
      max_depth        = max_depth,
      tighten_monotone = tighten_monotone,
      n_train          = n_yvar_train,
      Tm               = Tm
    )
  )

  # Tag as a rule-harvest ledger without interfering with other list-based tools
  class(out) <- c("boosted_harvest", class(out))
  return(out)

}


