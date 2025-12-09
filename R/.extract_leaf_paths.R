#' Title
#'
#' @param native_leaf_ids
#' @param Tm
#' @param tdt
#' @param helpers
#' @param trees_per_batch
#'
#' @return
#' @keywords internal
#'
#' @examples
.extract_leaf_paths <- function(native_leaf_ids,
                                Tm,
                                tdt,
                                helpers,
                                trees_per_batch = 250L) {
  # tdt contains the per-tree node tables (Tree, ID, Feature, Split, Yes, No)
  # We key by (Tree, ID) so later lookups are cheap.
  data.table::setkey(tdt, Tree, ID)

  build_raw_steps_batch <- function(batch_start, batch_end) {
    # total_steps counts the total number of path steps across all leaves
    # in this batch of trees. We use it to preallocate vectors.
    total_steps <- 0L

    # meta will store basic metadata per tree in this batch:
    # tt = 1-based tree index
    # tr = 0-based tree index (used by xgboost helpers)
    # uleaf = vector of native leaf_ids for this tree
    meta <- vector("list", batch_end - batch_start + 1L)
    mi <- 1L

    # First pass over trees in the batch:
    # - determine which leaves exist in each tree
    # - count total path length across all leaves (sum of number of nodes
    # along each leaf's path), so we know how big to make the output arrays.
    for (tt in batch_start:batch_end) {
      tr <- tt - 1L
      uleaf <- as.integer(native_leaf_ids[[tt]])
      meta[[mi]] <- list(tt = tt,
                         tr = tr,
                         uleaf = uleaf)

      # For each leaf in this tree, ask the helper for the sequence of
      # internal node IDs that lead to that leaf; each node contributes
      # one "step" in the final path.
      for (leaf_id in uleaf) {
        # path node IDs leading to leaf
        pn <- helpers$get_leaf_path(tr, leaf_id)$nodes
        if (length(pn)) {
          total_steps <- total_steps + length(pn)
        }
      }
      mi <- mi + 1L
    }

    # If there are no path steps at all in this batch, return an empty
    # table with the expected schema.
    if (total_steps == 0L) {
      return(
        data.table::data.table(
          Tree = integer(0),
          leaf_id = integer(0),
          depth = integer(0),
          feature = character(0),
          direction = character(0),
          split_val = numeric(0)
        )
      )
    }

    # Preallocate flat vectors for all path steps in this batch.
    # These will hold one row per (tree, leaf, depth) combination.
    Tree_v  <- integer(total_steps)
    Leaf_v  <- integer(total_steps)
    Depth_v <- integer(total_steps)
    Feat_v  <- character(total_steps)
    Dir_v   <- character(total_steps)
    Split_v <- numeric(total_steps)
    cursor <- 0L

    # Second pass over trees in the batch:
    # - build fast lookup arrays (yesA, noA, featA, spltA) from tdt for
    # this tree's node IDs
    # - for each leaf, expand its path into per-step (feature, direction,
    # split) entries and write into the preallocated vectors.
    for (m in meta) {
      tt <- m$tt
      tr <- m$tr
      uleaf <- m$uleaf

      # Subset node table for this tree (node IDs, children, feature, split)
      tt_dt <- tdt[.(tr)]
      if (!nrow(tt_dt))
        next

      # Build compact lookup arrays indexed by node ID + 1.
      # This lets us map node IDs -> (Yes, No, Feature, Split) in O(1).
      maxid <- max(tt_dt$ID, na.rm = TRUE)

      yesA <- rep.int(NA_integer_, maxid + 1L)
      yesA[tt_dt$ID + 1L] <- tt_dt$Yes

      noA <- rep.int(NA_integer_, maxid + 1L)
      noA [tt_dt$ID + 1L] <- tt_dt$No

      featA <- rep.int(NA_character_, maxid + 1L)
      featA[tt_dt$ID + 1L] <- as.character(tt_dt$Feature)

      spltA <- rep.int(NA_real_, maxid + 1L)
      spltA[tt_dt$ID + 1L] <-
        suppressWarnings(as.numeric(tt_dt$Split))

      # Loop over all leaves for this tree, and expand each leaf path
      # into its sequence of decision steps.
      for (leaf_id in uleaf) {
        # pn = internal node IDs along the path from root to this leaf
        pn <- helpers$get_leaf_path(tr, leaf_id)$nodes
        k <- length(pn)
        if (!k)
          next

        # child_along[i] = which child node we take after node pn[i]
        # along the path. For the last internal node, the "child" is
        # the leaf itself.
        child_along <- integer(k)
        if (k > 1L)
          child_along[1:(k - 1L)] <- pn[2:k]
        child_along[k] <- leaf_id

        # Use the node IDs to look up Yes/No children, feature, and split.
        parents <- pn + 1L
        yes     <- yesA [parents]
        no      <- noA [parents]
        feat    <- featA[parents]
        splt    <- spltA[parents]

        # Determine whether the path took the "yes" branch ("<") or
        # the "no" branch (">=") at each step. If neither matches,
        # label as "missing".
        dir <- ifelse(yes == child_along,
                      "<",
                      ifelse(no == child_along, ">=", "missing"))

        # Write these k steps into the preallocated vectors.
        # depth is 0,1,...,k-1 along the path for this leaf.
        rng <- (cursor + 1L):(cursor + k)

        Tree_v [rng] <- tr
        Leaf_v [rng] <- leaf_id
        Depth_v[rng] <- 0L:(k - 1L)
        Feat_v [rng] <- feat
        Dir_v  [rng] <- dir
        Split_v[rng] <- splt

        cursor <- cursor + k
      }
    }

    # Wrap batch-level vectors into a data.table.
    data.table::data.table(
      Tree      = Tree_v,
      leaf_id   = Leaf_v,
      depth     = Depth_v,
      feature   = Feat_v,
      direction = Dir_v,
      split_val = Split_v
    )
  }

  # Loop over all trees in batches to avoid large one-shot allocations.
  out_list <- vector("list", ceiling(Tm / trees_per_batch))
  oi <- 1L
  for (batch_start in seq(1L, Tm, by = trees_per_batch)) {
    batch_end      <- min(batch_start + trees_per_batch - 1L, Tm)
    out_list[[oi]] <- build_raw_steps_batch(batch_start, batch_end)
    oi <- oi + 1L
  }

  # Bind all batches into a single leaf-path table.
  LS <-
    data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(LS)) {
    data.table::setkey(LS, Tree, leaf_id)
    # Keep schema consistent with non-empty case
    return(LS[, .(Tree, leaf_id, depth, feature, direction, split_val = numeric())])
  }

  # Drop steps where feature is NA or empty and normalize types.
  LS <- LS[!is.na(feature) & nzchar(feature)]
  LS[, `:=`(
    Tree      = as.integer(Tree),
    leaf_id   = as.integer(leaf_id),
    depth     = as.integer(depth),
    feature   = as.character(feature),
    direction = as.character(direction),
    split_val = as.numeric(split_val)
  )]

  # Key by (Tree, leaf_id) for fast joins downstream.
  data.table::setkey(LS, Tree, leaf_id)
  LS[]
}
