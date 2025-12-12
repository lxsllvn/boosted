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
#' @import data.table
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
    # start with a modest capacity and grow as needed
    capacity <- 100000L

    Tree_v  <- integer(capacity)
    Leaf_v  <- integer(capacity)
    Depth_v <- integer(capacity)
    Feat_v  <- character(capacity)
    Dir_v   <- character(capacity)
    Split_v <- numeric(capacity)
    cursor  <- 0L

    # helper: ensure we have room for k more steps
    ensure_capacity <- function(k) {
      needed <- cursor + k
      if (needed <= capacity)
        return()

      new_cap <- max(capacity * 2L, needed)

      length(Tree_v)  <<- new_cap
      length(Leaf_v)  <<- new_cap
      length(Depth_v) <<- new_cap
      length(Feat_v)  <<- new_cap
      length(Dir_v)   <<- new_cap
      length(Split_v) <<- new_cap

      capacity <<- new_cap
    }

    # single pass over trees in the batch
    for (tt in batch_start:batch_end) {
      tr    <- tt - 1L
      uleaf <- native_leaf_ids[[tt]]
      uleaf <- as.integer(uleaf)

      # Subset node table for this tree (node IDs, children, feature, split)
      tt_dt <- tdt[.(tr)]

      # Build compact lookup arrays indexed by node ID + 1.
      maxid <- max(tt_dt$ID, na.rm = TRUE)
      maxid <- as.integer(maxid)

      yesA <- rep.int(NA_integer_, maxid + 1L)
      yesA[tt_dt$ID + 1L] <- tt_dt$Yes

      noA <- rep.int(NA_integer_, maxid + 1L)
      noA[tt_dt$ID + 1L] <- tt_dt$No

      featA <- rep.int(NA_character_, maxid + 1L)
      featA[tt_dt$ID + 1L] <- as.character(tt_dt$Feature)

      spltA <- rep.int(NA_real_, maxid + 1L)
      spltA[tt_dt$ID + 1L] <-
        suppressWarnings(as.numeric(tt_dt$Split))

      # Parent map + max depth are guaranteed to exist for trees
      # coming from make_boosted() / .build_tree_helpers().
      pmap <- helpers$get_parent_map(tr)
      max_depth <- as.integer(helpers$get_max_depth(tr))

      # buffer to hold a leaf→root chain (max length = max_depth)
      tmp_nodes <- integer(max_depth)

      # Loop over all leaves for this tree, and expand each leaf path
      # into its sequence of decision steps.
      for (leaf_id in uleaf) {
        # leaf → root walk using parent map
        k   <- 0L
        cur <- leaf_id

        while (!is.na(cur) && cur != 0L) {
          par <- pmap[cur + 1L]
          if (is.na(par)) break

          if (k < max_depth) {
            k <- k + 1L
            tmp_nodes[k] <- par
          } else {
            # if something goes weird with max_depth, bail on this leaf
            k <- 0L
            break
          }
          cur <- par
        }

        # If this leaf has no valid internal path, skip it
        if (k == 0L)
          next

        # nodes from root → leaf = reversed leaf→root chain
        pn <- tmp_nodes[seq.int(k, 1L)]

        ensure_capacity(k)

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

    used <- seq_len(cursor)

    data.table::data.table(
      Tree      = Tree_v [used],
      leaf_id   = Leaf_v [used],
      depth     = Depth_v[used],
      feature   = Feat_v [used],
      direction = Dir_v  [used],
      split_val = Split_v[used]
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
  LS <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(LS)) {
    stop("[.extract_leaf_paths] No leaf paths were extracted. Check xgboost model and .parse_xgboost_tree() output.")
  }

  LS <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)

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
