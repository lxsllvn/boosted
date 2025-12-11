#' Title
#'
#' @param tdt
#'
#' @return
#' @keywords internal
#'
#' @examples

.build_tree_helpers <- function(tdt) {

  tdt <- data.table::as.data.table(tdt)
  trees <- sort(unique(tdt$Tree))

  # 1. Build child -> parent lookup for each tree
  #
  # For each tree ID, we create an integer vector `pmap` such that:
  # pmap[child_id + 1] = parent_id
  # If a node has no parent (root), the entry stays NA.
  #
  # This lets us walk from a leaf back to the root using only array
  # lookups, which is much cheaper than repeated joins on `tdt`.

  parent_env <- new.env(parent = emptyenv())
  for (tr in trees) {
    tt     <- tdt[Tree == tr, .(ID, Yes, No, Missing)]
    max_id <- max(tt$ID, 0L)

    # child_id is 0-based, so allocate max_id + 1 and index at child_id + 1
    pmap <- rep(NA_integer_, max_id + 1L)

    # For each child type (Yes / No / Missing), record its parent ID
    for (col in c("Yes", "No", "Missing")) {
      ch <- tt[[col]]
      ok <- !is.na(ch)
      if (any(ok)) {
        # Example: if child = 5 and parent = 2, then pmap[5 + 1] = 2
        pmap[ch[ok] + 1L] <- tt$ID[ok]
      }
    }

    # Store the parent map keyed by tree ID (as character)
    parent_env[[as.character(tr)]] <- pmap
  }


  # 2. Precompute split-level gain and cover (for non-leaf nodes)
  #
  # We keep only rows where Leaf == FALSE, so each ID here represents
  # an internal split node. These will be used to annotate paths with
  # Gain and SplitCover when we reconstruct leaf-to-root paths.

  split_gc <-
    tdt[Leaf == FALSE, .(Tree, ID, Gain, SplitCover = Cover)]
  data.table::setkey(split_gc, Tree, ID)


  # 3. Caches for per-tree depth maps, max depth, and leaf paths
  #
  # depth_env : Tree -> integer vector of node depths (index by ID + 1)
  # maxd_env  : Tree -> single integer (maximum node depth for that tree)
  # leaf_cache: Tree -> (environment of leaf_id -> path info)
  #
  # These are filled on demand and then reused.

  depth_env  <- new.env(parent = emptyenv())
  maxd_env   <- new.env(parent = emptyenv())
  leaf_cache <- new.env(parent = emptyenv())


  # 4. get_depth_map(tr)
  #
  # For a given tree `tr`, returns an integer vector `dep` where
  # dep[node_id + 1] = depth of that node (root = 0).
  #
  # Depths are computed with a simple BFS over the parent map.
  # Results are cached in `depth_env` for reuse.

  get_depth_map <- function(tr) {
    key <- as.character(tr)

    # Return cached depth map if it already exists
    dep <- depth_env[[key]]
    if (!is.null(dep)) {
      return(dep)
    }

    # Look up the parent map for this tree
    pmap <- parent_env[[key]]
    if (is.null(pmap)) {
      return(NULL)
    }

    max_id <- length(pmap) - 1L
    dep    <- rep(NA_integer_, max_id + 1L)

    # Root node is assumed to be ID = 0 with depth 0
    dep[1] <- 0L

    # BFS queue over node IDs (0-based)
    q    <- integer(max_id + 1L)
    head <- 1L
    tail <- 1L
    q[1] <- 0L

    while (head <= tail) {
      id <- q[head]
      head <- head + 1L

      # Children are all indices where parent_map == id
      ch_idx <- which(pmap == id)
      if (length(ch_idx)) {
        new    <- ch_idx - 1L # convert back to 0-based IDs
        unseen <- is.na(dep[new + 1L])
        if (any(unseen)) {
          # Depth of child = depth(parent) + 1
          dep[new[unseen] + 1L] <- dep[id + 1L] + 1L

          # Push unseen children onto the queue
          n_new <- sum(unseen)
          q[tail + seq_len(n_new)] <- new[unseen]
          tail <- tail + n_new
        }
      }
    }

    # Cache and return
    depth_env[[key]] <- dep
    maxd_env[[key]]  <- max(dep, na.rm = TRUE)
    dep
  }


  # 5. get_max_depth(tr)
  #
  # Convenience wrapper: returns the maximum depth of any node in tree `tr`.
  # Uses the cached value if available, otherwise forces computation via
  # get_depth_map().

  get_max_depth <- function(tr) {
    key <- as.character(tr)

    # Return cached max depth if present
    md <- maxd_env[[key]]
    if (!is.null(md)) {
      return(md)
    }

    # Otherwise compute depth map, then take maximum
    dep <- get_depth_map(tr)
    if (is.null(dep)) {
      return(0L)
    }
    md <- max(dep, na.rm = TRUE)
    maxd_env[[key]] <- md
    md
  }


  # 6. get_leaf_path(tr, leaf_id)
  #
  # For a given tree `tr` and leaf node ID `leaf_id`, returns: list(nodes =
  # integer vector of internal node IDs along the path from root to the leaf (in
  # top-down order), gain = numeric vector of Gain values for those nodes,
  # gcover = numeric vector of SplitCover (Cover) values )
  #
  # The path is defined by walking from the leaf up to the root using
  # the child -> parent map, then reversing the order to get root->leaf.
  #
  # Results are cached per (tree, leaf_id) so repeated calls are cheap.

  get_leaf_path <- function(tr, leaf_id) {
    tkey <- as.character(tr)

    # Ensure there is a per-tree cache environment
    env <- leaf_cache[[tkey]]
    if (is.null(env)) {
      env <- new.env(parent = emptyenv())
      leaf_cache[[tkey]] <- env
    }

    # Check if we've already computed this leaf’s path
    hit <- env[[as.character(leaf_id)]]
    if (!is.null(hit)) {
      return(hit)
    }

    # If no parent map for this tree, return empty path
    pmap <- parent_env[[tkey]]
    if (is.null(pmap)) {
      out <- list(nodes  = integer(0),
                  gain   = numeric(0),
                  gcover = numeric(0))
      leaf_cache[[tkey]][[as.character(leaf_id)]] <- out
      return(out)
    }


    # 6a. Walk from leaf up to root to determine path length `k`

    cur <- leaf_id
    k <- 0L
    while (!is.na(cur) && cur != 0L) {
      par <- pmap[cur + 1L]
      if (is.na(par)) {
        break
      }
      k <- k + 1L
      cur <- par
    }

    # If the leaf has no valid parent chain, return empty path
    if (k == 0L) {
      out <- list(nodes  = integer(0),
                  gain   = numeric(0),
                  gcover = numeric(0))
      leaf_cache[[tkey]][[as.character(leaf_id)]] <- out
      return(out)
    }


    # 6b. Reconstruct the full sequence of internal node IDs on the path
    #
    # We walk leaf -> root again, fill `nodes` in reverse, so that the
    # final vector is ordered from root to leaf.

    nodes <- integer(k)
    cur <- leaf_id
    for (pos in k:1) {
      par <- pmap[cur + 1L]
      nodes[pos] <- par
      cur <- par
      if (is.na(cur) || cur == 0L) {
        break
      }
    }

    # Join in Gain and SplitCover for these internal nodes
    sgc <- split_gc[.(tr, nodes), .(Gain, SplitCover)]

    out <- list(
      nodes  = nodes,
      gain   = as.numeric(sgc$Gain),
      gcover = as.numeric(sgc$SplitCover)
    )

    # Cache the result for (tree, leaf_id)
    leaf_cache[[tkey]][[as.character(leaf_id)]] <- out
    out
  }

  # Return a small list of helper functions used elsewhere
  list(
    get_depth_map = get_depth_map,
    get_max_depth = get_max_depth,
    get_leaf_path = get_leaf_path
  )
}
