#' Find the maximum realized split depth across all trees in an `xgboost` model
#'
#' Performs a breadth-first traversal of every tree in \code{tdt}, tracking
#' the number of split-steps from the root to each node, and returns the
#' maximum depth observed across all leaves in all trees. Used by
#' \code{\link{make_boosted}} to set \code{boosted$max_depth}, which serves as
#' the default prefix depth cap in \code{\link{prepare_harvest}}.
#'
#' @param tdt A \code{data.table} returned by \code{\link{.parse_xgb_tree}},
#'   with columns \code{Tree}, \code{ID}, \code{Yes}, \code{No}, and
#'   \code{Leaf}.
#'
#' @return Integer scalar: the maximum number of split-steps from root to leaf
#'   observed across the entire ensemble.
#' @keywords internal
#' @examples
.infer_max_depth <- function(tdt) {
  trees <- sort(unique(as.integer(tdt$Tree)))
  max_depth <- 0L

  for (tr in trees) {
    dt <- tdt[Tree == tr, .(
      ID   = as.integer(ID),
      Yes  = as.integer(Yes),
      No   = as.integer(No),
      Leaf = as.logical(Leaf)
    )]
    if (!nrow(dt))
      next

    ids <- dt$ID
    max_id <- max(ids, na.rm = TRUE)

    row_of <- rep.int(NA_integer_, max_id + 1L)
    row_of[ids + 1L] <- seq_len(nrow(dt))

    root <- 0L
    if (is.na(row_of[root + 1L]))
      root <- ids[[1L]]

    depth <- rep.int(NA_integer_, max_id + 1L)
    depth[root + 1L] <- 0L

    q <- integer(nrow(dt))
    head <- 1L
    tail <- 1L
    q[tail] <- root

    while (head <= tail) {
      nid <- q[head]
      head <- head + 1L
      i <- row_of[nid + 1L]
      if (is.na(i))
        next

      d <- depth[nid + 1L]

      # Only expand internal nodes
      if (!isTRUE(dt$Leaf[[i]])) {
        for (child in c(dt$Yes[[i]], dt$No[[i]])) {
          if (is.na(child))
            next
          if (child < 0L || (child + 1L) > length(depth))
            next
          if (is.na(depth[child + 1L])) {
            depth[child + 1L] <- d + 1L
            tail <- tail + 1L
            if (tail > length(q))
              q <- c(q, integer(length(q)))
            q[tail] <- child
          }
        }
      }
    }

    max_depth <- max(max_depth, max(depth, na.rm = TRUE))
  }
  as.integer(max_depth)
}
