#' Find maximum realized tree depth in an xgboost model.
#'
#' @param tdt
#'
#' @return
#' @export
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
