#' Build a rule cache from a parsed and binned `xgboost` tree table
#'
#' Traverses every tree in \code{tdt} depth-first, recording the sequence of
#' split conditions encountered on the path from root to each leaf. For every
#' leaf, all prefix rules of length \eqn{1, \ldots, \min(\text{max\_depth},
#' \text{leaf\_depth})}{} are materialized as character strings and stored
#' alongside the leaf's \code{(Tree, leaf_id)} identity. Split thresholds are
#' taken from the \code{Split_bin} column (the binned values added by
#' \code{\link{prepare_harvest}}) rather than the raw \code{Split} values,
#' so that rules are expressed in terms of stable, interpretable bin
#' midpoints. When \code{tighten_monotone = TRUE}, multiple constraints on the
#' same feature and direction within a single prefix are collapsed to the
#' tightest bound, removing logically redundant clauses.
#'
#' @param tdt A \code{data.table} produced by \code{\link{.parse_xgb_tree}}
#'   and subsequently annotated by \code{\link{prepare_harvest}} with a
#'   \code{Split_bin} column. Required columns: \code{Tree}, \code{ID},
#'   \code{Feature}, \code{Split_bin}, \code{Yes}, \code{No}, and \code{Leaf}.
#' @param max_depth Positive integer. Maximum number of split-steps to include
#'   in any rule prefix. Rules for leaves deeper than \code{max_depth} are
#'   truncated at \code{max_depth} splits; \code{leaf_depth} is still recorded
#'   faithfully for all leaves regardless of the cap.
#' @param tighten_monotone Logical. If \code{TRUE} (default), redundant
#'   monotone constraints on the same feature and direction within a prefix are
#'   collapsed to the tightest bound. For example, \code{x < 5 | x < 3}
#'   becomes \code{x < 3}.
#'
#' @return A \code{data.table} with one row per \code{(Tree, leaf_id,
#'   rule_len)} combination, indexed (not keyed) on \code{rule_str},
#'   \code{(Tree, leaf_id)}, and \code{(rule_len, rule_str)} for fast joins.
#'   Columns:
#' \describe{
#'   \item{\code{Tree}}{Integer. Zero-based tree index, matching
#'     \code{boosted$tdt$Tree}.}
#'   \item{\code{leaf_id}}{Integer. Native `xgboost` node ID of the leaf,
#'     matching \code{boosted$tdt$ID} for leaf rows.}
#'   \item{\code{leaf_depth}}{Integer. Total number of split-steps from the
#'     root to this leaf (regardless of \code{max_depth}).}
#'   \item{\code{rule_len}}{Integer. Number of split-steps included in this
#'     rule prefix (\eqn{1 \le \text{rule\_len} \le \min(\text{max\_depth},
#'     \text{leaf\_depth})}{}).}
#'   \item{\code{n_clauses}}{Integer. Number of distinct split conditions in
#'     \code{rule_str} after monotone tightening (may be less than
#'     \code{rule_len} when tightening collapses redundant conditions).}
#'   \item{\code{rule_str}}{Character. The rule string: an ordered,
#'     pipe-delimited sequence of split conditions of the form
#'     \code{"feature >= threshold"} or \code{"feature < threshold"}, where
#'     thresholds are binned midpoints formatted in scientific notation.}
#' }
#' @keywords internal
.build_rule_cache <- function(tdt,
                              max_depth,
                              tighten_monotone = TRUE) {
  # Helper: numeric formatting for rule strings
  .fmt <- function(x) formatC(as.numeric(x), format = "e", digits = 2)

  # Tighten a rule string with monotonic split conditions.
  # If we encounter multiple constraints on the same feature in the same
  # direction, we keep only the most restrictive condition because weaker
  # conditions in the conjunction are redundant. Ordering is stable by first
  # appearance.
  tighten_prefix <- function(feat, dir, val) {
    key <-
      paste(feat, dir, sep = "\r") # stable key for (feature, direction)

    # Map key -> position in output vectors
    pos_map <- new.env(parent = emptyenv())

    first_pos <- integer(0)
    best_val <- numeric(0)
    best_feat <- character(0)
    best_dir <- character(0)

    for (i in seq_along(key)) {
      k <- key[[i]]
      j <- pos_map[[k]]

      if (is.null(j)) {
        j <- length(best_val) + 1L
        pos_map[[k]] <- j
        first_pos[j] <- i
        best_val[j] <- val[[i]]
        best_feat[j] <- feat[[i]]
        best_dir[j] <- dir[[i]]
      } else {
        # Update only if more restrictive in the same direction
        if (best_dir[[j]] == "<") {
          if (val[[i]] < best_val[[j]])
            best_val[[j]] <- val[[i]]
        } else if (best_dir[[j]] == ">=") {
          if (val[[i]] > best_val[[j]])
            best_val[[j]] <- val[[i]]
        } else {
          # Non-monotone direction (shouldn't happen in xgb-derived rules); keep first.
        }
      }
    }
    ord <- order(first_pos)
    paste0(best_feat[ord], " ", best_dir[ord], " ", .fmt(best_val[ord]))
  }

  # Build cache
  leaf_rule_list <- vector("list", 4096L) # pre-allocate per-leaf/per-prefix rows
  n_out <- 0L

  trees <- sort(unique(tdt$Tree))

  for (tt in trees) {
    tdt_t <- tdt[Tree == tt]
    if (!nrow(tdt_t))
      next

    # Map node ID -> row index
    max_id <- max(tdt_t$ID, na.rm = TRUE)
    row_of <- rep.int(NA_integer_, max_id + 1L)
    row_of[tdt_t$ID + 1L] <- seq_len(nrow(tdt_t))

    root_id <- 0L
    if (!(root_id %in% tdt_t$ID)) {
      root_id <- tdt_t$ID[which.min(tdt_t$ID)]
    }

    # Path buffers for rule representation: we store only the first max_depth
    # split-steps along a path; deeper splits are ignored for rule strings, but
    # still traversed to compute leaf_depth.
    path_feat <- character(max_depth)
    path_dir <- character(max_depth)
    path_val <- numeric(max_depth)

    walk <- function(node_id, depth_steps) {
      ridx <- row_of[node_id + 1L]
      if (is.na(ridx))
        return(invisible(NULL))

      if (isTRUE(tdt_t$Leaf[[ridx]])) {
        leaf_depth <- as.integer(depth_steps)

        max_keep <- min(max_depth, leaf_depth)

        if (max_keep >= 1L) {
          for (k in seq_len(max_keep)) {
            f <- path_feat[seq_len(k)]
            d <- path_dir[seq_len(k)]
            v <- path_val[seq_len(k)]

            if (isTRUE(tighten_monotone)) {
              conds <- tighten_prefix(f, d, v)
            } else {
              conds <- paste0(f, " ", d, " ", .fmt(v))
            }

            if (!length(conds))
              next
            rs <- paste(conds, collapse = " | ")

            n_out <<- n_out + 1L
            if (n_out > length(leaf_rule_list)) {
              leaf_rule_list <<-
                c(leaf_rule_list, vector("list", length(leaf_rule_list)))
            }

            leaf_rule_list[[n_out]] <<- data.table::data.table(
              Tree = as.integer(tt),
              leaf_id = as.integer(node_id),
              leaf_depth = as.integer(leaf_depth),
              rule_len = as.integer(k),
              n_clauses = as.integer(length(conds)),
              rule_str = as.character(rs)
            )
          }
        }
        return(invisible(NULL))
      }
      # Internal node: binned split
      feat <- tdt_t$Feature[[ridx]]
      split_val <- tdt_t$Split_bin[[ridx]]

      yes_id <- tdt_t$Yes[[ridx]]
      no_id <- tdt_t$No[[ridx]]

      # Next split-step depth (counting split conditions)
      next_depth <- depth_steps + 1L
      pos <- depth_steps + 1L # 1-indexed position in path buffers for this split-step

      # Store this step only if it's within the requested max_depth (prefix cap)
      if (pos <= max_depth) {
        path_feat[[pos]] <<- feat
        path_val[[pos]] <<- split_val
      }

      # Yes branch: "<"
      if (pos <= max_depth)
        path_dir[[pos]] <<- "<"
      if (is.finite(yes_id)) {
        walk(as.integer(yes_id), next_depth)
      }

      # No branch: ">="
      if (pos <= max_depth)
        path_dir[[pos]] <<- ">="
      if (is.finite(no_id)) {
        walk(as.integer(no_id), next_depth)
      }
      invisible(NULL)
    }
    walk(root_id, 0L)
  }

  leaf_rule_cache <-
    if (n_out) {
      data.table::rbindlist(leaf_rule_list[seq_len(n_out)], use.names = TRUE)
    } else {
      warning(sprintf(
        "[%s] failed to build rule cache.", "prepare_harvest"))
      return(invisible(NULL))
    }

  # Indexes (not keys) for fast joins without reordering
  data.table::setindex(leaf_rule_cache, rule_str)
  data.table::setindex(leaf_rule_cache, Tree, leaf_id)
  data.table::setindex(leaf_rule_cache, rule_len, rule_str)

  leaf_rule_cache[]
}


.build_rule_cache <- function(tdt,
                               max_depth,
                               tighten_monotone = TRUE) {
  # Helper: numeric formatting for rule strings
  .fmt <- function(x) formatC(as.numeric(x), format = "e", digits = 2)

  # Tighten a rule string with monotonic split conditions.
  # If we encounter multiple constraints on the same feature in the same
  # direction, we keep only the most restrictive condition because weaker
  # conditions in the conjunction are redundant. Ordering is stable by first
  # appearance.
  tighten_prefix <- function(feat, dir, val) {
    key <-
      paste(feat, dir, sep = "\r") # stable key for (feature, direction)

    # Map key -> position in output vectors
    pos_map <- new.env(parent = emptyenv())

    first_pos <- integer(0)
    best_val <- numeric(0)
    best_feat <- character(0)
    best_dir <- character(0)

    for (i in seq_along(key)) {
      k <- key[[i]]
      j <- pos_map[[k]]

      if (is.null(j)) {
        j <- length(best_val) + 1L
        pos_map[[k]] <- j
        first_pos[j] <- i
        best_val[j] <- val[[i]]
        best_feat[j] <- feat[[i]]
        best_dir[j] <- dir[[i]]
      } else {
        # Update only if more restrictive in the same direction
        if (best_dir[[j]] == "<") {
          if (val[[i]] < best_val[[j]])
            best_val[[j]] <- val[[i]]
        } else if (best_dir[[j]] == ">=") {
          if (val[[i]] > best_val[[j]])
            best_val[[j]] <- val[[i]]
        } else {
          # Non-monotone direction (shouldn't happen in xgb-derived rules); keep first.
        }
      }
    }
    ord <- order(first_pos)
    paste0(best_feat[ord], " ", best_dir[ord], " ", .fmt(best_val[ord]))
  }

  # Build cache
  leaf_rule_list <- vector("list", 4096L) # pre-allocate per-leaf/per-prefix rows
  n_out <- 0L

  trees <- sort(unique(tdt$Tree))

  for (tt in trees) {
    tdt_t <- tdt[Tree == tt]
    if (!nrow(tdt_t))
      next

    # Map node ID -> row index
    max_id <- max(tdt_t$ID, na.rm = TRUE)
    row_of <- rep.int(NA_integer_, max_id + 1L)
    row_of[tdt_t$ID + 1L] <- seq_len(nrow(tdt_t))

    root_id <- 0L
    if (!(root_id %in% tdt_t$ID)) {
      root_id <- tdt_t$ID[which.min(tdt_t$ID)]
    }

    # Path buffers for rule representation: we store only the first max_depth
    # split-steps along a path; deeper splits are ignored for rule strings, but
    # still traversed to compute leaf_depth.
    path_feat <- character(max_depth)
    path_dir <- character(max_depth)
    path_val <- numeric(max_depth)

    walk <- function(node_id, depth_steps) {
      ridx <- row_of[node_id + 1L]
      if (is.na(ridx))
        return(invisible(NULL))

      if (isTRUE(tdt_t$Leaf[[ridx]])) {
        leaf_depth <- as.integer(depth_steps)

        max_keep <- min(max_depth, leaf_depth)

        if (max_keep >= 1L) {
          for (k in seq_len(max_keep)) {
            f <- path_feat[seq_len(k)]
            d <- path_dir[seq_len(k)]
            v <- path_val[seq_len(k)]

            if (isTRUE(tighten_monotone)) {
              conds <- tighten_prefix(f, d, v)
            } else {
              conds <- paste0(f, " ", d, " ", .fmt(v))
            }

            if (!length(conds))
              next
            rs <- paste(conds, collapse = " | ")

            n_out <<- n_out + 1L
            if (n_out > length(leaf_rule_list)) {
              leaf_rule_list <<-
                c(leaf_rule_list, vector("list", length(leaf_rule_list)))
            }

            leaf_rule_list[[n_out]] <<- data.table::data.table(
              Tree = as.integer(tt),
              leaf_id = as.integer(node_id),
              leaf_depth = as.integer(leaf_depth),
              rule_len = as.integer(k),
              n_clauses = as.integer(length(conds)),
              rule_str = as.character(rs)
            )
          }
        }
        return(invisible(NULL))
      }
      # Internal node: binned split
      feat <- tdt_t$Feature[[ridx]]
      split_val <- tdt_t$Split_bin[[ridx]]

      yes_id <- tdt_t$Yes[[ridx]]
      no_id <- tdt_t$No[[ridx]]

      # Next split-step depth (counting split conditions)
      next_depth <- depth_steps + 1L
      pos <- depth_steps + 1L # 1-indexed position in path buffers for this split-step

      # Store this step only if it's within the requested max_depth (prefix cap)
      if (pos <= max_depth) {
        path_feat[[pos]] <<- feat
        path_val[[pos]] <<- split_val
      }

      # Yes branch: "<"
      if (pos <= max_depth)
        path_dir[[pos]] <<- "<"
      if (is.finite(yes_id)) {
        walk(as.integer(yes_id), next_depth)
      }

      # No branch: ">="
      if (pos <= max_depth)
        path_dir[[pos]] <<- ">="
      if (is.finite(no_id)) {
        walk(as.integer(no_id), next_depth)
      }
      invisible(NULL)
    }
    walk(root_id, 0L)
  }

  leaf_rule_cache <-
    if (n_out) {
      data.table::rbindlist(leaf_rule_list[seq_len(n_out)], use.names = TRUE)
    } else {
        warning(sprintf(
          "[%s] failed to build rule cache.", "prepare_harvest"))
        return(invisible(NULL))
    }

  # Indexes (not keys) for fast joins without reordering
  data.table::setindex(leaf_rule_cache, rule_str)
  data.table::setindex(leaf_rule_cache, Tree, leaf_id)
  data.table::setindex(leaf_rule_cache, rule_len, rule_str)

  leaf_rule_cache[]
}

