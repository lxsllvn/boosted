#' Build a rule cache from a parsed `xgboost` tree table
#'
#' Traverses every tree in \code{tdt} depth-first, recording the sequence of
#' split conditions encountered on the path from root to each leaf. For every
#' leaf, all prefix rules of length \eqn{1, \ldots, \min(\text{max\_depth},
#' \text{leaf\_depth})}{} are materialized as character strings and stored
#' alongside the leaf's \code{(Tree, leaf_id)} identity. Split thresholds are
#' taken directly from the \code{Split} column,
#' formatted to 4 significant figures in scientific notation. When
#' \code{tighten_monotone = TRUE}, multiple constraints on the same feature
#' and direction within a single prefix are collapsed to the tightest bound,
#' removing logically redundant clauses.
#'
#' @param tdt A \code{data.table} produced by \code{\link{.parse_xgb_tree}}.
#'   Required columns: \code{Tree}, \code{ID}, \code{Feature}, \code{Split},
#'   \code{Yes}, \code{No}, and \code{Leaf}.
#' @param max_depth Positive integer. Maximum number of split-steps to include
#'   in any rule prefix. Rules for leaves deeper than \code{max_depth} are
#'   truncated at \code{max_depth} splits; \code{leaf_depth} is still recorded
#'   faithfully for all leaves regardless of the cap.
#' @param tighten_monotone Logical. If \code{TRUE} (default), redundant
#'   monotone constraints on the same feature and direction within a prefix are
#'   collapsed to the tightest bound. For example, \code{x < 5 & x < 3}
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
#'     ampersand-delimited sequence of split conditions of the form
#'     \code{"feature >= threshold"} or \code{"feature < threshold"}, where
#'     thresholds are raw `xgboost` split values formatted to 4 significant
#'     figures in scientific notation.}
#' }
#' @keywords internal
.build_rule_cache <- function(tdt,
                              max_depth,
                              tighten_monotone = TRUE) {
  # Helper: numeric formatting for rule strings.
  # 4 significant figures in scientific notation gives enough precision to
  # distinguish distinct raw xgboost splits without making strings unreadable.
  .fmt <- function(x) formatC(as.numeric(x), format = "e", digits = 4)

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
    tightened_feat <- best_feat[ord]
    tightened_dir  <- best_dir[ord]
    tightened_val  <- best_val[ord]

    # Feasibility check: for any feature appearing with both a >= and a <
    # bound, the rule is satisfiable only if the >= bound is strictly less
    # than the < bound. Binning can collapse distinct raw thresholds to the
    # same midpoint, producing contradictions like f >= -1.5 & f < -1.5.
    # Return character(0) to signal an infeasible rule; the caller skips it
    # via the existing if (!length(conds)) next guard.
    lt_feats <- tightened_feat[tightened_dir == "<"]
    ge_feats <- tightened_feat[tightened_dir == ">="]
    shared   <- intersect(lt_feats, ge_feats)
    if (length(shared)) {
      for (f in shared) {
        val_lt <- tightened_val[tightened_feat == f & tightened_dir == "<"]
        val_ge <- tightened_val[tightened_feat == f & tightened_dir == ">="]
        if (val_ge >= val_lt)
          return(character(0))
      }
    }

    paste0(tightened_feat, " ", tightened_dir, " ", .fmt(tightened_val))
  }

  # Hoist the tighten_monotone branch out of the inner loop. build_conds
  # is defined once here and closes over the constant tighten_monotone value,
  # returning the vector of individual clause strings so the caller can both
  # count them (for n_clauses) and collapse them into a rule string.
  build_conds <- if (isTRUE(tighten_monotone)) {
    function(f, d, v) tighten_prefix(f, d, v)
  } else {
    function(f, d, v) paste0(f, " ", d, " ", .fmt(v))
  }

  # Pre-allocate typed output vectors with an initial capacity of 4096 rows,
  # doubling when full. This replaces the per-row data.table() construction
  # and the rbindlist() at the end with a single data.table() call.
  cap             <- 4096L
  n_out           <- 0L
  out_Tree        <- integer(cap)
  out_leaf_id     <- integer(cap)
  out_leaf_depth  <- integer(cap)
  out_rule_len    <- integer(cap)
  out_n_clauses   <- integer(cap)
  out_rule_str    <- character(cap)

  .grow <- function() {
    cap        <<- cap * 2L
    out_Tree       <<- c(out_Tree,       integer(cap / 2L))
    out_leaf_id    <<- c(out_leaf_id,    integer(cap / 2L))
    out_leaf_depth <<- c(out_leaf_depth, integer(cap / 2L))
    out_rule_len   <<- c(out_rule_len,   integer(cap / 2L))
    out_n_clauses  <<- c(out_n_clauses,  integer(cap / 2L))
    out_rule_str   <<- c(out_rule_str,   character(cap / 2L))
  }

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
        max_keep   <- min(max_depth, leaf_depth)

        if (max_keep >= 1L) {
          for (k in seq_len(max_keep)) {
            idx   <- seq_len(k)
            conds <- build_conds(path_feat[idx], path_dir[idx], path_val[idx])

            if (!length(conds))
              next
            rs <- paste(conds, collapse = " & ")

            n_out <<- n_out + 1L
            if (n_out > cap)
              .grow()

            out_Tree[n_out]       <<- as.integer(tt)
            out_leaf_id[n_out]    <<- as.integer(node_id)
            out_leaf_depth[n_out] <<- leaf_depth
            out_rule_len[n_out]   <<- as.integer(k)
            out_n_clauses[n_out]  <<- as.integer(length(conds))
            out_rule_str[n_out]   <<- rs
          }
        }
        return(invisible(NULL))
      }
      # Internal node: raw xgboost split threshold
      feat      <- tdt_t$Feature[[ridx]]
      split_val <- tdt_t$Split[[ridx]]
      yes_id    <- tdt_t$Yes[[ridx]]
      no_id     <- tdt_t$No[[ridx]]

      # 1-indexed position in path buffers for this split-step
      pos        <- depth_steps + 1L
      within_cap <- pos <= max_depth

      # Store this step only if it's within the requested max_depth (prefix cap)
      if (within_cap) {
        path_feat[[pos]] <<- feat
        path_val[[pos]]  <<- split_val
      }

      # Yes branch: "<"
      if (within_cap)
        path_dir[[pos]] <<- "<"
      if (is.finite(yes_id))
        walk(as.integer(yes_id), pos)

      # No branch: ">="
      if (within_cap)
        path_dir[[pos]] <<- ">="
      if (is.finite(no_id))
        walk(as.integer(no_id), pos)
      invisible(NULL)
    }
    walk(root_id, 0L)
  }

  if (!n_out) {
    warning(sprintf("[%s] failed to build rule cache.", "prepare_harvest"))
    return(invisible(NULL))
  }

  # Build output table in one shot from the pre-allocated vectors
  leaf_rule_cache <- data.table::data.table(
    Tree       = out_Tree[seq_len(n_out)],
    leaf_id    = out_leaf_id[seq_len(n_out)],
    leaf_depth = out_leaf_depth[seq_len(n_out)],
    rule_len   = out_rule_len[seq_len(n_out)],
    n_clauses  = out_n_clauses[seq_len(n_out)],
    rule_str   = out_rule_str[seq_len(n_out)]
  )

  # Indexes (not keys) for fast joins without reordering
  data.table::setindex(leaf_rule_cache, rule_str)
  data.table::setindex(leaf_rule_cache, Tree, leaf_id)
  data.table::setindex(leaf_rule_cache, rule_len, rule_str)

  leaf_rule_cache[]
}
