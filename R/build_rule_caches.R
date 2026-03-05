#' Build a rule cache from a parsed xgboost tree table
#'
#' @description
#'
#' @param tdt the parsed `xgboost` tree table produced by `boosted::prepare_harvest`, with columns
#'   `Tree`, `ID`, `Feature`, `Split`, `Yes`, `No`, `Missing`, `Leaf`, `Split_bin`.
#' @param max_depth maximum number of decision nodes to include in rule strings. Must be <=
#' @param tighten_monotone If `TRUE`, collapse repeated monotone bounds per
#'   feature and direction to the tightest bound.
#'
#' @return A `data.table` with one row per `(Tree, leaf_id, rule_len)`, where `rule_len` is an integer ranging from 1, ..., min(max_depth, leaf_depth).
#' * `Tree`: unique integer ID denoting each tree in the model (zero-based
#'     index)
#'
#' * `leaf_id`: unique integer ID of each leaf in a given `Tree`.
#'
#' * `leaf_depth`: the actual number of decision nodes traversed from the root to this leaf.
#'
#' * `rule_len`: the number of decision nodes represented by this rule.
#'
#' * `n_clauses`: the realized length of the rule after tightening monotonic constraints.
#'
#' * `rule_str`: a character string of the ordered split conditions comprising this rule. Rule strings follow the format of "feature {< | => } threshold" and multiple conditions are joined by "|".
#'
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

