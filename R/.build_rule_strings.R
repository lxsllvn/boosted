#' Title
#'
#' @param leaf_paths
#' @param bin_spec
#' @param max_depth
#' @param tighten_monotone
#'
#' @return
#' @keywords internal
#'
#' @examples

.build_rule_strings <- function(leaf_paths,
                                bin_spec,
                                max_depth = NULL,
                                tighten_monotone = TRUE) {
  LS <- data.table::as.data.table(leaf_paths)

  # Map split_val -> thresh_bin using bin_spec$mids / bin_spec$breaks
  LS[, thresh_bin := {
    f  <- feature[1L]
    br <- bin_spec$breaks[[f]]
    md <- bin_spec$mids[[f]]

    if (is.null(br) || is.null(md) || !is.finite(split_val[1L])) {
      rep(NA_real_, .N)
    } else {
      idx <- findInterval(split_val, br, all.inside = TRUE)
      as.numeric(md[idx])
    }
  }, by = feature]

  LS[, `:=`(
    Tree       = as.integer(Tree),
    leaf_id    = as.integer(leaf_id),
    depth      = as.integer(depth),
    feature    = as.character(feature),
    direction  = as.character(direction),
    thresh_bin = as.numeric(thresh_bin)
  )]

  LS <- LS[depth < as.integer(max_depth)]
  data.table::setorder(LS, Tree, leaf_id, depth)

  if (isTRUE(tighten_monotone)) {
    eps <- 1e-12

    PATHS <- LS[, {
      # We will keep at most one "<" and one ">=" per feature.
      out_cond    <- character(0)
      out_depth   <- integer(0)
      # strongest ">=" bound per feature
      best_ge_val <- new.env(parent = emptyenv())
      # position in out_cond/out_depth
      best_ge_pos <- new.env(parent = emptyenv())
      # strongest "<" bound per feature
      best_lt_val <- new.env(parent = emptyenv())
      # position in out_cond/out_depth
      best_lt_pos <- new.env(parent = emptyenv())

      for (i in seq_len(.N)) {
        f <- feature[i]
        d <- direction[i]
        b <- thresh_bin[i]
        if (!is.finite(b))
          next

        cond_str <- paste0(f, " ", d, " ", .format_num(b))

        if (d == ">=") {
          cur_val <- best_ge_val[[f]]
          cur_pos <- best_ge_pos[[f]]

          if (is.null(cur_val)) {
            # First ">=" constraint on this feature: append
            out_cond  <- c(out_cond, cond_str)
            out_depth <- c(out_depth, depth[i])
            pos <- length(out_cond)
            best_ge_val[[f]] <- b
            best_ge_pos[[f]] <- pos
          } else if (b > cur_val + eps) {
            # Stronger lower bound: overwrite previous constraint
            pos <- cur_pos
            out_cond [pos] <- cond_str
            out_depth[pos] <- depth[i]
            best_ge_val[[f]] <- b
          } else {
            # Weaker or equal ">=" bound; skip
            next
          }

        } else if (d == "<") {
          cur_val <- best_lt_val[[f]]
          cur_pos <- best_lt_pos[[f]]

          if (is.null(cur_val)) {
            # First "<" constraint on this feature: append
            out_cond  <- c(out_cond, cond_str)
            out_depth <- c(out_depth, depth[i])
            pos <- length(out_cond)
            best_lt_val[[f]] <- b
            best_lt_pos[[f]] <- pos
          } else if (b < cur_val - eps) {
            # Stronger upper bound: overwrite previous constraint
            pos <- cur_pos
            out_cond [pos] <- cond_str
            out_depth[pos] <- depth[i]
            best_lt_val[[f]] <- b
          } else {
            # Weaker or equal "<" bound; skip
            next
          }

        } else {
          # Non-monotone direction (e.g. "missing"): always keep as-is
          out_cond  <- c(out_cond, cond_str)
          out_depth <- c(out_depth, depth[i])
        }
      }

      if (!length(out_cond)) {
        NULL
      } else {
        .(cond = out_cond, depth = out_depth)
      }
    }, by = .(Tree, leaf_id)]

  } else {
    # No tightening: keep full sequence as-is
    PATHS <-
      LS[, .(cond = paste0(feature, " ", direction, " ", .format_num(thresh_bin)),
             depth = depth), by = .(Tree, leaf_id)]
  }

  PATHS_out <- PATHS[, .(
    rule_len  = .N,
    rule_str  = paste(cond, collapse = " | "),
    depth_max = max(depth)
  ), by = .(Tree, leaf_id)]

  data.table::setkey(PATHS_out, Tree, leaf_id)
  PATHS_out[]
}
