#' Title
#'
#' @param scores
#' @param n
#' @param extr_idx
#' @param grid
#'
#' @return
#' @keywords internal
#'
#' @examples

.gain_curve <- function(scores,
                        n,
                        extr_idx,
                        grid) {
  # Precompute targets
  n_screened <- pmax(1L, pmin(n, ceiling(grid * n)))

  # Fast full sort (radix); NAs sink to end
  si <- sort.int(
    scores,
    decreasing   = TRUE,
    na.last      = TRUE,
    method       = "radix",
    index.return = TRUE
  )
  sc_sorted <- si$x
  ord       <- si$ix

  # Class mask in sorted order
  is_pos           <- logical(n)
  is_pos[extr_idx] <- TRUE
  is_pos_sorted    <- is_pos[ord]

  # Cumulative true positives and derived metrics
  cum_tp     <- cumsum(is_pos_sorted)
  tp         <- cum_tp[n_screened]
  recall     <- tp / length(extr_idx)
  lift_curve <- recall / grid

  # Score threshold at each cutoff
  score_threshold <- sc_sorted[n_screened]

  data.table::data.table(
    frac_screened   = grid,
    n_screened      = n_screened,
    score_threshold = score_threshold,
    recall          = recall,
    lift_curve      = lift_curve
  )
}
