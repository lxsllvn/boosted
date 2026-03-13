#' Compute a gain curve for a scored set of SNPs
#'
#' Sorts SNPs by score descending (with \code{NA}s ranked last), then
#' evaluates the gain curve at each fraction in \code{grid}: for the top
#' \code{ceiling(f * n)} SNPs by score, how many extreme SNPs are recovered
#' (recall) relative to random screening at the same fraction? The ratio of
#' recall to screening fraction is the lift at that point.
#'
#' @param scores Numeric vector of length \code{n} containing per-SNP scores.
#'   \code{NA}s are treated as the lowest possible score.
#' @param n Integer scalar: total number of SNPs (length of \code{scores}).
#' @param extr_idx Integer vector of 1-based indices of the extreme (positive)
#'   SNPs within \code{scores}.
#' @param grid Numeric vector of screening fractions in \code{(0, 1]} at
#'   which to evaluate the gain curve; typically the validated
#'   \code{gain_grid} from the calling function.
#'
#' @return A \code{data.table} with one row per element of \code{grid} and
#'   the following columns:
#' \describe{
#'   \item{\code{frac_screened}}{Numeric. The screening fraction from
#'     \code{grid}.}
#'   \item{\code{n_screened}}{Integer. Number of SNPs in the top-scoring
#'     fraction (\code{ceiling(frac_screened * n)}, clipped to
#'     \code{[1, n]}).}
#'   \item{\code{score_threshold}}{Numeric. The score of the last SNP
#'     included at this cutoff.}
#'   \item{\code{recall}}{Numeric. Fraction of extreme SNPs recovered among
#'     the top \code{n_screened} SNPs.}
#'   \item{\code{lift_curve}}{Numeric. \code{recall / frac_screened}; equals
#'     1 for a random ranker.}
#' }
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
