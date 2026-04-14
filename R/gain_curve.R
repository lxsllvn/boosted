#' Compute recall and labeled precision along a gain curve for scored SNPs
#'
#' Sorts SNPs by score in descending order (with \code{NA}s ranked last), then
#' evaluates performance at each screening fraction in \code{grid}. For the top
#' \code{ceiling(f * n)} SNPs by score, the function computes: \itemize{ \item
#' \code{recall}: the fraction of labeled extreme SNPs recovered; and \item
#' \code{precision_vs_bg}: the fraction of labeled SNPs recovered that belong to
#' the extreme set rather than the labeled background set. } Unlabeled SNPs are
#' included in the ranking and may occupy screened positions, but they are
#' excluded from labeled-set counts. Thus, recall is computed only with respect
#' to the labeled extreme set, and precision_vs_bg is computed only among
#' labeled extreme and labeled background SNPs recovered at each cutoff.
#'
#' @param scores Numeric vector of length \code{n} containing per-SNP scores.
#'   \code{NA}s are treated as the lowest possible score.
#' @param n Integer scalar: total number of SNPs (length of \code{scores}).
#' @param extr_idx Integer vector of 1-based indices of the labeled extreme SNPs
#'   within \code{scores}.
#' @param bg_idx Integer vector of 1-based indices of the labeled background
#'   SNPs within \code{scores}.
#' @param grid Numeric vector of screening fractions in \code{(0, 1]} at which
#'   to evaluate the gain curve; typically the validated \code{gain_grid} from
#'   the calling function.
#'
#' @return A \code{data.table} with one row per element of \code{grid} and the
#'   following columns: \describe{ \item{\code{frac_screened}}{Numeric. The
#'   screening fraction from \code{grid}.} \item{\code{n_screened}}{Integer.
#'   Number of SNPs in the top-scoring fraction (\code{ceiling(frac_screened *
#'   n)}, clipped to \code{[1, n]}).} \item{\code{score_threshold}}{Numeric. The
#'   score of the last SNP included at this cutoff.} \item{\code{tp}}{Integer.
#'   Number of labeled extreme SNPs recovered among the top \code{n_screened}
#'   SNPs.} \item{\code{fp_bg}}{Integer. Number of labeled background SNPs
#'   recovered among the top \code{n_screened} SNPs.}
#'   \item{\code{recall}}{Numeric. Fraction of labeled extreme SNPs recovered
#'   among the top \code{n_screened} SNPs.}
#'   \item{\code{precision_vs_bg}}{Numeric. Among recovered SNPs that belong to
#'   either the labeled extreme or labeled background set, the fraction
#'   belonging to the extreme set: \code{tp / (tp + fp_bg)}. Returns \code{NA}
#'   when no labeled SNPs have been recovered at that cutoff.}
#'   \item{\code{lift_curve}}{Numeric. \code{recall / frac_screened}; equals 1
#'   for a random ranker with respect to recovery of the extreme set.} }
#' @keywords internal

.gain_curve <- function(scores,
                        n,
                        extr_idx,
                        bg_idx,
                        grid) {
  # Precompute targets
  n_screened <- pmax(1L, pmin(n, ceiling(grid * n)))

  # Sort scores descending; NAs sink to end
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
  is_extr <- logical(n)
  is_extr[extr_idx] <- TRUE

  is_bg <- logical(n)
  if (length(bg_idx)) {
    is_bg[bg_idx] <- TRUE
  }

  # Sorted masks
  is_extr_sorted <- is_extr[ord]
  is_bg_sorted   <- is_bg[ord]

  # Cumulative counts
  cum_tp <- cumsum(is_extr_sorted)
  cum_fp <- cumsum(is_bg_sorted)

  tp <- cum_tp[n_screened]
  fp <- cum_fp[n_screened]

  recall <- tp / length(extr_idx)

  denom_labeled <- tp + fp
  precision_vs_bg <- ifelse(denom_labeled > 0, tp / denom_labeled, NA_real_)

  lift_curve <- recall / grid
  score_threshold <- sc_sorted[n_screened]

  data.table::data.table(
    frac_screened   = grid,
    n_screened      = n_screened,
    score_threshold = score_threshold,
    tp              = tp,
    fp_bg           = fp,
    recall          = recall,
    precision_vs_bg = precision_vs_bg,
    lift_curve      = lift_curve
  )
}
