#' Title
#'
#' @param n_extreme
#' @param n_bg
#' @param N_extr_total
#' @param N_bg_total
#' @param alpha
#'
#' @return
#' @keywords internal
#'
#' @examples

.rule_llr <- function(n_extreme,
                      n_bg,
                      N_extr_total,
                      N_bg_total,
                      alpha = 0.5) {
  # Jeffreys prior smoothing on rule-level contingency counts
  pE <- (n_extreme + alpha) / (N_extr_total + 2 * alpha)
  pB <- (n_bg      + alpha) / (N_bg_total   + 2 * alpha)
  log(pE / pB)
}
