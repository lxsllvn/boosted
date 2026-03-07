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
  # Near-empirical LLRs (alpha == 0) with tiny epsilon to avoid log(0)
  if (alpha == 0) {
    eps <- 1e-12
    pE_raw <- (n_extreme + alpha) / (N_extr_total + 2 * alpha)
    pB_raw <- (n_bg      + alpha) / (N_bg_total   + 2 * alpha)
    # Epsilon floor to avoid zeros
    pE <- pmax(pE_raw, eps)
    pB <- pmax(pB_raw, eps)
    log(pE / pB)
  } else {
    # Jeffreys prior smoothing on rule-level contingency counts
    pE <- (n_extreme + alpha) / (N_extr_total + 2 * alpha)
    pB <- (n_bg      + alpha) / (N_bg_total   + 2 * alpha)
    log(pE / pB)}
}
