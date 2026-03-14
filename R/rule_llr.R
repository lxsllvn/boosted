#' Compute a rule-level log-likelihood ratio with Jeffreys smoothing
#'
#' Computes \eqn{\log(p_E / p_B)}{} where \eqn{p_E}{} and \eqn{p_B}{} are the
#' smoothed probabilities that an extreme or background SNP (respectively) falls
#' within a rule's SNP bucket. The same Jeffreys-prior scheme used in
#' \code{\link{.leaf_llrs}} is applied: counts are incremented by \code{alpha}
#' and denominators by \code{2 * alpha}, so that a rule capturing no labelled
#' SNPs of either class returns a LLR near zero rather than \code{-Inf} or
#' \code{Inf}.
#'
#' @param n_extreme Integer scalar: number of extreme SNPs in the rule's bucket.
#' @param n_bg Integer scalar: number of background SNPs in the rule's bucket.
#' @param N_extr_total Integer scalar: total number of extreme SNPs in the
#'   relevant partition (the denominator for \eqn{p_E}{}).
#' @param N_bg_total Integer scalar: total number of background SNPs in the
#'   relevant partition (the denominator for \eqn{p_B}{}).
#' @param alpha Numeric scalar \eqn{\ge 0}: Jeffreys-prior concentration.
#'   When \code{alpha = 0} a small epsilon floor is applied; see
#'   \code{\link{.boosted_params}}.
#'
#' @return Numeric scalar: the smoothed log-likelihood ratio for the rule.
#'   Positive values indicate enrichment for extreme SNPs; negative values
#'   indicate depletion.
#'
#' @keywords internal

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
