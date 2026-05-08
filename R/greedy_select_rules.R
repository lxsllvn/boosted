#' Greedy selection of a small non-redundant rule set
#'
#' Iteratively selects rules from the antichain produced by
#' \code{\link{build_rule_dag}} to maximize extreme SNP coverage while
#' penalizing background SNP inclusion. At each step the rule with the
#' highest marginal score is selected, where the score is
#' \eqn{|\Delta E| - \lambda \cdot |\Delta B|}{} with \eqn{|\Delta E|}{} and
#' \eqn{|\Delta B|}{} being the number of extreme and background SNPs newly
#' covered by that rule (i.e. not already covered by previously selected
#' rules). Selection stops when any of the supplied stopping criteria is met.
#'
#' Because coverage is tracked for both \eqn{E}{} and \eqn{B}{} incrementally,
#' the background penalty correctly counts only the unique marginal background
#' added at each step — double-counting of background SNPs shared between
#' selected rules does not occur.
#'
#' @param dag_result Named list returned by \code{\link{build_rule_dag}}.
#' @param lambda Non-negative numeric scalar. Penalty weight on marginal
#'   background SNP coverage. A value of \code{1} treats each background SNP
#'   equally with each extreme SNP. Values below \code{1} prioritize recall;
#'   values above \code{1} prioritize precision. Default \code{1.0}.
#' @param max_rules Positive integer or \code{NULL}. Maximum number of rules
#'   to select. \code{NULL} (default) places no explicit limit.
#' @param min_marginal_ext Non-negative integer. Selection stops when the best
#'   candidate would add no more than this many new extreme SNPs. Default
#'   \code{0L}: stop only when no new extreme SNPs can be added.
#' @param target_coverage_frac Numeric scalar in \code{(0, 1]}. Selection
#'   stops once this fraction of the total extreme SNP set (\code{A_n}) has
#'   been covered. Default \code{1.0}.
#'
#' @return \code{NULL} (invisibly) if no rules were selected, otherwise a
#'   named list with:
#' \describe{
#'   \item{\code{selected_rules}}{Character vector of selected rule strings,
#'     in selection order.}
#'   \item{\code{ledger}}{A \code{data.table} with one row per selection step
#'     and columns: \code{step}, \code{rule_str}, \code{marginal_ext}
#'     (new extreme SNPs covered), \code{marginal_bg} (new background SNPs
#'     added), \code{cumulative_ext} (total extreme SNPs covered so far),
#'     \code{coverage_frac} (\code{cumulative_ext / A_n}),
#'     \code{cumulative_bg} (total background SNPs covered so far), and
#'     \code{score} (the greedy objective value for this step).}
#'   \item{\code{n_selected}}{Integer. Number of rules selected.}
#'   \item{\code{coverage_frac}}{Numeric. Final fraction of extreme SNPs
#'     covered.}
#'   \item{\code{A_n}, \code{B_n}}{Integers. Total extreme and background SNP
#'     set sizes.}
#'   \item{\code{lambda}}{The penalty weight used.}
#' }
#' @export
greedy_select_rules <- function(dag_result,
                                lambda               = 1.0,
                                max_rules            = NULL,
                                min_marginal_ext     = 0L,
                                target_coverage_frac = 1.0) {
  FUN <- "greedy_select_rules"

  # Validate dag_result
  req <- c("antichain_rule_ids", "antichain_sets", "n_antichain")
  missing <- setdiff(req, names(dag_result))
  if (length(missing)) {
    stop(sprintf(
      "[%s] dag_result missing elements: %s. Re-run build_rule_dag().",
      FUN, paste(missing, collapse = ", ")
    ))
  }

  cand_A <- dag_result$antichain_sets$A_sets   # list; cand_A[[j]] = compact extreme IDs
  cand_B <- dag_result$antichain_sets$B_sets   # list; cand_B[[j]] = compact bg IDs
  A_n    <- dag_result$antichain_sets$A_n
  B_n    <- dag_result$antichain_sets$B_n
  rule_ids_antichain <- dag_result$antichain_rule_ids

  K <- dag_result$n_antichain
  if (!K) stop(sprintf("[%s] antichain is empty; nothing to select.", FUN))
  if (A_n == 0L) stop(sprintf("[%s] A_n = 0; no extreme SNPs to cover.", FUN))

  # Resolve stopping criteria
  effective_max   <- if (is.null(max_rules)) K else as.integer(max_rules)
  target_count    <- ceiling(target_coverage_frac * A_n)
  min_marg        <- as.integer(min_marginal_ext)

  # State: which extreme/background SNPs are already covered
  covered_ext <- logical(A_n)
  covered_bg  <- logical(B_n)

  # Candidates not yet selected (1..K positions into antichain)
  remaining <- seq_len(K)

  # Pre-allocate ledger rows
  ledger_rows <- vector("list", effective_max)
  n_selected  <- 0L
  selected    <- integer(effective_max)

  message(sprintf(
    "[%s] start: %d antichain candidates, A_n=%d, lambda=%.3g",
    FUN, K, A_n, lambda
  ))

  for (iter in seq_len(effective_max)) {
    if (!length(remaining)) break

    # Score each remaining candidate by marginal extreme coverage minus
    # lambda * marginal background coverage.
    scores <- vapply(remaining, function(j) {
      mext <- if (length(cand_A[[j]])) sum(!covered_ext[cand_A[[j]]]) else 0L
      mbg  <- if (length(cand_B[[j]])) sum(!covered_bg[cand_B[[j]]])  else 0L
      mext - lambda * mbg
    }, numeric(1L))

    best_local <- which.max(scores)
    best_j     <- remaining[best_local]

    marg_ext <- if (length(cand_A[[best_j]])) sum(!covered_ext[cand_A[[best_j]]]) else 0L
    marg_bg  <- if (length(cand_B[[best_j]])) sum(!covered_bg[cand_B[[best_j]]])  else 0L

    # Stopping: best candidate adds no new extreme SNPs (at or below threshold)
    if (marg_ext <= min_marg) {
      message(sprintf(
        "[%s] stopped at step %d: marginal extreme coverage %d <= min_marginal_ext %d.",
        FUN, iter, marg_ext, min_marg
      ))
      break
    }

    # Accept selection; update covered sets
    n_selected           <- n_selected + 1L
    selected[n_selected] <- best_j
    remaining            <- remaining[-best_local]

    if (length(cand_A[[best_j]])) covered_ext[cand_A[[best_j]]] <- TRUE
    if (length(cand_B[[best_j]])) covered_bg[cand_B[[best_j]]]  <- TRUE

    cum_ext <- sum(covered_ext)
    cum_bg  <- sum(covered_bg)

    ledger_rows[[n_selected]] <- data.table::data.table(
      step           = as.integer(iter),
      rule_str       = rule_ids_antichain[best_j],
      marginal_ext   = as.integer(marg_ext),
      marginal_bg    = as.integer(marg_bg),
      cumulative_ext = as.integer(cum_ext),
      coverage_frac  = cum_ext / A_n,
      cumulative_bg  = as.integer(cum_bg),
      score          = as.numeric(scores[best_local])
    )

    if (cum_ext >= target_count) {
      message(sprintf(
        "[%s] target coverage %.1f%% reached at step %d.",
        FUN, 100 * target_coverage_frac, iter
      ))
      break
    }
  }

  if (!n_selected) {
    warning(sprintf("[%s] no rules selected.", FUN))
    return(invisible(NULL))
  }

  ledger <- data.table::rbindlist(ledger_rows[seq_len(n_selected)])

  message(sprintf(
    "[%s] selected %d rules; extreme coverage %.1f%% (%d / %d).",
    FUN, n_selected, 100 * sum(covered_ext) / A_n, sum(covered_ext), A_n
  ))

  list(
    selected_rules = rule_ids_antichain[selected[seq_len(n_selected)]],
    ledger         = ledger[],
    n_selected     = n_selected,
    coverage_frac  = sum(covered_ext) / A_n,
    A_n            = A_n,
    B_n            = B_n,
    lambda         = lambda
  )
}
