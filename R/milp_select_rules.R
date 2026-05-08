#' MILP-based selection of an optimal rule set
#'
#' Finds the globally optimal subset of antichain rules (produced by
#' \code{\link{build_rule_dag}}) that maximises extreme SNP coverage while
#' penalising background SNP inclusion, using mixed-integer linear programming
#' via \pkg{Rglpk}. The formulation introduces one binary variable \eqn{z_i}{}
#' per extreme SNP (is it covered?) and one binary variable \eqn{x_k}{} per
#' candidate rule (is it selected?), with coverage constraints enforcing
#' \eqn{z_i \le \sum_{k : i \in A_k} x_k}{}.
#'
#' **Background penalty note.** The objective penalises background using total
#' incidences: \eqn{\lambda \sum_k n\_bg_k \cdot x_k}{}. This counts
#' background SNPs covered by multiple selected rules multiple times, which
#' overpenalises when selected rules share background coverage. The
#' overpenalisation is conservative and becomes small when the antichained rule
#' set has low pairwise background overlap (which well-applied dominance
#' pruning tends to produce). If precise unique background coverage is
#' required, use \code{\link{greedy_select_rules}}, which tracks marginal
#' background exactly.
#'
#' @param dag_result Named list returned by \code{\link{build_rule_dag}}.
#' @param lambda Non-negative numeric scalar. Penalty weight on background
#'   incidences per selected rule. The objective is
#'   \eqn{\sum_i z_i - \lambda \sum_k n\_bg_k \cdot x_k}{}. Default
#'   \code{1.0}.
#' @param max_rules Positive integer or \code{NULL}. Optional cardinality
#'   constraint: at most \code{max_rules} rules may be selected. \code{NULL}
#'   (default) imposes no limit.
#'
#' @return \code{NULL} (invisibly) if the solver selects no rules, otherwise
#'   a named list with:
#' \describe{
#'   \item{\code{selected_rules}}{Character vector of selected rule strings.
#'     Order is determined by a post-hoc greedy simulation of the selected
#'     rules (decreasing marginal extreme coverage) so that the ledger is
#'     directly comparable with \code{\link{greedy_select_rules}} output.}
#'   \item{\code{ledger}}{A \code{data.table} with one row per selected rule
#'     in the simulated greedy order, with the same columns as the ledger
#'     returned by \code{\link{greedy_select_rules}} (\code{step},
#'     \code{rule_str}, \code{marginal_ext}, \code{marginal_bg},
#'     \code{cumulative_ext}, \code{coverage_frac}, \code{cumulative_bg}).
#'     The \code{score} column is \code{NA} (the MILP selects jointly, not
#'     greedily).}
#'   \item{\code{n_selected}}{Integer. Number of rules selected.}
#'   \item{\code{coverage_frac}}{Numeric. Fraction of extreme SNPs covered
#'     by the selected rule set.}
#'   \item{\code{objective}}{Numeric. The optimal objective value (on the
#'     maximisation scale: \eqn{\sum z_i - \lambda \sum n\_bg_k x_k}{}).}
#'   \item{\code{solver_status}}{Integer. Status code from
#'     \code{Rglpk_solve_LP}: 0 indicates optimal.}
#'   \item{\code{A_n}, \code{B_n}}{Integers. Total extreme and background
#'     SNP set sizes.}
#'   \item{\code{lambda}}{The penalty weight used.}
#' }
#' @export
milp_select_rules <- function(dag_result,
                               lambda    = 1.0,
                               max_rules = NULL) {
  FUN <- "milp_select_rules"

  if (!requireNamespace("Rglpk", quietly = TRUE)) {
    stop(sprintf(
      "[%s] Rglpk is required. Install with install.packages('Rglpk').", FUN
    ))
  }

  # Validate dag_result
  req <- c("antichain_rule_ids", "antichain_sets", "antichain_summary",
           "n_antichain")
  missing_el <- setdiff(req, names(dag_result))
  if (length(missing_el)) {
    stop(sprintf(
      "[%s] dag_result missing elements: %s. Re-run build_rule_dag().",
      FUN, paste(missing_el, collapse = ", ")
    ))
  }

  cand_A   <- dag_result$antichain_sets$A_sets
  cand_B   <- dag_result$antichain_sets$B_sets
  A_n      <- dag_result$antichain_sets$A_n
  B_n      <- dag_result$antichain_sets$B_n
  rule_ids <- dag_result$antichain_rule_ids
  cand_bg  <- dag_result$antichain_summary$n_bg   # length K

  K <- dag_result$n_antichain
  if (!K) stop(sprintf("[%s] antichain is empty; nothing to select.", FUN))
  if (A_n == 0L) stop(sprintf("[%s] A_n = 0; no extreme SNPs to cover.", FUN))

  # -------------------------------------------------------------------
  # Variable layout: [z_1, ..., z_{A_n}, x_1, ..., x_K]
  #
  # z_i in {0,1}: is extreme SNP i covered by the selected rule set?
  # x_k in {0,1}: is antichain rule k selected?
  #
  # Objective (minimise): -sum(z_i) + lambda * sum(n_bg_k * x_k)
  #   => maximise  sum(z_i) - lambda * sum(n_bg_k * x_k)
  #
  # Coverage constraints (one per extreme SNP i):
  #   z_i - sum_{k: i in A_k} x_k <= 0
  #   i.e. z_i can be 1 only if at least one covering rule is selected.
  #
  # Optional cardinality constraint:
  #   sum_k x_k <= max_rules
  # -------------------------------------------------------------------
  n_vars <- A_n + K

  # Objective vector
  obj <- c(rep(-1.0, A_n), lambda * as.numeric(cand_bg))

  # Variable types (all binary)
  types <- rep("B", n_vars)

  # -------------------------------------------------------------------
  # Constraint matrix (sparse triplet)
  # -------------------------------------------------------------------

  # z-variable diagonal (coverage constraint rows, z_i column)
  z_rows <- seq_len(A_n)
  z_cols <- seq_len(A_n)
  z_vals <- rep(1.0, A_n)

  # x-variable entries: for each rule k, for each extreme SNP i in A_k,
  # row = i, col = A_n + k, value = -1
  lens   <- lengths(cand_A)
  x_rows <- unlist(cand_A,  use.names = FALSE)   # extreme SNP IDs (1..A_n)
  x_cols <- rep(seq_len(K) + A_n, times = lens)  # column = A_n + k
  x_vals <- rep(-1.0, length(x_rows))

  all_rows <- c(z_rows, x_rows)
  all_cols <- c(z_cols, x_cols)
  all_vals <- c(z_vals, x_vals)

  n_constr <- A_n
  rhs      <- rep(0.0, A_n)
  dir_vec  <- rep("<=", A_n)

  # Optional cardinality row: sum_k x_k <= max_rules
  if (!is.null(max_rules)) {
    n_constr  <- A_n + 1L
    card_cols <- seq_len(K) + A_n
    all_rows  <- c(all_rows, rep(n_constr, K))
    all_cols  <- c(all_cols, card_cols)
    all_vals  <- c(all_vals, rep(1.0, K))
    rhs       <- c(rhs, as.numeric(max_rules))
    dir_vec   <- c(dir_vec, "<=")
  }

  # Build slam sparse triplet matrix
  mat <- slam::simple_triplet_matrix(
    i    = all_rows,
    j    = all_cols,
    v    = all_vals,
    nrow = n_constr,
    ncol = n_vars
  )

  message(sprintf(
    "[%s] solving MILP: %d vars (%d z + %d x), %d constraints, lambda=%.3g%s",
    FUN, n_vars, A_n, K, n_constr, lambda,
    if (!is.null(max_rules)) sprintf(", max_rules=%d", max_rules) else ""
  ))

  # Solve
  sol <- Rglpk::Rglpk_solve_LP(
    obj    = obj,
    mat    = mat,
    dir    = dir_vec,
    rhs    = rhs,
    types  = types,
    max    = FALSE
  )

  if (sol$status != 0L) {
    warning(sprintf(
      "[%s] solver status %d (0 = optimal). Solution may be sub-optimal.",
      FUN, sol$status
    ))
  }

  # Extract binary solutions (round to nearest integer to handle near-zero/one)
  x_sol <- round(sol$solution[(A_n + 1L):n_vars])
  z_sol <- round(sol$solution[seq_len(A_n)])

  selected <- which(x_sol > 0.5)

  if (!length(selected)) {
    warning(sprintf("[%s] MILP selected no rules.", FUN))
    return(invisible(NULL))
  }

  n_covered_ext <- sum(z_sol > 0.5)

  # -------------------------------------------------------------------
  # Build a greedy-ordered ledger for the selected rules.
  # The MILP selects rules simultaneously, so there is no natural step
  # order. We simulate a greedy ordering over the selected subset
  # (always pick the rule with the highest marginal extreme coverage)
  # to produce a ledger that is directly comparable with
  # greedy_select_rules() output.
  # -------------------------------------------------------------------
  sel_A <- cand_A[selected]
  sel_B <- cand_B[selected]
  m     <- length(selected)

  covered_ext    <- logical(A_n)
  covered_bg     <- logical(B_n)
  remaining_pos  <- seq_len(m)
  order_seq      <- integer(m)

  for (step in seq_len(m)) {
    marg_exts <- vapply(remaining_pos, function(j) {
      if (length(sel_A[[j]])) sum(!covered_ext[sel_A[[j]]]) else 0L
    }, integer(1L))

    best_local         <- which.max(marg_exts)
    best_j             <- remaining_pos[best_local]
    order_seq[step]    <- best_j
    remaining_pos      <- remaining_pos[-best_local]

    if (length(sel_A[[best_j]])) covered_ext[sel_A[[best_j]]] <- TRUE
    if (length(sel_B[[best_j]])) covered_bg[sel_B[[best_j]]]  <- TRUE
  }

  # Rebuild coverage state for ledger (one pass in order_seq order)
  covered_ext2 <- logical(A_n)
  covered_bg2  <- logical(B_n)
  ledger_rows  <- vector("list", m)

  for (step in seq_len(m)) {
    j       <- order_seq[step]
    k_orig  <- selected[j]

    marg_ext <- if (length(sel_A[[j]])) sum(!covered_ext2[sel_A[[j]]]) else 0L
    marg_bg  <- if (length(sel_B[[j]])) sum(!covered_bg2[sel_B[[j]]])  else 0L

    if (length(sel_A[[j]])) covered_ext2[sel_A[[j]]] <- TRUE
    if (length(sel_B[[j]])) covered_bg2[sel_B[[j]]]  <- TRUE

    ledger_rows[[step]] <- data.table::data.table(
      step           = as.integer(step),
      rule_str       = rule_ids[k_orig],
      marginal_ext   = as.integer(marg_ext),
      marginal_bg    = as.integer(marg_bg),
      cumulative_ext = as.integer(sum(covered_ext2)),
      coverage_frac  = sum(covered_ext2) / A_n,
      cumulative_bg  = as.integer(sum(covered_bg2)),
      score          = NA_real_
    )
  }

  ledger <- data.table::rbindlist(ledger_rows)

  # Return selected rules in the same greedy simulation order as the ledger
  selected_ordered <- rule_ids[selected[order_seq]]

  message(sprintf(
    "[%s] optimal solution: %d rules selected, extreme coverage %.1f%% (%d / %d), status=%d.",
    FUN, m, 100 * n_covered_ext / A_n, n_covered_ext, A_n, sol$status
  ))

  list(
    selected_rules = selected_ordered,
    ledger         = ledger[],
    n_selected     = m,
    coverage_frac  = n_covered_ext / A_n,
    objective      = -sol$optimum,    # flip sign back to maximisation scale
    solver_status  = sol$status,
    A_n            = A_n,
    B_n            = B_n,
    lambda         = lambda
  )
}
