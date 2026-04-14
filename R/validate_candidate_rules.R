#' Resolve candidate rule strings to \code{(Tree, leaf_id)} pairs
#'
#' Validates and resolves a set of candidate rule strings against the rule
#' cache, returning the unique rule strings alongside a \code{data.table} of
#' every \code{(Tree, leaf_id)} pair in the ensemble where each rule occurs.
#' When \code{candidate_rules = NULL}, the terminal rule of each leaf (the
#' longest prefix stored for that leaf, up to \code{max_depth} splits) is used
#' as the candidate set. When an explicit character vector is supplied, each
#' string is checked for validity (non-empty, no leading/trailing whitespace,
#' no \code{NA}s) and matched against \code{leaf_rule_cache}; any strings not
#' found in the cache cause an informative error that names the missing rules
#' and hints at common causes such as mismatched \code{max_depth} or binning
#' parameters.
#'
#' @param leaf_rule_cache A \code{data.table} produced by
#'   \code{\link{.build_rule_cache}} and stored as
#'   \code{boosted$leaf_rule_cache} after \code{\link{prepare_harvest}}.
#' @param candidate_rules Optional character vector of rule strings to
#'   evaluate. Each string must be present in \code{leaf_rule_cache$rule_str}.
#'   If \code{NULL} (default), the terminal rule of every leaf in the ensemble
#'   is used.
#' @param caller Character scalar: name of the calling function, included in
#'   error messages.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{candidate_rules}}{Character vector of deduplicated rule
#'     strings, one element per unique rule that will be evaluated.}
#'   \item{\code{pairs_all}}{A \code{data.table} with one row per
#'     \code{(rule_str, Tree, leaf_id)} combination, containing all rows
#'     from \code{leaf_rule_cache} that match any element of
#'     \code{candidate_rules}. Columns: \code{rule_str}, \code{Tree},
#'     \code{leaf_id}, \code{rule_len}, \code{n_clauses}, \code{leaf_depth}.
#'     An index on \code{rule_str} is set for fast downstream joins.}
#' }
#' @keywords internal
.validate_candidate_rules <- function(leaf_rule_cache,
                                      candidate_rules = NULL,
                                      caller = ".validate_candidate_rules") {
  if (is.null(candidate_rules)) {
    # If no rules provided, use the terminal rule strings from the model.
    candidate_rules <-
      leaf_rule_cache[, .SD[which.max(rule_len)], by = .(Tree, leaf_id)]$rule_str
    uniq_rules <- unique(candidate_rules)

    # Map rule_str -> (Tree, leaf_id) pairs using the terminal rule cache.
    pairs_all <- leaf_rule_cache[data.table::data.table(rule_str = uniq_rules),
                            on = "rule_str",
                            nomatch = 0L,
                            .(rule_str, Tree, leaf_id, rule_len, n_clauses, leaf_depth)]
  } else {
    if (!is.character(candidate_rules)) {
      stop(sprintf("[%s] %s must be a character vector.", caller))
    }
    if (anyNA(candidate_rules)) {
      stop(sprintf("[%s] candidate_rules contain NA.", caller))
    }
    if (any(!nzchar(candidate_rules))) {
      stop(sprintf("[%s] candidate_rules contain empty strings.", caller))
    }
    if (any(candidate_rules != trimws(candidate_rules))) {
      stop(sprintf(
        "[%s] candidate_rules contains leading/trailing whitespace.",
        caller
      ))
    }
    uniq_rules <- unique(candidate_rules)
    pairs_all <- leaf_rule_cache[data.table::data.table(rule_str = uniq_rules),
                            on = "rule_str",
                            nomatch = 0L,
                            .(rule_str, Tree, leaf_id, rule_len, n_clauses, leaf_depth)]

    missing <- setdiff(uniq_rules, unique(pairs_all$rule_str))

    if (length(missing)) {
      ex <- paste(utils::head(missing, 10L), collapse = "\n  ")
      stop(
        sprintf(
          "[%s] %d/%d candidate rules not present in leaf cache
            (maybe different max_depth/binning/tighten_monotone?).\n  %s",
          caller,
          length(missing),
          length(uniq_rules),
          ex
        )
      )
    }
  }
  list(candidate_rules = uniq_rules,
       pairs_all       = pairs_all[])

}
