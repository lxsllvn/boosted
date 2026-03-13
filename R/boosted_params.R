#' Shared parameter documentation for the boosted package
#'
#' This unexported function exists solely to hold canonical \code{@@param}
#' descriptions that are inherited by other functions via
#' \code{@@inheritParams .boosted_params}. It has no runtime behaviour.
#'
#' @param boosted A \code{boosted} object returned by \code{\link{make_boosted}},
#'   or a \code{boosted_harvest} object returned by \code{\link{prepare_harvest}}.
#' @param alpha Numeric scalar \eqn{\ge 0}. Jeffreys-prior concentration applied
#'   when smoothing leaf-level count ratios into log-likelihood ratios. When
#'   \code{alpha = 0} a small epsilon floor is used to avoid \code{log(0)};
#'   the default \code{alpha = 0.5} corresponds to the symmetric Jeffreys prior.
#' @param gain_grid Numeric vector of screening fractions in \code{(0, 1]} at
#'   which to evaluate the gain curve. Values are clipped to
#'   \code{[0.001, 0.999]} and deduplicated before use.
#' @param R Positive integer. Number of resampling or permutation iterations.
#' @param progress_every Optional positive integer. If supplied, a progress
#'   message is printed to the console after every \code{progress_every}
#'   iterations. \code{NULL} (default) suppresses progress messages.
#' @param verbose Logical. If \code{TRUE}, timestamped progress messages are
#'   printed at major internal checkpoints.
#' @param which Character string, one of \code{"train"} or \code{"test"}.
#'   Selects which data partition to use for SNP bucket assembly and scoring.
#' @param fold_indices Optional integer vector of SNP indices (1-based). If
#'   supplied, the analysis is restricted to these SNPs; extreme and background
#'   label sets are intersected with \code{fold_indices} before scoring. By
#'   default all SNPs in the selected partition are used.
#' @param candidate_rules Optional character vector of rule strings. Each
#'   string must be present in \code{boosted$leaf_rule_cache} (built by
#'   \code{\link{prepare_harvest}}). If \code{NULL}, the terminal (longest)
#'   rule of every leaf in the ensemble is used.
#' @param center Character string, one of \code{"mean"} or \code{"median"}.
#'   Measure of center used when standardising \code{yvar_train} to define
#'   extreme and background set boundaries.
#' @param scale Character string, one of \code{"sd"} or \code{"mad"}.
#'   Measure of spread used when standardising \code{yvar_train} to define
#'   extreme and background set boundaries.
#' @param lower_tail Logical. If \code{FALSE} (default), the extreme set is the
#'   upper tail (\eqn{\mu + k \cdot \sigma}); if \code{TRUE}, it is the lower
#'   tail (\eqn{\mu - k \cdot \sigma}).
#'
#' @keywords internal
.boosted_params <- function(boosted, alpha, gain_grid, R, progress_every,
                             verbose, which, fold_indices, candidate_rules,
                             center, scale, lower_tail) NULL
