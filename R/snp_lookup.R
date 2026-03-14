#' Collect unique SNP indices for a set of \code{(Tree, leaf_id)} pairs
#'
#' Given a table of \code{(Tree, leaf_id)} pairs, looks up the SNP indices
#' that fell in each leaf using the all-SNP inverse map and returns their
#' union (deduplicated) across all supplied pairs. Deduplication ensures that
#' a SNP counted by the same rule via multiple trees is not double-counted
#' when computing rule support. Currently only \code{snps_all_by_leaf} is
#' used in practice; downstream code (e.g. \code{\link{validate_rules}})
#' applies vector masks to \code{bucket_all} to separate extreme and
#' background SNPs rather than maintaining separate maps. The
#' \code{snps_ext_by_leaf} and \code{snps_bg_by_leaf} parameters are
#' retained for compatibility but are not used by any current call site.
#'
#' @param pairs A \code{data.frame} or \code{data.table} with at least two
#'   columns: \code{Tree} (0-based integer tree index, matching `xgboost`
#'   convention) and \code{leaf_id} (native integer leaf ID, matching the keys
#'   in the lookup maps).
#' @param snps_ext_by_leaf Legacy parameter, not currently used.
#'   List-of-lists mapping tree index (1-based) → native leaf ID (character)
#'   → integer vector of extreme SNP indices, or \code{NULL} (default) to
#'   skip.
#' @param snps_bg_by_leaf Legacy parameter, not currently used. As
#'   \code{snps_ext_by_leaf}, for background SNPs, or \code{NULL} (default)
#'   to skip.
#' @param snps_all_by_leaf List-of-lists mapping tree index (1-based) →
#'   native leaf ID (character) → integer vector of all SNP indices
#'   (labelled and unlabelled), as built by \code{\link{.build_lookup_R}} and
#'   stored in \code{boosted$snps_all_by_leaf_train} or
#'   \code{boosted$snps_all_by_leaf_test}. \code{NULL} to skip.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{bucket_all}}{Integer vector of unique SNP indices (labelled
#'     and unlabelled) covered by the supplied pairs. \code{integer(0)} if
#'     \code{snps_all_by_leaf} was \code{NULL} or no matching SNPs were
#'     found.}
#'   \item{\code{bucket_ext}}{Always \code{integer(0)} under current usage
#'     (legacy; see \code{snps_ext_by_leaf}).}
#'   \item{\code{bucket_bg}}{Always \code{integer(0)} under current usage
#'     (legacy; see \code{snps_bg_by_leaf}).}
#' }
#' @keywords internal

.snp_lookup_R <- function(pairs,
                          snps_ext_by_leaf = NULL,
                          snps_bg_by_leaf  = NULL,
                          snps_all_by_leaf = NULL) {

  bucket_ext <- integer(0)
  bucket_bg  <- integer(0)
  bucket_all <- integer(0)

  if (is.null(pairs) || !nrow(pairs)) {
    return(list(bucket_ext = bucket_ext,
                bucket_bg  = bucket_bg,
                bucket_all = bucket_all))
  }

  # TODO: snps_ext_by_leaf/snps_bg_by_leaf and their corresponding bucket
  # branches are legacy code from an earlier design that maintained three
  # separate lookup objects. All current call sites pass only snps_all_by_leaf
  # and use vector masks on bucket_all to separate extreme and background SNPs.
  # Consider simplifying the signature and return value to remove the dead
  # branches. Mirrors the equivalent dead code in .build_lookup_R.

  have_ext <- !is.null(snps_ext_by_leaf)
  have_bg  <- !is.null(snps_bg_by_leaf)
  have_all <- !is.null(snps_all_by_leaf)

  for (j in seq_len(nrow(pairs))) {
    tr0 <- pairs$Tree[j]                 # 0-based tree index (xgboost)
    tt  <- tr0 + 1L                      # 1-based list index (R)
    lid <- as.character(pairs$leaf_id[j]) # native leaf ID, used as list key

    if (have_ext) {
      re <- snps_ext_by_leaf[[tt]][[lid]]
      if (!is.null(re) && length(re)) bucket_ext <- c(bucket_ext, re)
    }
    if (have_bg) {
      rb <- snps_bg_by_leaf[[tt]][[lid]]
      if (!is.null(rb) && length(rb)) bucket_bg <- c(bucket_bg, rb)
    }
    if (have_all) {
      ra <- snps_all_by_leaf[[tt]][[lid]]
      if (!is.null(ra) && length(ra)) bucket_all <- c(bucket_all, ra)
    }
  }

  list(
    bucket_ext = unique(bucket_ext),
    bucket_bg  = unique(bucket_bg),
    bucket_all = unique(bucket_all)
  )
}

#' Internal dispatcher for \code{.snp_lookup_R}
#'
#' Calls the Rcpp implementation \code{.snp_lookup_rcpp} when available,
#' falling back to the pure-R \code{.snp_lookup_R}. See
#' \code{\link{.snp_lookup_R}} for full parameter and return documentation,
#' including notes on the legacy \code{snps_ext_by_leaf} and
#' \code{snps_bg_by_leaf} parameters.
#'
#' @keywords internal
.snp_lookup <- function(pairs,
                        snps_ext_by_leaf = NULL,
                        snps_bg_by_leaf  = NULL,
                        snps_all_by_leaf = NULL,
                        use_rcpp = TRUE) {

  if (use_rcpp &&
      exists(".snp_lookup_rcpp", mode = "function", inherits = TRUE)) {
    .snp_lookup_rcpp(pairs,
                     snps_ext_by_leaf,
                     snps_bg_by_leaf,
                     snps_all_by_leaf)
  } else {
    .snp_lookup_R(pairs,
                  snps_ext_by_leaf,
                  snps_bg_by_leaf,
                  snps_all_by_leaf)
  }
}
