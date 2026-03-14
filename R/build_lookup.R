#' Build a leaf-to-SNP inverse index map
#'
#' For each tree, inverts the SNP-to-leaf assignment recorded in
#' \code{dense_leaf_ids} to produce a per-tree map from native leaf ID to the
#' SNP indices that fell in that leaf. Currently only the all-SNP map
#' (\code{snps_all_by_leaf}) is built in practice; downstream code uses
#' vector masks on \code{\link{.snp_lookup}} output to separate extreme and
#' background SNPs rather than maintaining separate maps. The \code{extr_idx}
#' and \code{bg_idx} parameters are retained for compatibility but are not
#' used by any current call site.
#'
#' @param dense_leaf_ids List of length \code{Tm}. Element \code{t} is an
#'   integer vector of length \eqn{n}{} (total SNPs in the relevant partition)
#'   giving each SNP's 1-based dense leaf index in tree \code{t}, as produced
#'   by \code{\link{.build_train_leaf_map}} or \code{\link{.build_test_leaf_map}}.
#' @param native_leaf_ids List of length \code{Tm}. Element \code{t} is an
#'   integer vector of length \eqn{L_t}{} mapping each dense leaf index to its
#'   native `xgboost` leaf ID, as stored in
#'   \code{train_leaf_map$native_leaf_ids}.
#' @param extr_idx Legacy parameter, not currently used. Integer vector of
#'   1-based SNP indices for the extreme set, or \code{NULL} (default) to skip
#'   building the extreme map.
#' @param bg_idx Legacy parameter, not currently used. Integer vector of
#'   1-based SNP indices for the background set, or \code{NULL} (default) to
#'   skip building the background map.
#' @param all Pass any non-\code{NULL} value (e.g. \code{TRUE}) to build the
#'   all-SNP map; \code{NULL} (default) skips it. Only the presence or absence
#'   of this argument is checked — the value itself is ignored.
#'
#' @return A named list with three elements. Elements for which the
#'   corresponding input was \code{NULL} are returned as \code{NULL}:
#' \describe{
#'   \item{\code{snps_ext_by_leaf}}{List of length \code{Tm}. Element \code{t}
#'     is a named list keyed by native leaf ID (as character), where each
#'     entry is an integer vector of the extreme SNP indices that fell in that
#'     leaf in tree \code{t}.}
#'   \item{\code{snps_bg_by_leaf}}{As above, for background SNPs.}
#'   \item{\code{snps_all_by_leaf}}{As above, for all SNPs (labelled and
#'     unlabelled).}
#' }
#' @keywords internal

.build_lookup_R <- function(dense_leaf_ids,
                            native_leaf_ids,
                            extr_idx = NULL,
                            bg_idx   = NULL,
                            all      = NULL) {

  Tm <- length(dense_leaf_ids)

  # TODO: extr_idx/bg_idx and their corresponding map branches are legacy code
  # from an earlier design that maintained three separate lookup objects. The
  # current approach builds only snps_all_by_leaf and uses vector masks in the
  # calling code (e.g. in validate_rules) to separate extreme and background
  # SNPs. Consider simplifying the signature to remove these parameters and the
  # dead branches below.

  have_ext <- !is.null(extr_idx) && length(extr_idx) > 0L
  have_bg  <- !is.null(bg_idx)   && length(bg_idx)   > 0L
  have_all <- !is.null(all)

  snps_ext_by_leaf <- if (have_ext) vector("list", Tm) else NULL
  snps_bg_by_leaf  <- if (have_bg)  vector("list", Tm) else NULL
  snps_all_by_leaf <- if (have_all) vector("list", Tm) else NULL

  for (tt in seq_len(Tm)) {
    inv_t <- dense_leaf_ids[[tt]]          # length = n (train or test universe)
    native_ids_t <- native_leaf_ids[[tt]]  # length = n_leaves_t (dense -> native)

    # For each SNP, what native leaf did it end up in?
    leaf_native <- native_ids_t[inv_t]     # length = n; values are native leaf IDs

    if (have_ext) {
      # Extreme SNPs under each native leaf ID
      snps_ext_by_leaf[[tt]] <- split(
        extr_idx,
        leaf_native[extr_idx],
        drop = TRUE
      )
    }

    if (have_bg) {
      # Background SNPs under each native leaf ID
      snps_bg_by_leaf[[tt]] <- split(
        bg_idx,
        leaf_native[bg_idx],
        drop = TRUE
      )
    }

    if (have_all) {
      # All SNPs (labeled + unlabeled) under each native leaf ID
      snps_all_by_leaf[[tt]] <- split(
        seq_along(leaf_native),
        leaf_native,
        drop = TRUE
      )
    }
  }

  list(
    snps_ext_by_leaf = snps_ext_by_leaf,
    snps_bg_by_leaf  = snps_bg_by_leaf,
    snps_all_by_leaf = snps_all_by_leaf
  )
}

#' Internal dispatcher for \code{.build_lookup_R}
#'
#' Calls the Rcpp implementation \code{.build_lookup_rcpp} when available,
#' falling back to the pure-R \code{.build_lookup_R}. The \code{all} argument
#' is coerced to a logical before being passed to the Rcpp layer, which
#' expects a scalar \code{SEXP} rather than an arbitrary sentinel value.
#' See \code{\link{.build_lookup_R}} for full parameter and return
#' documentation.
#'
#' @keywords internal
.build_lookup <- function(dense_leaf_ids,
                          native_leaf_ids,
                          extr_idx = NULL,
                          bg_idx   = NULL,
                          all      = NULL,
                          use_rcpp = TRUE) {

  if (use_rcpp &&
      exists(".build_lookup_rcpp", mode = "function", inherits = TRUE)) {

    # Rcpp expects SEXP for extr/bg now; pass through NULLs.
    .build_lookup_rcpp(
      dense_leaf_ids,
      native_leaf_ids,
      extr_idx,
      bg_idx,
      all = !is.null(all)  # historical: any non-NULL triggers building "all"
    )

  } else {

    .build_lookup_R(
      dense_leaf_ids  = dense_leaf_ids,
      native_leaf_ids = native_leaf_ids,
      extr_idx        = extr_idx,
      bg_idx          = bg_idx,
      all             = all
    )
  }
}
