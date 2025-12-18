#' Leaf → SNP lookup maps for fast bucket assembly
#'
#' Builds inverse indices for each tree: native leaf ID -> SNP indices.
#' You can request any combination of extreme/background/all maps by supplying
#' the corresponding inputs:
#'   - extr_idx NULL or length 0 => skip/ext map returned as NULL
#'   - bg_idx   NULL or length 0 => skip/bg map returned as NULL
#'   - all NULL => skip/all map returned as NULL
#'
#' NOTE: "all" follows historical semantics in boosted: any non-NULL value
#' triggers building snps_all_by_leaf (including all SNPs).
#'
#' @param dense_leaf_ids list<IntegerVector> per tree, length n (SNPs), values 1..L_t (dense leaf IDs)
#' @param native_leaf_ids list<IntegerVector> per tree, length L_t, mapping dense leaf ID -> native leaf ID
#' @param extr_idx integer vector of SNP indices (1-based) for extremes (or NULL)
#' @param bg_idx integer vector of SNP indices (1-based) for background (or NULL)
#' @param all NULL (skip) or any non-NULL value (build all)
#' @param use_rcpp logical; use Rcpp implementation if available
#'
#' @return list(snps_ext_by_leaf, snps_bg_by_leaf, snps_all_by_leaf) where omitted maps are NULL
#' @keywords internal
NULL

.build_lookup_R <- function(dense_leaf_ids,
                            native_leaf_ids,
                            extr_idx = NULL,
                            bg_idx   = NULL,
                            all      = NULL) {

  Tm <- length(dense_leaf_ids)

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

#' Internal dispatcher for build_lookup
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
