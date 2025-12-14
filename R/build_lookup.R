#' Title
#'
#' @param dense_leaf_ids
#' @param native_leaf_ids
#' @param extr_idx
#' @param bg_idx
#' @param all
#'
#' @return
#' @keywords internal
#'
#' @examples

.build_lookup <- function(dense_leaf_ids,
                          native_leaf_ids,
                          extr_idx,
                          bg_idx,
                          all = NULL) {

  Tm <- length(dense_leaf_ids)

  snps_ext_by_leaf <- vector("list", Tm)
  snps_bg_by_leaf  <- vector("list", Tm)

  have_all <- !is.null(all)
  snps_all_by_leaf <- if (have_all)
    vector("list", Tm)
  else
    NULL

  for (tt in seq_len(Tm)) {
    # inv_t: dense leaf ID for each SNP
    inv_t <- dense_leaf_ids[[tt]]          # length = n_train (or n_test)
    native_ids_t <- native_leaf_ids[[tt]]  # length = n_leaves_t

    # For each SNP, what native leaf did it end up in?
    leaf_native <- native_ids_t[inv_t]     # row -> native leaf id

    # Extreme SNPs under each native leaf ID
    snps_ext_by_leaf[[tt]] <- split(
      extr_idx,
      leaf_native[extr_idx],
      drop = TRUE
    )

    # Background SNPs under each native leaf ID
    snps_bg_by_leaf[[tt]] <- split(
      bg_idx,
      leaf_native[bg_idx],
      drop = TRUE
    )

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
                          extr_idx,
                          bg_idx,
                          all = NULL,
                          use_rcpp = TRUE) {

  if (use_rcpp &&
      exists(".build_lookup_rcpp", mode = "function", inherits = TRUE)) {

    .build_lookup_rcpp(
      dense_leaf_ids,
      native_leaf_ids,
      extr_idx,
      bg_idx,
      all = !is.null(all)   # bool for C++
    )

  } else {

    .build_lookup_R(
      dense_leaf_ids,
      native_leaf_ids,
      extr_idx,
      bg_idx,
      all = all             # NULL/TRUE passthrough to R version
    )
  }
}

