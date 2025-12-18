#' SNP lookup for rule buckets
#'
#' Given a data.table/data.frame of (Tree, leaf_id) pairs, return the unique SNP
#' indices captured by those leaves, optionally split into extreme/background/all
#' universes depending on which lookup maps are provided.
#'
#' @param pairs data.frame/data.table with columns Tree (0-based) and leaf_id (native)
#' @param snps_ext_by_leaf list-of-lists mapping tree -> leaf_id -> integer SNP indices (or NULL)
#' @param snps_bg_by_leaf  list-of-lists mapping tree -> leaf_id -> integer SNP indices (or NULL)
#' @param snps_all_by_leaf list-of-lists mapping tree -> leaf_id -> integer SNP indices (or NULL)
#'
#' @return list(bucket_ext, bucket_bg, bucket_all)
#' @keywords internal
#'
#' @examples
#'

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

#' Internal dispatcher for snp_lookup
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
