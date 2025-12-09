#' Title
#'
#' @param pairs
#' @param snps_ext_by_leaf
#' @param snps_bg_by_leaf
#' @param snps_all_by_leaf
#'
#' @return
#' @keywords internal
#'
#' @examples

.snp_lookup <- function(pairs,
                        snps_ext_by_leaf,
                        snps_bg_by_leaf,
                        snps_all_by_leaf = NULL) {

  # Initialize containers
  bucket_ext <- integer(0)
  bucket_bg  <- integer(0)
  bucket_all <- integer(0)

  if (!nrow(pairs)) {
    return(list(
      bucket_ext = bucket_ext,
      bucket_bg  = bucket_bg,
      bucket_all = bucket_all
    ))
  }

  for (j in seq_len(nrow(pairs))) {
    tr0 <- pairs$Tree[j] # 0-based tree index
    tt  <- tr0 + 1L      # convert to 1-based list index
    lid <- as.character(pairs$leaf_id[j]) # native leaf ID, used as list key

    # Look up SNP indices for this leaf in this tree
    re <- snps_ext_by_leaf[[tt]][[lid]]
    rb <- snps_bg_by_leaf [[tt]][[lid]]
    if (!is.null(re))
      bucket_ext <- c(bucket_ext, re)
    if (!is.null(rb))
      bucket_bg <- c(bucket_bg, rb)

    if (!is.null(snps_all_by_leaf)) {
      ra <- snps_all_by_leaf[[tt]][[lid]]
      if (!is.null(ra))
        bucket_all <- c(bucket_all, ra)
    }
  }

  # Deduplicate: a SNP can appear in multiple (Tree, leaf_id) pairs, but it only
  # counts once.
  bucket_ext <- unique(bucket_ext)
  bucket_bg  <- unique(bucket_bg)
  bucket_all <- unique(bucket_all)

  # Return the bucket of SNPs captured by the (Tree, leaf_id) pairs
  list(bucket_ext = bucket_ext,
       bucket_bg  = bucket_bg,
       bucket_all = bucket_all)
}
