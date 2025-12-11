#' Title
#'
#' @param test_leaves
#' @param train_leaf_map
#'
#' @return
#' @keywords internal
#'
#' @examples

.build_test_leaf_map <- function(test_leaves,
                                 train_leaf_map) {
  # number of trees
  Tm <- ncol(test_leaves)
  # Initialize container for results
  pos_list <- vector("list", Tm)

  # For each tree t in 1..Tm,
  for (t in seq_len(Tm)) {
    # Fetch unique native leaf IDs
    ids <- train_leaf_map$native_leaf_ids[[t]]
    # Get the positional index of each SNP's native leaf ID in ids
    pos <- match(as.integer(test_leaves[, t]), ids)

    # If any test leaf ID was not seen in train, fail loudly.
    if (any(is.na(pos))) {
      bad <- unique(as.integer(test_leaves[is.na(pos), t]))
      stop(
        sprintf(
          "[build_test_leaf_map] Tree %d contains leaf IDs in test not present in training leaves: %s",
          t,
          paste(bad, collapse = ", ")
        )
      )
    }
    pos_list[[t]] <- as.integer(pos)
  }

  # Return map
  list(dense_leaf_ids = pos_list)
}
