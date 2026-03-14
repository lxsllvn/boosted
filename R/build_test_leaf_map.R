#' Build a dense leaf map for the test partition
#'
#' Maps each test SNP's native `xgboost` leaf assignment to its position in the
#' training leaf vocabulary established by \code{\link{.build_train_leaf_map}}.
#' A test SNP whose leaf was observed in training receives a positive dense ID;
#' a leaf unseen in training raises an error immediately, since score lookup
#' would be undefined for such SNPs.
#'
#' @param test_leaves Integer matrix of dimensions
#'   \eqn{n_\text{test} \times T_m}{}, as returned by
#'   \code{predict(model, predleaf = TRUE)} on the test feature matrix.
#' @param train_leaf_map Named list returned by
#'   \code{\link{.build_train_leaf_map}}, used to look up the dense position
#'   of each test leaf within the training vocabulary.
#'
#' @return A named list with a single element:
#' \describe{
#'   \item{\code{dense_leaf_ids}}{List of length \code{Tm}. Element \code{t}
#'     is an integer vector of length \eqn{n_\text{test}}{} giving each test
#'     SNP's dense leaf position in tree \code{t} (i.e. its 1-based index
#'     into \code{train_leaf_map$native_leaf_ids[[t]]}). Raises an error if
#'     any test leaf ID was not present among the training leaves.}
#' }
#' @keywords internal

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
