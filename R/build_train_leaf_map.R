#' Title
#'
#' @param train_leaves
#'
#' @return
#' @keywords internal
#'
#' @examples

.build_train_leaf_map <- function(train_leaves) {
  # number of trees
  Tm <- ncol(train_leaves)

  # Initialize containers for results
  # native IDs
  native_leaf_ids <- vector("list", Tm)
  # 1..length(native_leaf_ids[[t]]) dense mapping
  dense_leaf_ids <- vector("list", Tm)
  # number of unique leaves per tree
  n_leaves <- integer(Tm)

  # For each tree t in 1..Tm,
  for (t in seq_len(Tm)) {
    # Fetch native assignments for each SNP
    lt  <- as.integer(train_leaves[, t])
    # vector of unique leaf IDs
    ids <- sort(unique(lt))
    # Stash unique native leaf IDs
    native_leaf_ids[[t]] <- ids
    # Get the positional index of each SNP's native leaf ID in ids
    dense_leaf_ids[[t]] <- match(lt, ids)
    # number of leaves in this tree
    n_leaves[t] <- length(ids)
  }

  # Return map
  list(
    # unique native IDs per tree
    native_leaf_ids = native_leaf_ids,
    # dense 1...length(ids) leaf IDs for each SNP in each tree
    dense_leaf_ids = dense_leaf_ids,
    # number of unique leaves per tree
    n_leaves = n_leaves
  )
}
