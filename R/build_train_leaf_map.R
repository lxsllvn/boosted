#' Build a dense leaf map for the training partition
#'
#' Converts a `xgboost` leaf assignment matrix (native integer leaf IDs,
#' one column per tree) into a compact representation that supports fast
#' tabulation in downstream scoring. For each tree, the set of unique native
#' leaf IDs is sorted and stored; each SNP's assignment is then re-expressed
#' as a 1-based position in that sorted set (the dense ID). This allows
#' downstream code to use \code{tabulate()} on a contiguous integer sequence
#' rather than on arbitrary native IDs.
#'
#' @param train_leaves Integer matrix of dimensions
#'   \eqn{n_\text{train} \times T_m}{}, as returned by
#'   \code{predict(model, predleaf = TRUE)}, with column \code{t} giving the
#'   native leaf ID assigned to each training SNP by tree \code{t}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{native_leaf_ids}}{List of length \code{Tm}. Element \code{t}
#'     is a sorted integer vector of the unique native leaf IDs observed in
#'     tree \code{t} among the training SNPs.}
#'   \item{\code{dense_leaf_ids}}{List of length \code{Tm}. Element \code{t}
#'     is an integer vector of length \eqn{n_\text{train}}{} giving each
#'     training SNP's position (1-based) within
#'     \code{native_leaf_ids[[t]]}.}
#'   \item{\code{n_leaves}}{Integer vector of length \code{Tm} giving the
#'     number of unique leaves in each tree.}
#' }
#' @keywords internal

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
