#' Precompute a sparse SNP-by-leaf membership matrix for permutation hot loops
#'
#' Builds a sparse binary matrix \code{A} of dimensions
#' \code{sum(Lvec) x n_train} from the training leaf map, where leaves from
#' all trees are stacked vertically in tree order. Inside the permutation loop
#' in \code{\link{.leaf_llrs_fast}}, one or more sparse matrix-vector products
#' replace repeated per-tree \code{tabulate()} calls. The pre-computation cost
#' is paid once before the loop begins.
#'
#' @param train_leaf_map Named list returned by
#'   \code{\link{.build_train_leaf_map}}, containing \code{dense_leaf_ids}
#'   (per-tree compact leaf assignments for all training SNPs) and
#'   \code{n_leaves} (leaf count per tree).
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{A}{Sparse binary matrix of dimensions
#'       \code{sum(Lvec) x n_train} (package \pkg{Matrix} class
#'       \code{dgCMatrix}). Row \code{i} corresponds to a single dense leaf
#'       in a single tree; column \code{j} corresponds to a training SNP.
#'       \code{A[i, j] = 1} if SNP \code{j} falls in leaf \code{i},
#'       0 otherwise. Leaves from all trees are stacked vertically in
#'       tree order.}
#'     \item{offsets}{Integer vector of length \code{Tm + 1}. The rows of
#'       \code{A} belonging to tree \code{t} are
#'       \code{(offsets[t] + 1):offsets[t + 1]}.}
#'     \item{Lvec}{Integer vector of length \code{Tm} giving the number of
#'       dense leaves in each tree. Copied from
#'       \code{train_leaf_map$n_leaves}.}
#'   }
#'
#' @import Matrix
#' @keywords internal

.build_leaf_matrix <- function(train_leaf_map) {
  dense_ids <- train_leaf_map$dense_leaf_ids
  Lvec      <- train_leaf_map$n_leaves
  Tm        <- length(dense_ids)
  n_train   <- length(dense_ids[[1]])

  # Row offset for each tree's block
  offsets <- c(0L, cumsum(Lvec))

  # Build triplet (row, col, value) for the full stacked matrix
  rows <- unlist(mapply(function(inv, offset) {
    inv + offset
  }, dense_ids, offsets[-length(offsets)], SIMPLIFY = FALSE))
  cols <- rep(seq_len(n_train), times = Tm)

  A <- sparseMatrix(
    i = rows,
    j = cols,
    x = 1,
    dims = c(offsets[Tm + 1], n_train)
  )

  list(
    A       = A,
    offsets = offsets,
    Lvec    = Lvec
  )
}
