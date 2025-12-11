#' Title
#'
#' @param n_yvar_train
#' @param extr_idx_train
#' @param bg_idx_train
#' @param n_yvar_test
#' @param extr_idx_test
#' @param bg_idx_test
#' @param caller
#'
#' @return
#' @keywords internal
#'
#' @examples

.validate_indices <- function(n_yvar_train,
                              extr_idx_train,
                              bg_idx_train,
                              n_yvar_test,
                              extr_idx_test,
                              bg_idx_test,
                              caller = ".validate_indices") {
  check_idx <- function(idx, n, label) {
    if (is.null(idx))
      stop(sprintf("[%s] %s is NULL.", caller, label))
    if (length(idx) == 0L)
      stop(sprintf("[%s] %s is empty.", caller, label))
    if (anyNA(idx))
      stop(sprintf("[%s] %s contains NA.", caller, label))
    if (!is.numeric(idx))
      stop(sprintf("[%s] %s must be numeric/integerish.", caller, label))
    if (any(idx %% 1 != 0))
      stop(sprintf("[%s] %s has non-integer values.", caller, label))
    idx <- as.integer(idx)
    if (any(idx < 1L | idx > n))
      stop(sprintf("[%s] %s out of range (1..%d).", caller, label, n))
    if (anyDuplicated(idx))
      stop(sprintf("[%s] %s has duplicates.", caller, label))
    sort(idx)
  }

  # Validate training index sets
  ei_tr <- check_idx(extr_idx_train, n_yvar_train, "extr_idx_train")
  bi_tr <- check_idx(bg_idx_train,   n_yvar_train, "bg_idx_train")
  if (length(intersect(ei_tr, bi_tr)))
    stop(sprintf("[%s] extr_idx_train and bg_idx_train overlap.", caller))

  # Validate test index sets
  ei_te <- check_idx(extr_idx_test, n_yvar_test, "extr_idx_test")
  bi_te <- check_idx(bg_idx_test,   n_yvar_test, "bg_idx_test")
  if (length(intersect(ei_te, bi_te)))
    stop(sprintf("[%s] extr_idx_test and bg_idx_test overlap.", caller))

  # Return SNP indices and set sizes for training and test data
  list(
    extr_idx_train = ei_tr,
    bg_idx_train   = bi_tr,
    extr_idx_test  = ei_te,
    bg_idx_test    = bi_te,
    N_extr_train   = length(ei_tr),
    N_bg_train     = length(bi_tr),
    N_extr_test    = length(ei_te),
    N_bg_test      = length(bi_te),
    N_index_train  = length(c(ei_tr, bi_tr)),
    N_index_test   = length(c(ei_te, ei_te))
  )
}
