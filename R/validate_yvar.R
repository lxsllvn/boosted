#' Check a y vector for compatibility with an xgboost model and return center and scale
#'
#' Checks if a y-variable vector is 1) numeric with finite center, 2) has finite and positive scale, and 3) is the length expected from the feature matrix. If valid, returns the mean, median, sd, and mad of the vector. Called by `make_boosted`.
#'
#' @param y vector of response values used either to train or test an xgboost model
#' @param features matrix of feature data
#' @param caller an optional string to customize error messages
#'
#' @return
#' @keywords internal
#'
#' @examples
.validate_yvar <- function(y,
                           features,
                           caller = ".validate_yvar") {
  if (!is.numeric(y))
    stop(sprintf("[%s] y must be numeric.", caller))
  if (!is.matrix(features))
    stop(sprintf("[%s] features must be a matrix.", caller))
  if (length(y) != nrow(features))
    stop(sprintf("[%s] length(y) != nrow(features).", caller))

  m    <- mean(y, na.rm = TRUE)
  md   <- stats::median(y, na.rm = TRUE)
  sdv  <- stats::sd(y, na.rm = TRUE)
  madv <- stats::mad(y, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(m))
    stop(sprintf("[%s] mean(y) is not finite.", caller))
  if (!is.finite(md))
    stop(sprintf("[%s] median(y) is not finite.", caller))
  if (!is.finite(sdv) || sdv <= 0)
    stop(sprintf("[%s] sd(y) is not finite/positive.", caller))
  if (!is.finite(madv) || madv <= 0)
    stop(sprintf("[%s] mad(y) is not finite/positive.", caller))
  invisible(TRUE)
}
