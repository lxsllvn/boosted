#' Validate a y-variable vector for compatibility with a feature matrix
#'
#' @description
#' Checks that \code{y} is numeric with a finite, positive-scale distribution
#' and has the same length as \code{nrow(features)}. Stops with an informative
#' error if any check fails. Called by \code{\link{make_boosted}} for both the
#' training and test y-variables before any leaf predictions are made.
#'
#' @param y Numeric vector of response values (training or test).
#' @param features Numeric matrix of feature values; \code{nrow(features)} must equal \code{length(y)}.
#' @param caller an optional string to customize error messages
#' @param caller Character scalar: name of the calling function, included in error messages. Defaults to \code{".validate_yvar"}.
#'
#' @return \code{invisible(TRUE)}, called for its side-effect of stopping on
#'   invalid input.
#' @keywords internal

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
