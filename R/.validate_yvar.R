#' Title
#'
#' @param y
#' @param features
#' @param caller
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
