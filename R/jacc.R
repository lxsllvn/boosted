#' Drop entirely NA columns from data table
#'
#' @param DT
#'
#' @return
#' @keywords internal
#'
#' @examples
.drop_all_na_cols <- function(DT) {
  keep <- !vapply(DT, function(x)
    all(is.na(x)), logical(1))
  DT[,names(DT)[keep], with = FALSE]
}


.jacc <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  u <- length(union(a, b))
  if (u == 0L)
    return(0)
  length(intersect(a, b)) / u
}

