#' Drop entirely-NA columns from a data.table
#'
#' @param DT A \code{data.table}.
#'
#' @return The input \code{data.table} with any columns that are entirely
#'   \code{NA} removed. Columns containing at least one non-\code{NA} value
#'   are retained.
#'
#' @examples
.drop_all_na_cols <- function(DT) {
  keep <- !vapply(DT, function(x)
    all(is.na(x)), logical(1))
  DT[,names(DT)[keep], with = FALSE]
}

#' Jaccard similarity between two index vectors
#'
#' @param a Integer or numeric vector (duplicates are removed internally).
#' @param b Integer or numeric vector (duplicates are removed internally).
#'
#' @return Numeric scalar in \code{[0, 1]}: \code{|intersect(a, b)| /
#'   |union(a, b)|}. Returns \code{0} when both vectors are empty.
#' @keywords internal
.jacc <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  u <- length(union(a, b))
  if (u == 0L)
    return(0)
  length(intersect(a, b)) / u
}

