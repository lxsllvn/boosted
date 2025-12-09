.jacc <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  u <- length(union(a, b))
  if (u == 0L)
    return(0)
  length(intersect(a, b)) / u
}
