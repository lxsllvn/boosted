#' Title
#'
#' @param yvar_train
#' @param extreme_k
#' @param lower_tail
#' @param bg_band_k
#' @param center
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
set_snp_class <- function(yvar_train,
                          extreme_k,
                          lower_tail = FALSE,
                          bg_band_k,
                          center = c("mean", "median"),
                          scale  = c("sd", "mad")) {
  # Parse input arguments
  center <- match.arg(center)
  scale  <- match.arg(scale)

  # Compute center and scale of the training data
  mu <- switch(
    center,
    mean   = mean(yvar_train, na.rm = TRUE),
    median = median(yvar_train, na.rm = TRUE)
  )

  sig <- switch(
    scale,
    sd  = stats::sd(yvar_train, na.rm = TRUE),
    mad = stats::mad(yvar_train, constant = 1.4826, na.rm = TRUE)
  )

  # Check for pathological scale
  if (!is.finite(sig) || sig <= 0) {
    stop("set_snp_class(): scale estimate is zero or non-finite.")
  }

  # Define fixed thresholds for background set
  # Background is symmetric: mu ± bg_band_k * sigma
  bg_low  <- mu - bg_band_k * sig
  bg_high <- mu + bg_band_k * sig

  # Define fixed threshold for extreme set (upper or lower
  # depending on lower_tail)
  if (lower_tail) {
    extr_thr <- mu - extreme_k * sig
  } else {
    extr_thr <- mu + extreme_k * sig
  }

  # Construct function to create sets
  apply_fun <- function(y) {
    if (lower_tail) {
      extr_idx <- which(y <= extr_thr)
    } else {
      extr_idx <- which(y >= extr_thr)
    }

    bg_idx <- which(y >= bg_low & y <= bg_high)

    list(extr_idx = extr_idx,
         bg_idx   = bg_idx)
  }

  # Return function and the value of center, scale, and thresholds
  list(
    mu        = mu,
    sig       = sig,
    extr_thr  = extr_thr,
    bg_low    = bg_low,
    bg_high   = bg_high,
    apply     = apply_fun
  )
}
