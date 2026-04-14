#' Define extreme and background SNP class boundaries from training data
#'
#' Estimates a center and spread from \code{yvar_train} and computes fixed
#' thresholds that partition any y-variable vector into an extreme set (a
#' single tail beyond \code{extreme_k} standard units from center) and a
#' background set (a symmetric band within \code{bg_band_k} standard units of
#' center). Returns a named list containing the estimated parameters, the
#' derived thresholds, and a reusable \code{$apply} closure. Calling
#' \code{$apply} on the training and test y-variables separately ensures that
#' class boundaries are always derived from training data alone.
#'
#' @param yvar_train Numeric vector of response values for the training SNPs,
#'   used to estimate the center (\code{mu}) and spread (\code{sig}) that
#'   define the class boundaries.
#' @param extreme_k Positive numeric scalar. Distance from the training center
#'   in units of \code{sig} at which the extreme set boundary is placed.
#' @param bg_band_k Positive numeric scalar, must be less than
#'   \code{extreme_k}. Half-width of the symmetric background band around the
#'   training center, in units of \code{sig}. SNPs with
#'   \eqn{|\text{y} - \mu| \le \text{bg\_band\_k} \cdot \sigma}{} fall in the
#'   background set.
#' @inheritParams .boosted_params
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{mu}}{Numeric scalar. Estimated center of \code{yvar_train}
#'     (mean or median according to \code{center}).}
#'   \item{\code{sig}}{Numeric scalar. Estimated spread of \code{yvar_train}
#'     (SD or MAD according to \code{scale}).}
#'   \item{\code{extr_thr}}{Numeric scalar. Extreme set threshold on the
#'     original y-scale.}
#'   \item{\code{bg_low}, \code{bg_high}}{Numeric scalars. Lower and upper
#'     boundaries of the background band on the original y-scale.}
#'   \item{\code{apply}}{A function \code{function(y)} that accepts any
#'     numeric vector and returns a list with integer vectors
#'     \code{extr_idx} and \code{bg_idx} giving the 1-based positions of
#'     extreme and background SNPs, respectively, using the frozen thresholds.}
#' }
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

    bg_idx <- which(y > bg_low & y < bg_high)

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
