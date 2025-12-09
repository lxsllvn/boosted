#' Title
#'
#' @param boosted
#' @param target_bins
#' @param min_per_bin
#' @param winsor_prob
#' @param method
#'
#' @return
#' @export
#'
#' @examples

prepare_harvest <- function(boosted,
                            target_bins = 10L,
                            min_per_bin = 50L,
                            winsor_prob = 0.01,
                            method = c("fd", "quantile")) {
  # Signature & basic checks
  FUN <- "prepare_harvest"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted")) {
    stop(sprintf(
      "[%s]Input must be an object of class 'boosted' (from make_boosted())",
      FUN
    ))
  }

  # `leaf_paths` are the raw per-leaf splits (feature, direction, split_val)
  # precomputed in make_boosted() via .extract_leaf_paths().
  leaf_paths <- boosted$leaf_paths
  method     <- match.arg(method)

  # Only use finite splits (drop NA/inf) and identify all features
  LS_split <- leaf_paths[is.finite(split_val)]
  feats    <- sort(unique(LS_split$feature))

  # Per-feature breakpoints and bin midpoints are stored in environments
  # keyed by feature name, then attached to boosted as harvest_bins.
  breaks_env <- new.env(parent = emptyenv())
  mids_env   <- new.env(parent = emptyenv())

  # Propose initial breakpoints for a single feature's split values
  propose_breaks <- function(v,
                             max_bins,
                             method,
                             winsor) {
    # Clean and sort values
    v <- sort(v[is.finite(v)])
    if (length(v) <= 1L) {
      # Degenerate: just fall back to [min, max] range
      r   <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)))
        r <- c(0, 0)
      return(unique(c(r[1], r[2])))
    }

    # Winsorize to trim extreme tails before binning
    lo <- suppressWarnings(stats::quantile(
      x     = v,
      probs = winsor,
      names = FALSE
    ))
    hi <- suppressWarnings(stats::quantile(
      x     = v,
      probs = 1 - winsor,
      names = FALSE
    ))
    v <- v[v >= lo & v <= hi]

    # If everything collapses after winsorization, again fall back to range
    if (length(v) <= 1L) {
      r <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)))
        r <- c(0, 0)
      return(unique(c(r[1], r[2])))
    }

    # Choose initial breaks either via Freedman–Diaconis or quantile-based bins
    if (method == "fd") {
      # FD rule for histogram width, with fallback if IQR is zero / degenerate
      h   <- 2 * stats::IQR(v) / (length(v) ^ (1 / 3))
      if (!is.finite(h) || h <= 0) {
        h <- (max(v) - min(v)) / max(2L, max_bins)
      }
      nb <- ceiling((max(v) - min(v)) / max(h, .Machine$double.eps))
      nb <- min(max(nb, 2L), max_bins)
      br <- pretty(v, n = nb)
    } else {
      # Quantile-based: aim for `max_bins` bins, then deduplicate / fall back
      nb    <- max(2L, as.integer(max_bins))
      probs <- seq(0, 1, length.out = nb + 1L)
      br    <- unique(stats::quantile(
        x     = v,
        probs = probs,
        names = FALSE,
        type  = 7
      ))
      if (length(br) < 3L) {
        br <- pretty(v, n = min(max_bins, 3L))
      }
    }

    # Clean up and ensure we have at least two distinct boundaries
    br <- sort(unique(as.numeric(br)))
    if (length(br) < 2L) {
      br <- unique(c(min(v), max(v)))
    }
    br
  }

  # Iteratively merge bins until:
  # - number of bins <= target_bins, and
  # - every non-empty bin has at least min_per_bin splits.
  enforce_min_bin <- function(v,
                              br,
                              target_bins,
                              min_per_bin) {
    if (length(br) < 2L)
      return(br)

    repeat {
      idx <- findInterval(v, br, all.inside = TRUE)
      tab <- tabulate(idx, nbins = length(br) - 1L)

      # Stop if: number of bins is acceptable AND all non-empty bins are large
      # enough
      if ((length(tab) <= target_bins) &&
          all(tab >= min_per_bin | tab == 0L)) {
        break
      }

      # Otherwise, remove one boundary adjacent to the smallest bin
      k <- if (length(tab))
        which.min(tab)
      else
        1L
      if (length(tab) <= 1L)
        break

      rm_pos <- if (k == 1L) {
        2L
      } else if (k == length(tab)) {
        length(tab)
      } else if (tab[k - 1L] <= tab[k + 1L]) {
        k
      } else {
        k + 1L
      }

      br <- br[-rm_pos]
      if (length(br) < 2L) {
        # Force at least two boundaries to avoid completely collapsing
        br <- br[1:2]
        break
      }
    }
    br
  }

  # Build per-feature breakpoints and midpoints across all trees
  for (f in feats) {
    # v = all split thresholds in the model for this feature
    v <- LS_split[feature == f, split_val]
    if (!length(v))
      next

    # Initial proposal for bin boundaries
    br <- propose_breaks(
      v        = v,
      max_bins = target_bins,
      method   = method,
      winsor   = winsor_prob
    )

    # Enforce minimum bin counts and cap on number of bins
    br <- enforce_min_bin(
      v           = v,
      br          = br,
      target_bins = target_bins,
      min_per_bin = min_per_bin
    )

    # Final safety check: always keep at least a min/max pair
    if (length(br) < 2L) {
      r <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)))
        r <- c(0, 0)
      br <- unique(c(r[1], r[2]))
    }

    # Store bin midpoints; these become the numeric thresholds used in rule
    # strings
    mids            <- (br[-1L] + br[-length(br)]) / 2
    breaks_env[[f]] <- br
    mids_env [[f]]  <- mids
  }

  # Attach binning spec to boosted object
  boosted$harvest_bins <- list(breaks = breaks_env,
                               mids   = mids_env)

  # Add second class tag so harvest_rules() can detect that binning is done
  if (!inherits(boosted, "boosted_binned")) {
    class(boosted) <- c("boosted_binned", class(boosted))
  }
  boosted
}
