#' Prepare a boosted object for rule harvesting and overlap analysis
#'
#' Extends a \code{boosted} object with all the cached data structures needed
#' by \code{\link{validate_rules}}, \code{\link{analyze_rule_depth}}, and
#' \code{\link{analyze_rule_overlap}}. Specifically, this function: (1) builds
#' inverse leaf-to-SNP lookup tables for both the training and test partitions;
#' (2) constructs a feature binning specification from the split thresholds in
#' \code{boosted$tdt}, discretising continuous split values into stable
#' interpretable bins; (3) annotates \code{boosted$tdt} with the binned split
#' value (\code{Split_bin}); and (4) builds a rule cache mapping every
#' \code{(Tree, leaf_id)} to all of its prefixes up to \code{max_depth} splits,
#' with optional tightening of redundant monotone constraints. On completion,
#' the returned object gains class \code{"boosted_harvest"} and can be passed
#' directly to all rule-level functions.
#'
#' @param boosted A \code{boosted} object returned by
#'   \code{\link{make_boosted}}.
#' @param target_bins Positive integer. Target maximum number of bins per
#'   feature when constructing the binning specification.
#' @param min_per_bin Numeric scalar. Minimum bin mass (total split weight)
#'   required for a bin to be retained; bins below this threshold are merged
#'   with their neighbour.
#' @param winsor_prob Numeric scalar in \code{[0, 0.5)}. Tail probability used
#'   for winsorising split values before proposing bin boundaries, making the
#'   bins robust to extreme outlier splits.
#' @param method Character string, one of \code{"fd"} or \code{"quantile"}.
#'   Algorithm used to propose the initial bin boundaries: \code{"fd"} uses a
#'   weighted Freedman-Diaconis rule; \code{"quantile"} places boundaries at
#'   equally spaced weighted quantiles.
#' @param bin_weight Character string, one of \code{"node"}, \code{"leaf"}, or
#'   \code{"cover"}. How split thresholds are weighted when constructing bins:
#'   \code{"node"} weights each split equally; \code{"leaf"} weights by the
#'   number of descendant leaves; \code{"cover"} weights by the xgboost node
#'   Cover statistic (falls back to equal weighting if Cover is absent).
#' @param max_depth \code{NULL} or a positive integer. Maximum number of
#'   split-steps to include in a rule prefix. \code{NULL} (default) uses
#'   \code{boosted$max_depth}. Must not exceed the realized model depth.
#' @param tighten_monotone Logical. If \code{TRUE} (default), redundant
#'   monotone constraints on the same feature and direction within a rule
#'   prefix are collapsed to the tightest bound, shortening rule strings and
#'   making them easier to interpret.
#' @inheritParams .boosted_params
#'
#' @return The input \code{boosted} object with the following additional slots
#'   attached, and class \code{"boosted_harvest"} prepended:
#' \describe{
#'   \item{\code{harvest_bins}}{A named list with elements \code{breaks}
#'     (environment keyed by feature name, each holding a numeric breakpoint
#'     vector), \code{mids} (environment of bin midpoints), and \code{params}
#'     (the binning parameters used).}
#'   \item{\code{tdt}}{The original \code{tdt} \code{data.table}, now with an
#'     additional \code{Split_bin} column giving the binned split value for
#'     each internal node.}
#'   \item{\code{snps_all_by_leaf_train}, \code{snps_all_by_leaf_test}}{
#'     List-of-lists (tree → native leaf ID → integer SNP indices) mapping
#'     each leaf to all SNPs (labelled and unlabelled) that fell into it in
#'     the training and test partitions, respectively.}
#'   \item{\code{leaf_rule_cache}}{A \code{data.table} with one row per
#'     \code{(Tree, leaf_id, rule_len)} produced by \code{.build_rule_cache}.
#'     See that function for the full column schema.}
#'   \item{\code{tighten_monotone}}{Logical. The value of
#'     \code{tighten_monotone} used when building the cache.}
#'   \item{\code{harvest_max_depth}}{Integer. The effective \code{max_depth}
#'     used when building the cache.}
#' }
#' @export

prepare_harvest <- function(boosted,
                            target_bins = 10L,
                            min_per_bin = 50L,
                            winsor_prob = 0.01,
                            method  = c("fd", "quantile"),
                            bin_weight = c("node", "leaf", "cover"),
                            verbose = FALSE,
                            max_depth = NULL,
                            tighten_monotone = TRUE) {
  # Signature & basic checks
  FUN <- "prepare_harvest"
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] start: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }

  if (!inherits(boosted, "boosted")) {
    stop(sprintf(
      "[%s]Input must be an object of class 'boosted' (from make_boosted())",
      FUN
    ))
  }

  # Parse arguments
  method     <- match.arg(method)
  bin_weight <- match.arg(bin_weight)

  # Tabular model dump
  tdt <- boosted$tdt

  # Pull native and dense leaf IDs
  native_ids_all       <- boosted$train_leaf_map$native_leaf_ids
  dense_leaf_ids_train <- boosted$train_leaf_map$dense_leaf_ids
  dense_leaf_ids_test  <- boosted$test_leaf_map$dense_leaf_ids

  # If max_depth = NULL, pull max_depth from boosted. If provided, sanity-check
  # to make sure it's a single integer value no greater than the maximum depth
  # of the trees.
  if (is.null(max_depth)) {
    max_depth <- boosted$max_depth
  } else {
    if (length(max_depth) != 1L) {
      stop(sprintf("[%s] max_depth must be a single value.", FUN))
    }
    md_num <- suppressWarnings(as.numeric(max_depth))
    if (is.na(md_num) ||
        !is.finite(md_num) ||
        md_num < 1 || md_num != floor(md_num)) {
      stop(sprintf("[%s] max_depth must be a positive integer.", FUN))
    }
    max_depth <- as.integer(md_num)
    model_md <- boosted$max_depth
    if (max_depth > model_md) {
      stop(sprintf(
        "[%s] max_depth must be <= realized model max depth (%d).",
        FUN,
        model_md
      ))
    }
  }

  # -----------------------------------
  # Build inverse leaf maps (leaf → SNP list)
  # -----------------------------------
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] building SNP lookups: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }

  # For each tree tt and native leaf_id:
  # snps_all_by_leaf[[tt]][[leaf_id]] = all SNP indices (labeled + unlabeled)
  lookup_test <- .build_lookup(
    dense_leaf_ids  = dense_leaf_ids_test,
    native_leaf_ids = native_ids_all,
    all             = TRUE
  )

  lookup_train <- .build_lookup(
    dense_leaf_ids  = dense_leaf_ids_train,
    native_leaf_ids = native_ids_all,
    all             = TRUE
  )

  # -----------------------------------
  # Build feature binning specification
  # -----------------------------------
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] binning features: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }

  # Per-feature breakpoints and bin midpoints are stored in environments
  # keyed by feature name, then attached to boosted as harvest_bins.
  # Only internal nodes with finite splits
  splits <- tdt[!Leaf & is.finite(Split) & nzchar(Feature) & Feature != "Leaf"]
  feats  <- sort(unique(splits$Feature))

  # Environments keyed by feature
  breaks_env <- new.env(parent = emptyenv())
  mids_env   <- new.env(parent = emptyenv())

  # Weighted quantile helper (type=7-like interpolation)
  wquantile <- function(x, w, probs) {
    x <- as.numeric(x)
    w <- as.numeric(w)
    ok <- is.finite(x) & is.finite(w) & w > 0
    x <- x[ok]
    w <- w[ok]
    if (!length(x))
      return(rep(NA_real_, length(probs)))
    o <- order(x)
    x <- x[o]
    w <- w[o]
    cw <- cumsum(w)
    tot <- sum(w)
    if (!is.finite(tot) ||
        tot <= 0)
      return(rep(NA_real_, length(probs)))
    q <- numeric(length(probs))
    for (i in seq_along(probs)) {
      p <- probs[[i]]
      if (!is.finite(p)) {
        q[[i]] <- NA_real_
        next
      }
      if (p <= 0) {
        q[[i]] <- x[[1L]]
        next
      }
      if (p >= 1) {
        q[[i]] <- x[[length(x)]]
        next
      }
      tgt <- p * tot
      j <- which(cw >= tgt)[1L]
      q[[i]] <- x[[j]]
    }
    q
  }

  # If leaf-weighting, compute descendant leaf counts once per tree.
  if (identical(bin_weight, "leaf")) {
    tdt[, desc_leaves := as.integer(NA)]

    trees <- sort(unique(tdt$Tree))
    for (tt in trees) {
      dtt <- tdt[Tree == tt]
      if (!nrow(dtt))
        next

      max_id <- max(dtt$ID, na.rm = TRUE)
      row_of <- rep.int(NA_integer_, max_id + 1L)
      row_of[dtt$ID + 1L] <- seq_len(nrow(dtt))

      root_id <- 0L
      if (!(root_id %in% dtt$ID))
        root_id <- dtt$ID[which.min(dtt$ID)]

      # Postorder
      stack <- as.integer(root_id)
      order <- integer(0)
      seen  <- new.env(parent = emptyenv())

      while (length(stack)) {
        nid <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        key <- as.character(nid)
        if (!is.null(seen[[key]]))
          next
        seen[[key]] <- TRUE
        ridx <- row_of[nid + 1L]
        if (is.na(ridx))
          next
        order <- c(order, nid)
        if (!isTRUE(dtt$Leaf[[ridx]])) {
          y <- dtt$Yes[[ridx]]
          n <- dtt$No[[ridx]]
          if (is.finite(y))
            stack <- c(stack, as.integer(y))
          if (is.finite(n))
            stack <- c(stack, as.integer(n))
        }
      }

      order <- rev(order)
      leaf_count <- new.env(parent = emptyenv())

      for (nid in order) {
        ridx <- row_of[nid + 1L]
        if (is.na(ridx))
          next
        if (isTRUE(dtt$Leaf[[ridx]])) {
          leaf_count[[as.character(nid)]] <- 1L
        } else {
          y <- dtt$Yes[[ridx]]
          n <- dtt$No[[ridx]]
          ly <- 0L
          if (is.finite(y)) {
            tmp <- leaf_count[[as.character(y)]]
            if (!is.null(tmp))
              ly <- as.integer(tmp)
          }

          ln <- 0L
          if (is.finite(n)) {
            tmp <- leaf_count[[as.character(n)]]
            if (!is.null(tmp))
              ln <- as.integer(tmp)
          }

          leaf_count[[as.character(nid)]] <- as.integer(ly + ln)
        }
      }

      # Write back to tdt (Tree == tt) in-place
      keys <- ls(leaf_count, all.names = TRUE)
      if (length(keys)) {
        vals <- mget(keys, leaf_count, inherits = FALSE)
        map  <- as.integer(unlist(vals, use.names = FALSE))
        names(map) <- keys
        tdt[Tree == tt, desc_leaves := map[as.character(ID)]]
      }
    }
  }


  # Recompute internal split table after any leaf-weight precomputation
  splits <-
    tdt[!Leaf & is.finite(Split) & nzchar(Feature) & Feature != "Leaf"]
  feats  <- sort(unique(splits$Feature))

  # Split weights per internal node
  splits[, split_weight := {
    if (identical(bin_weight, "node")) {
      rep(1, .N)
    } else if (identical(bin_weight, "cover")) {
      if (!("Cover" %in% names(splits))) {
        rep(1, .N)
      } else {
        w <- Cover
        w[!is.finite(w) | w <= 0] <- 1
        w
      }
    } else {
      # leaf
      w <- as.numeric(desc_leaves)
      w[!is.finite(w) | w <= 0] <- 1
      w
    }
  }]

  # Drop scratch column if we created it.
  if ("desc_leaves" %in% names(tdt)) {
    tdt[, desc_leaves := NULL]
  }

  # Propose bin boundaries
  propose_breaks <- function(v, w, max_bins, method, winsor) {
    v <- as.numeric(v)
    w <- as.numeric(w)
    ok <- is.finite(v) & is.finite(w) & w > 0
    v <- v[ok]
    w <- w[ok]
    if (length(v) <= 1L) {
      r <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)))
        r <- c(0, 0)
      return(unique(c(r[1], r[2])))
    }

    # Weighted winsorization
    lo <- wquantile(v, w, winsor)
    hi <- wquantile(v, w, 1 - winsor)
    keep <- (v >= lo) & (v <= hi)
    v <- v[keep]
    w <- w[keep]

    if (length(v) <= 1L) {
      r <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)))
        r <- c(0, 0)
      return(unique(c(r[1], r[2])))
    }

    if (method == "quantile") {
      nb <- max(2L, as.integer(max_bins))
      probs <- seq(0, 1, length.out = nb + 1L)
      br <- unique(wquantile(v, w, probs))
      if (length(br) < 3L) {
        br <- pretty(v, n = min(max_bins, 3L))
      }
    } else {
      # FD rule using weighted IQR for robustness
      q25 <- wquantile(v, w, 0.25)
      q75 <- wquantile(v, w, 0.75)
      iqr <- q75 - q25
      h <- 2 * iqr / (length(v) ^ (1 / 3))
      if (!is.finite(h) || h <= 0) {
        h <- (max(v) - min(v)) / max(2L, max_bins)
      }
      nb <- ceiling((max(v) - min(v)) / max(h, .Machine$double.eps))
      nb <- min(max(nb, 2L), max_bins)
      br <- pretty(v, n = nb)
    }

    br <- sort(unique(as.numeric(br)))
    if (length(br) < 2L) {
      br <- unique(c(min(v), max(v)))
    }
    br
  }

  # Merge bins until constraints satisfied (using bin mass = sum(weights))
  enforce_min_bin <- function(v, w, br, target_bins, min_mass) {
    if (length(br) < 2L)
      return(br)

    repeat {
      nb <- length(br) - 1L
      idx <- findInterval(v, br, all.inside = TRUE)

      mass <- numeric(nb)
      # aggregate weights by bin
      agg <- tapply(w, idx, sum)
      if (length(agg)) {
        bins <- as.integer(names(agg))
        mass[bins] <- as.numeric(agg)
      }

      if ((nb <= target_bins) &&
          all(mass >= min_mass | mass == 0))
        break
      if (nb <= 1L)
        break

      k <- which.min(mass) # includes empty bins

      rm_pos <- if (k == 1L) {
        2L
      } else if (k == nb) {
        nb
      } else if (mass[k - 1L] <= mass[k + 1L]) {
        k
      } else {
        k + 1L
      }

      br <- br[-rm_pos]
      if (length(br) < 2L) {
        br <- sort(unique(br))
        break
      }
    }
    br
  }

  for (f in feats) {
    v <- splits[Feature == f, Split]
    w <- splits[Feature == f, split_weight]
    if (!length(v))
      next

    br <- propose_breaks(
      v        = v,
      w        = w,
      max_bins = target_bins,
      method   = method,
      winsor   = winsor_prob
    )

    br <- enforce_min_bin(
      v           = v,
      w           = w,
      br          = br,
      target_bins = target_bins,
      min_mass    = min_per_bin
    )

    if (length(br) < 2L) {
      r <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)))
        r <- c(0, 0)
      br <- unique(c(r[1], r[2]))
    }

    mids            <- (br[-1L] + br[-length(br)]) / 2
    breaks_env[[f]] <- br
    mids_env[[f]]   <- mids
  }

  # Attach binning spec
  boosted$harvest_bins <- list(
    breaks = breaks_env,
    mids   = mids_env,
    params = list(
      target_bins = as.integer(target_bins),
      min_per_bin = as.numeric(min_per_bin),
      winsor_prob = as.numeric(winsor_prob),
      method      = method,
      bin_weight  = bin_weight
    )
  )

  # Add bins to boosted$tdt
  if (!("Split_bin" %in% names(tdt))) {
    tdt[, Split_bin := as.numeric(NA_real_)]
  }
  tdt[!Leaf & is.finite(Split) & nzchar(Feature) & Feature != "Leaf",
      Split_bin := {
        f0 <- Feature[1L]
        br <- breaks_env[[f0]]
        md <- mids_env[[f0]]
        if (is.null(br) || is.null(md)) {
          Split
        } else {
          ii <- findInterval(Split, br, all.inside = TRUE)
          as.numeric(md[ii])
        }
      },
      by = Feature]

  boosted$tdt <- tdt

  # -----------------------------------
  # Attach SNP lookup and binning spec
  # -----------------------------------
  boosted$snps_all_by_leaf_train <- lookup_train$snps_all_by_leaf
  boosted$snps_all_by_leaf_test  <- lookup_test$snps_all_by_leaf

  # -----------------------------------
  # Cache per-leaf rule strings and prefix→leaf map
  # -----------------------------------
  if (isTRUE(verbose)) {
    message(sprintf("[%s] building rule cache: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  }

  cache <- .build_rule_cache(
    tdt              = tdt,
    max_depth        = max_depth,
    tighten_monotone = isTRUE(tighten_monotone)
  )

  boosted$leaf_rule_cache    <- cache
  boosted$tighten_monotone    <- isTRUE(tighten_monotone)
  boosted$harvest_max_depth   <- as.integer(max_depth)

  # Tag object as harvest-ready
  if (!inherits(boosted, "boosted_harvest")) {
    class(boosted) <- c("boosted_harvest", class(boosted))
  }
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] completed: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  boosted
}
