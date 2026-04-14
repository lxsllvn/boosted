#' Prepare a boosted object for rule harvesting and overlap analysis
#'
#' Extends a \code{boosted} object with all the cached data structures needed
#' by \code{\link{validate_rules}}, \code{\link{analyze_rule_depth}}, and
#' \code{\link{analyze_rule_overlap}}. Specifically, this function: (1) builds
#' inverse leaf-to-SNP lookup tables for both the training and test partitions;
#' (2) builds a rule cache mapping every \code{(Tree, leaf_id)} to all of its
#' prefixes up to \code{max_depth} splits, with optional tightening of
#' redundant monotone constraints. Rule strings use raw xgboost split
#' thresholds formatted to 4 significant figures in scientific notation. On
#' completion, the returned object gains class \code{"boosted_harvest"} and can
#' be passed directly to all rule-level functions.
#'
#' @param boosted A \code{boosted} object returned by
#'   \code{\link{make_boosted}}.
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
#'   \item{\code{tdt}}{The original \code{tdt} \code{data.table}, unchanged.}
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
      "[%s] Input must be an object of class 'boosted' (from make_boosted())",
      FUN
    ))
  }

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

  boosted$tdt <- tdt

  # -----------------------------------
  # Attach SNP lookup
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
