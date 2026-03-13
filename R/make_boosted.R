#' Assemble a `boosted` object from a trained `xgboost` model and labelled SNP index sets
#'
#' @description  Validates all model inputs, predicts SNP leaf assignments for test and training data, builds compact dense leaf maps for fast downstream tabulation, parses
#' the `xgboost` tree structure into a tidy \code{data.table}, and bundles
#' everything into a \code{boosted} object.
#'
#' The `boosted` object is the central data structure passed to all subsequent functions. An
#' optional \code{\link{set_snp_class}} step can be run beforehand to derive
#' the extreme and background index sets; alternatively, they can be supplied
#' directly.
#'
#' Typical usage does not require manipulating any of the `boosted` elements directly, but all internal functions (prefaced by a `.`) are documented by help pages and commented code for interested users.
#'
#' @param model A trained object of class \code{xgb.Booster}.
#' @param features_train Numeric matrix of feature values for the training
#'   SNPs (\eqn{n_\text{train} \times p}{}). Column names must match
#'   \code{model$feature_names} exactly.
#' @param features_test Numeric matrix of feature values for the test SNPs
#'   (\eqn{n_\text{test} \times p}{}). Must have identical column names to
#'   \code{features_train}.
#' @param yvar_train Numeric vector of length \eqn{n_\text{train}} containing
#'   the response values for training SNPs.
#' @param yvar_test Numeric vector of length \eqn{n_\text{test}} containing
#'   the response values for test SNPs.
#' @param extr_idx_train Integer vector of 1-based row indices into
#'   \code{yvar_train} identifying the extreme (positive-class) training SNPs.
#'   Must not overlap with \code{bg_idx_train}.
#' @param bg_idx_train Integer vector of 1-based row indices into
#'   \code{yvar_train} identifying the background (negative-class) training
#'   SNPs. Must not overlap with \code{extr_idx_train}.
#' @param extr_idx_test Integer vector of 1-based row indices into
#'   \code{yvar_test} identifying the extreme test SNPs. Must not overlap with
#'   \code{bg_idx_test}.
#' @param bg_idx_test Integer vector of 1-based row indices into
#'   \code{yvar_test} identifying the background test SNPs. Must not overlap
#'   with \code{extr_idx_test}.
#' @inheritParams .boosted_params
#'
#' @return An object of class \code{"boosted"}, a named list with the following
#'   elements:
#' \describe{
#'   \item{\code{yvar_train}, \code{yvar_test}}{The response vectors supplied
#'     by the caller, retained for use by sensitivity and rule analysis
#'     functions.}
#'   \item{\code{extr_idx_train}, \code{bg_idx_train}, \code{extr_idx_test},
#'     \code{bg_idx_test}}{Validated and sorted integer index vectors for the
#'     extreme and background sets in each partition.}
#'   \item{\code{N_extr_train}, \code{N_bg_train}, \code{N_extr_test},
#'     \code{N_bg_test}}{Integer scalars giving the cardinality of each
#'     labelled set.}
#'   \item{\code{N_index_train}, \code{N_index_test}}{Integer scalars: total
#'     number of labelled SNPs in the training and test partitions,
#'     respectively (\eqn{N_\text{extr} + N_\text{bg}}{}).}
#'   \item{\code{train_leaves}, \code{test_leaves}}{Integer matrices of
#'     dimensions \eqn{n \times T_m}{} containing the native xgboost leaf
#'     assignment for every SNP in every tree.}
#'   \item{\code{train_leaf_map}}{Named list produced by
#'     \code{.build_train_leaf_map}: \code{native_leaf_ids} (unique leaf IDs
#'     per tree), \code{dense_leaf_ids} (compact 1-based re-indexing for each
#'     SNP), and \code{n_leaves} (leaf count per tree).}
#'   \item{\code{test_leaf_map}}{Named list produced by
#'     \code{.build_test_leaf_map}: \code{dense_leaf_ids} mapping each test
#'     SNP to its position in the training leaf vocabulary (0 if unseen,
#'     though unseen leaves raise an error).}
#'   \item{\code{tdt}}{\code{data.table} with one row per \code{(Tree, node)}
#'     produced by \code{\link[boosted]{.parse_xgb_tree}}. See that function
#'     for the full column schema.}
#'   \item{\code{n_yvar_train}, \code{n_yvar_test}}{Integer scalars: total
#'     number of SNPs in each partition (labelled and unlabelled).}
#'   \item{\code{base_rate_train}, \code{base_rate_test}}{Numeric scalars:
#'     proportion of labelled SNPs that are extreme in each partition
#'     (\eqn{N_\text{extr} / (N_\text{extr} + N_\text{bg})}{}).}
#'   \item{\code{Tm}}{Integer scalar: number of trees in the model.}
#'   \item{\code{max_depth}}{Integer scalar: maximum number of split-steps
#'     from root to leaf observed across all trees, as inferred by
#'     \code{.infer_max_depth}.}
#' }
#' @export
#'
#' @examples
make_boosted <- function(model,
                         features_train,
                         features_test,
                         yvar_train,
                         yvar_test,
                         extr_idx_train,
                         bg_idx_train,
                         extr_idx_test,
                         bg_idx_test,
                         verbose = FALSE) {
  # Signature
  FUN <- "make_boosted"
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] start: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }

  # Validate xgboost model
  if (!inherits(model, "xgb.Booster"))
    stop(sprintf("[%s] Must provide a trained xgboost model.", FUN))

  # Validate feature matrices
  if (!is.matrix(features_train) || !is.matrix(features_test))
    stop(sprintf("[%s] features_train and features_test must be matrices.", FUN))
  if (!identical(colnames(features_train), colnames(features_test)))
    stop(sprintf(
      "[%s] Training and test feature sets have mismatched columns.",
      FUN
    ))
  model_feats <- model$feature_names
  if (!identical(model_feats, colnames(features_train)))
    stop(sprintf(
      "[%s] Feature names in model and feature matrices do not match.",
      FUN
    ))

  # Validate y-vars
  .validate_yvar(y        = yvar_train,
                 features = features_train,
                 caller   = FUN)
  .validate_yvar(y        = yvar_test,
                 features = features_test,
                 caller   = FUN)

  # Validate extreme and background index sets
  idxs <- .validate_indices(
    n_yvar_train   = length(yvar_train),
    extr_idx_train = extr_idx_train,
    bg_idx_train   = bg_idx_train,
    n_yvar_test    = length(yvar_test),
    extr_idx_test  = extr_idx_test,
    bg_idx_test    = bg_idx_test
  )

  # Predict leaf assignments
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] predicting leaves: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  features_train <- as.matrix(features_train)
  features_test  <- as.matrix(features_test)

  train_leaves <- predict(model, features_train, predleaf = TRUE)
  test_leaves  <- predict(model, features_test,  predleaf = TRUE)
  storage.mode(train_leaves) <- "integer"
  storage.mode(test_leaves)  <- "integer"

  # Get number of trees
  Tm <- ncol(train_leaves)

  # Precompute dense leaf maps for fast tabulation
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] building dense leaf maps: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  train_leaf_map <-
    .build_train_leaf_map(train_leaves  = train_leaves)
  test_leaf_map  <-
    .build_test_leaf_map(test_leaves    = test_leaves,
                         train_leaf_map = train_leaf_map)

  # Parse model into a numeric, per-node data.table
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] parsing tree dump: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  tdt <- .parse_xgb_tree(model = model)

  # Release large inputs we no longer need
  rm(model, features_train, features_test)
  gc()

  # Find the maximum realized tree depth (number of of splits)
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] finding maximum depth: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  max_depth <- .infer_max_depth(tdt)

  # Bundle
  boosted <- list(
    # y variables and indices
    yvar_train      = yvar_train,
    yvar_test       = yvar_test,
    extr_idx_train  = idxs$extr_idx_train,
    bg_idx_train    = idxs$bg_idx_train,
    extr_idx_test   = idxs$extr_idx_test,
    bg_idx_test     = idxs$bg_idx_test,
    N_extr_train    = idxs$N_extr_train,
    N_bg_train      = idxs$N_bg_train,
    N_extr_test     = idxs$N_extr_test,
    N_bg_test       = idxs$N_bg_test,
    N_index_train   = idxs$N_index_train,
    N_index_test    = idxs$N_index_test,

    # leaves and maps
    train_leaves    = train_leaves,
    test_leaves     = test_leaves,
    train_leaf_map  = train_leaf_map,
    test_leaf_map   = test_leaf_map,

    # model structure
    tdt             = tdt,

    # data sizes and metadata
    n_yvar_train    = length(yvar_train),
    n_yvar_test     = length(yvar_test),
    base_rate_train = idxs$N_extr_train / (idxs$N_extr_train + idxs$N_bg_train),
    base_rate_test  = idxs$N_extr_test / (idxs$N_extr_test + idxs$N_bg_test),
    Tm              = Tm,
    max_depth       = max_depth
  )
  class(boosted) <- "boosted"

  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] completed: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  boosted
}
