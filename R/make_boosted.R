#' Title
#'
#' @param model
#' @param features_train
#' @param features_test
#' @param yvar_train
#' @param yvar_test
#' @param extr_idx_train
#' @param bg_idx_train
#' @param extr_idx_test
#' @param bg_idx_test
#' @param verbose
#'
#' @return
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

  # Construct parent/child maps and depth/path lookup functions for tree
  # navigation
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] building helpers: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  helpers <- .build_tree_helpers(tdt = tdt)

  # Generate table of per-leaf decision paths
  if (isTRUE(verbose)) {
    message(sprintf(
      "[%s] building path table: %s",
      FUN,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  }
  leaf_paths <-
    .extract_leaf_paths(
      native_leaf_ids = train_leaf_map$native_leaf_ids,
      Tm              = Tm,
      tdt             = tdt,
      helpers         = helpers,
      # trees_per_batch limits RAM usage
      trees_per_batch = 250L
    )

  #Find number of splits along the deepest path
  max_depth <- as.integer(max(leaf_paths$depth) + 1L) # zero-based to 1-based

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

    # model structure and helpers
    tdt             = tdt,
    helpers         = helpers,
    leaf_paths      = leaf_paths,

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
