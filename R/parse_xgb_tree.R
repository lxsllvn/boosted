#' Parse an `xgboost` model tree dump into a tidy \code{data.table}
#'
#' A wrapper around \code{xgboost::xgb.model.dt.tree} that standardizes the
#' output into a consistent schema for downstream use by
#' \code{\link{make_boosted}} and \code{\link{prepare_harvest}}. Child-node
#' references are coerced to integer IDs; leaf prediction values and split
#' gain are extracted into dedicated columns; and a logical \code{Leaf} column
#' is added to distinguish decision nodes from leaf nodes. The function
#' handles variability in `xgboost` dump column names across versions (e.g.,
#' \code{Quality} vs. \code{Gain}) and fills missing columns with \code{NA}
#' so that the returned schema is always identical regardless of model version.
#'
#' @param model An object of class \code{xgb.Booster}.
#'
#' @return A \code{data.table} with one row per \code{(Tree, node)}, keyed
#'   and ordered by \code{(Tree, ID)}, with the following columns:
#' \describe{
#'   \item{\code{Tree}}{Integer. Zero-based index of the tree within the
#'     ensemble.}
#'   \item{\code{ID}}{Integer. Zero-based node index within the tree (the
#'     \code{Node} column from the raw `xgboost` dump, renamed for consistency
#'     with downstream code).}
#'   \item{\code{Feature}}{Character. For decision nodes, the name of the
#'     feature used at this split. For leaf nodes, the literal string
#'     \code{"Leaf"}.}
#'   \item{\code{Split}}{Numeric. For decision nodes, the raw split threshold
#'     on the original feature scale. \code{NA} for leaf nodes.}
#'   \item{\code{Yes}}{Integer. For decision nodes, the \code{ID} of the child
#'     node reached when the feature value is strictly less than \code{Split}.
#'     \code{NA} for leaf nodes.}
#'   \item{\code{No}}{Integer. For decision nodes, the \code{ID} of the child
#'     node reached when the feature value is greater than or equal to
#'     \code{Split}. \code{NA} for leaf nodes.}
#'   \item{\code{Missing}}{Integer. For decision nodes, the \code{ID} of the
#'     child node reached when the feature value is missing. \code{NA} for
#'     leaf nodes.}
#'   \item{\code{Leaf}}{Logical. \code{TRUE} for leaf nodes, \code{FALSE} for
#'     decision nodes.}
#'   \item{\code{LeafVal}}{Numeric. The raw `xgboost` leaf prediction value for
#'     leaf nodes; \code{NA} for decision nodes.}
#'   \item{\code{Gain}}{Numeric. The impurity reduction (loss function gain)
#'     achieved by the split at a decision node; \code{NA} for leaf nodes.}
#'   \item{\code{Cover}}{Numeric. The sum of instance weights seen by a
#'     decision node, or collected at a leaf node, as reported by `xgboost`.}
#' }
#' @keywords internal

.parse_xgb_tree <- function(model) {
  dt <- xgboost::xgb.model.dt.tree(model = model)
  data.table::setDT(dt)

  # Helper: extract integer child node IDs from strings like "yes=4" or "0-4".
  to_int_child <- function(x) {
    if (is.null(x))
      return(NA_integer_)
    y <- suppressWarnings(as.integer(sub(".*-", "", x)))
    y[is.na(x)] <- NA_integer_
    y
  }

  # Helper: safe numeric conversion with NA fallback.
  to_num <- function(x)
    suppressWarnings(as.numeric(x))

  # Standardize ID columns
  dt[, `:=`(Tree = as.integer(Tree),
            ID = as.integer(Node))]

  # Clean and standardize child pointers ("Yes", "No", "Missing").
  # If missing, we explicitly fill with NA_integer_ so the table
  # has a consistent schema for all models.
  for (col in c("Yes", "No", "Missing")) {
    if (col %in% names(dt))
      dt[, (col) := to_int_child(get(col))]
    else
      dt[, (col) := NA_integer_]
  }

  # Split value and Cover should be numeric. If missing, fill with NA.
  if (!"Split" %in% names(dt))
    dt[, Split := NA_real_]
  if (!"Cover" %in% names(dt))
    dt[, Cover := NA_real_]
  dt[, `:=`(Split = to_num(Split),
            Cover = to_num(Cover))]

  # Identify leaf nodes. XGBoost marks Feature="Leaf" at leaf entries.
  dt[, Leaf := (Feature == "Leaf")]

  # The gain column in tree dumps could be named either `Quality` or `Gain.`
  # Leaf values are given in one of the `Quality`,`Gain`, or `Split` columns,
  # with the others NA.
  has_quality <- "Quality" %in% names(dt)
  has_gaincol <- "Gain" %in% names(dt)

  Q <- if (has_quality)
    to_num(dt$Quality)
  else
    rep(NA_real_, nrow(dt))

  G <- if (has_gaincol)
    to_num(dt$Gain)
  else
    rep(NA_real_, nrow(dt))

  # LeafVal = leaf prediction. Prefer Quality if given; fall back to Split.
  dt[, LeafVal := ifelse(Leaf, Q, NA_real_)]
  dt[Leaf & is.na(LeafVal), LeafVal := Split]

  # Gain = impurity reduction. Prefer explicit Gain column; fall back to Quality.
  dt[, Gain := ifelse(!Leaf, ifelse(!is.na(G), G, Q), NA_real_)]

  # Drop columns that are either duplicated or irrelevant after processing.
  drop_cols <- intersect(c("Node", "Quality"), names(dt))
  if (length(drop_cols))
    dt[, (drop_cols) := NULL]

  # Final ordering for predictable downstream joins
  data.table::setkey(dt, Tree, ID)
  data.table::setorder(dt, Tree, ID)

  # Return canonical schema:
  # Tree, ID, Feature, Split, Yes, No, Missing, Leaf, LeafVal, Gain, Cover
  dt[, .(Tree,
         ID,
         Feature,
         Split,
         Yes,
         No,
         Missing,
         Leaf,
         LeafVal,
         Gain,
         Cover)]
}
