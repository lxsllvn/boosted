#' Title
#'
#' @param model
#' @param caller
#'
#' @return
#' @keywords internal
#'
#' @examples
.parse_xgb_tree <- function(model,
                            caller = ".parse_xgboost_tree") {
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
