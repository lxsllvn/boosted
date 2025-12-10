#' Title
#'
#' @param boosted
#' @param harvest
#' @param top_n
#' @param progress_every
#'
#' @return
#' @export
#'
#' @examples

analyze_rule_overlap <- function(boosted,
                                 harvest,
                                 top_n = 500L,
                                 progress_every = NULL) {
  # Signature & basic checks
  FUN <- "analyze_rule_overlap"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted")) {
    stop(
      sprintf(
        "[%s]  Input must be an object of class 'boosted' (from make_boosted()).",
        FUN
      )
    )
  }
  if (!inherits(harvest, "boosted_harvest")) {
    stop(
      sprintf(
        "[%s] harvest must be the object returned by harvest/validate_rules(..., return_ledger = TRUE).",
        FUN
      )
    )
  }

  # Pull data from boosted
  extr_idx_train <- boosted$extr_idx_train
  bg_idx_train   <- boosted$bg_idx_train
  Tm             <- boosted$Tm

  train_leaf_map <- boosted$train_leaf_map
  dense_leaf_ids <- train_leaf_map$dense_leaf_ids
  native_ids_all <- train_leaf_map$native_leaf_ids

  # Harvest ledger pieces: rule table + mapping from rules to (Tree, leaf_id)
  R_tbl     <- data.table::as.data.table(harvest$R)
  pairs_all <- data.table::as.data.table(harvest$pairs_all)

  # Order rules exactly as in harvest_rules() so everything lines up
  ord_cols <-
    intersect(c("lift", "odds_ratio", "recall", "precision"),
              names(R_tbl))
  if (length(ord_cols)) {
    data.table::setorderv(R_tbl,
                          cols  = ord_cols,
                          order = rep.int(-1L, length(ord_cols)))
  }

  # Restrict to the top_n rules (if requested)
  if (is.finite(top_n) && top_n > 0L && top_n < nrow(R_tbl)) {
    R_use <- R_tbl[seq_len(top_n)]
  } else {
    R_use <- R_tbl
  }

  K <- nrow(R_use)
  if (K < 2L) {
    warning(sprintf("[%s] fewer than two rules selected; nothing to compare.", FUN))
    return(invisible(NULL))
  }

  rule_ids <- R_use$rule_str

  # Build per-tree leaf → SNP maps, reused for all rules and prefixes.
  # For each tree tt and native leaf_id:
  # snps_ext_by_leaf[[tt]][[leaf_id]] = extreme SNP indices
  # snps_bg_by_leaf [[tt]][[leaf_id]] = background SNP indices
  # snps_all_by_leaf[[tt]][[leaf_id]] = all SNP indices (labeled + unlabeled)
  lookup <- .build_lookup(
    dense_leaf_ids  = dense_leaf_ids,
    native_leaf_ids = native_ids_all,
    extr_idx        = extr_idx_train,
    bg_idx          = bg_idx_train,
    all             = TRUE
  )
  snps_ext_by_leaf <- lookup$snps_ext_by_leaf
  snps_bg_by_leaf  <- lookup$snps_bg_by_leaf
  snps_all_by_leaf <- lookup$snps_all_by_leaf

  # For each rule, reconstruct which SNPs it captures on the training set
  buckets_ext <- vector("list", K)
  buckets_bg  <- vector("list", K)
  buckets_all <- vector("list", K)

  for (k in seq_len(K)) {
    rs <- rule_ids[k]

    # All (Tree, leaf_id) pairs where this rule held in the ensemble
    pr <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
    if (!nrow(pr)) {
      buckets_ext[[k]] <- integer(0)
      buckets_bg [[k]] <- integer(0)
      buckets_all[[k]] <- integer(0)
      next
    }

    # Look up SNP indices covered by these (Tree, leaf_id) pairs.
    # .snp_lookup() handles deduplication (no double-counting SNPs
    # that land in the same rule via multiple trees).
    buckets <- .snp_lookup(
      pairs = pr,
      snps_ext_by_leaf = snps_ext_by_leaf,
      snps_bg_by_leaf  = snps_bg_by_leaf,
      snps_all_by_leaf = snps_all_by_leaf
    )
    buckets_ext[[k]] <- buckets$bucket_ext
    buckets_bg [[k]] <- buckets$bucket_bg
    buckets_all[[k]] <- buckets$bucket_all
  }

  # Optional progress over rule buckets
  if (!is.null(progress_every) && progress_every > 0L &&
      (k %% progress_every == 0L)) {
    message(sprintf("[%s] built SNP buckets for %d / %d rules",
                    FUN, k, K))
  }
}
# Rule-level summary: just reuse the harvest output
summary_cols <- intersect(
  c(
    "rule_str",
    "rule_len",
    "n_extreme",
    "n_bg",
    "support",
    "support_all",
    "precision",
    "recall",
    "lift",
    "odds_ratio",
    "enrichment",
    # if present
    "med_y_extreme",
    "med_y_bg",
    "med_y_overall"
  ),
  names(R_use)
)
summary_tbl <- R_use[, ..summary_cols]

# Initialize Jaccard matrices over extreme / background / all SNPs
labs <- paste0("r", seq_len(K))
J_ext <-
  matrix(
    NA_real_,
    nrow = K,
    ncol = K,
    dimnames = list(labs, labs)
  )
J_bg <-
  matrix(
    NA_real_,
    nrow = K,
    ncol = K,
    dimnames = list(labs, labs)
  )
J_all <-
  matrix(
    NA_real_,
    nrow = K,
    ncol = K,
    dimnames = list(labs, labs)
  )

# Also build a rich long-form overlap table:
# for each (i,j) we store intersection sizes and directional proportions.
overlap_list <- vector("list", K * (K - 1L) / 2L)
oi <- 1L

for (i in seq_len(K)) {
  Ai_ext <- buckets_ext[[i]]
  Ai_bg  <- buckets_bg [[i]]
  Ai_all <- buckets_all[[i]]

  nA_ext <- length(Ai_ext)
  nA_bg  <- length(Ai_bg)
  nA_all <- length(Ai_all)

  for (j in seq_len(K)) {
    Aj_ext <- buckets_ext[[j]]
    Aj_bg  <- buckets_bg [[j]]
    Aj_all <- buckets_all[[j]]

    nB_ext <- length(Aj_ext)
    nB_bg  <- length(Aj_bg)
    nB_all <- length(Aj_all)

    # Jaccard similarities for extremes, background, and all SNPs
    J_ext[i, j] <- .jacc(Ai_ext, Aj_ext)
    J_bg [i, j] <- .jacc(Ai_bg,  Aj_bg)
    J_all[i, j] <- .jacc(Ai_all, Aj_all)

    # For the "rich" view, only store one triangle (i < j)
    if (j <= i)
      next

    # Intersection counts
    inter_ext <- length(intersect(Ai_ext, Aj_ext))
    inter_bg  <- length(intersect(Ai_bg,  Aj_bg))
    inter_all <- length(intersect(Ai_all, Aj_all))

    # Directional overlap proportions for all-SNP buckets:
    # prop_all_i_in_j = fraction of A's bucket that B also catches
    # prop_all_j_in_i = fraction of B's bucket that A also catches
    prop_all_i_in_j <-
      if (nA_all > 0L)
        inter_all / nA_all
    else
      NA_real_
    prop_all_j_in_i <-
      if (nB_all > 0L)
        inter_all / nB_all
    else
      NA_real_

    # Unique counts in all-SNP buckets
    unique_i_all <- nA_all - inter_all
    unique_j_all <- nB_all - inter_all

    overlap_list[[oi]] <- data.table::data.table(
      i_index = i,
      j_index = j,
      rule_i  = rule_ids[i],
      rule_j  = rule_ids[j],

      # bucket sizes (all SNPs)
      n_all_i = nA_all,
      n_all_j = nB_all,
      n_all_intersect = inter_all,
      n_all_unique_i  = unique_i_all,
      n_all_unique_j  = unique_j_all,

      # extremes & background intersection sizes
      n_ext_i = nA_ext,
      n_ext_j = nB_ext,
      n_bg_i  = nA_bg,
      n_bg_j  = nB_bg,

      n_bg_intersect  = inter_bg,
      n_ext_intersect = inter_ext,

      # directional overlap for all-SNP buckets
      prop_all_i_in_j = prop_all_i_in_j,
      prop_all_j_in_i = prop_all_j_in_i,

      # Jaccard (duplicated from matrices, but handy here)
      jacc_ext = J_ext[i, j],
      jacc_bg  = J_bg [i, j],
      jacc_all = J_all[i, j]
    )
    oi <- oi + 1L
  }
}

# Optional progress over overlap matrix rows
if (!is.null(progress_every) && progress_every > 0L &&
    (i %% progress_every == 0L)) {
  message(sprintf("[%s] computed overlaps for %d / %d focal rules",
                  FUN, i, K))
}
}

overlap_tbl <-
  data.table::rbindlist(overlap_list, use.names = TRUE, fill = TRUE)

list(
  summary     = summary_tbl[],
  jaccard_ext = J_ext,
  jaccard_bg  = J_bg,
  jaccard_all = J_all,
  overlap     = overlap_tbl[],
  rule_ids    = rule_ids
)
}
