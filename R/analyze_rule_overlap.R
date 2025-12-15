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

  if (!data.table::haskey(pairs_all) ||
      !identical(data.table::key(pairs_all), "rule_str")) {
    data.table::setkey(pairs_all, rule_str)
  }

  for (k in seq_len(K)) {
    rs <- rule_ids[k]

    # All (Tree, leaf_id) pairs where this rule held in the ensemble
    pr <- unique(pairs_all[J(rs), .(Tree, leaf_id)])
    if (!nrow(pr)) {
      buckets_ext[[k]] <- integer(0)
      buckets_bg [[k]] <- integer(0)
      buckets_all[[k]] <- integer(0)
      next
    } else {
      # Look up SNP indices covered by these (Tree, leaf_id) pairs.
      # .snp_lookup() handles deduplication (no double-counting SNPs
      # that land in the same rule via multiple trees).
      buckets <- .snp_lookup(
        pairs = pr,
        snps_ext_by_leaf = snps_ext_by_leaf,
        snps_bg_by_leaf  = snps_bg_by_leaf,
        snps_all_by_leaf = snps_all_by_leaf
      )
      buckets_ext[[k]] <- sort.int(unique(buckets$bucket_ext))
      buckets_bg [[k]] <- sort.int(unique(buckets$bucket_bg))
      buckets_all[[k]] <- sort.int(unique(buckets$bucket_all))
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

  # Compute pairwise intersections in compiled code via sparse incidence
  # matrices and tcrossprod(). This replaces the O(K^2) pairwise
  # set-intersection loops.
  #
  # Note: SNP indices in buckets_* are row positions in the train/test leaf
  # matrices, so they are not contiguous. We remap them to a compact 1..N_used
  # index first.

  # Build a compact SNP universe over the selected rules (all buckets)
  snps_universe <-
    sort.int(unique(unlist(buckets_all, use.names = FALSE)))
  N_used <- length(snps_universe)

  # Guard: if all selected rules are empty, return empty structures
  if (N_used == 0L) {
    labs <- paste0("r", seq_len(K))

    J_ext <-
      matrix(0,
             nrow = K,
             ncol = K,
             dimnames = list(labs, labs))
    J_bg  <-
      matrix(0,
             nrow = K,
             ncol = K,
             dimnames = list(labs, labs))
    J_all <-
      matrix(0,
             nrow = K,
             ncol = K,
             dimnames = list(labs, labs))

    overlap_tbl <- data.table::data.table()

    return(
      list(
        summary     = summary_tbl[],
        jaccard_ext = J_ext,
        jaccard_bg  = J_bg,
        jaccard_all = J_all,
        overlap     = overlap_tbl[],
        rule_ids    = rule_ids
      )
    )
  }

  # Remap each bucket to 1..N_used (based on snps_universe)
  buckets_all_m <-
    lapply(buckets_all, function(b)
      match(b, snps_universe))
  buckets_ext_m <-
    lapply(buckets_ext, function(b)
      match(b, snps_universe))
  buckets_bg_m  <-
    lapply(buckets_bg,  function(b)
      match(b, snps_universe))

  # Build sparse incidence matrices: rows = rules, cols = SNPs (remapped 1..N_used)
  .build_incidence <- function(buckets_m) {
    lens <- lengths(buckets_m)
    if (!any(lens)) {
      return(Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = 1L,
        dims = c(K, N_used)
      ))
    }
    i <- rep.int(seq_len(K), lens)
    j <- unlist(buckets_m, use.names = FALSE)
    Matrix::sparseMatrix(
      i = i,
      j = j,
      x = 1L,
      dims = c(K, N_used)
    )
  }

  M_all <- .build_incidence(buckets_all_m)
  M_ext <- .build_incidence(buckets_ext_m)
  M_bg  <- .build_incidence(buckets_bg_m)

  # Intersection count matrices (K x K)
  C_all <- Matrix::tcrossprod(M_all)
  C_ext <- Matrix::tcrossprod(M_ext)
  C_bg  <- Matrix::tcrossprod(M_bg)

  # Convert to dense matrices for simple downstream arithmetic (K is small; e.g., 2000 → ~32 MB per matrix)
  I_all <- as.matrix(C_all)
  I_ext <- as.matrix(C_ext)
  I_bg  <- as.matrix(C_bg)

  # Sizes per rule
  n_all <- diag(I_all)
  n_ext <- diag(I_ext)
  n_bg  <- diag(I_bg)

  # Jaccard matrices over extreme / background / all SNPs
  # Use the same convention as before: if union is 0, Jaccard = 0.
  .jaccard_from_intersections <- function(I, nA) {
    U <- outer(nA, nA, "+") - I
    J <- matrix(0, nrow = K, ncol = K)
    ok <- (U > 0)
    J[ok] <- I[ok] / U[ok]
    J
  }

  labs <- paste0("r", seq_len(K))
  J_all <- .jaccard_from_intersections(I_all, n_all)
  J_ext <- .jaccard_from_intersections(I_ext, n_ext)
  J_bg  <- .jaccard_from_intersections(I_bg,  n_bg)

  dimnames(J_all) <- list(labs, labs)
  dimnames(J_ext) <- list(labs, labs)
  dimnames(J_bg)  <- list(labs, labs)

  # Also build a rich long-form overlap table:
  # for each (i,j) we store intersection sizes and directional proportions.
  n_pairs <- K * (K - 1L) / 2L

  i_index <- integer(n_pairs)
  j_index <- integer(n_pairs)

  # Fill (i,j) indices for i < j without allocating large row()/col() matrices
  pos <- 1L
  for (i in seq_len(K - 1L)) {
    nj <- K - i
    idx <- pos:(pos + nj - 1L)
    i_index[idx] <- i
    j_index[idx] <- (i + 1L):K
    pos <- pos + nj
  }

  # Extract intersections for i < j
  inter_all <- I_all[cbind(i_index, j_index)]
  inter_ext <- I_ext[cbind(i_index, j_index)]
  inter_bg  <- I_bg [cbind(i_index, j_index)]

  n_all_i <- as.integer(n_all[i_index])
  n_all_j <- as.integer(n_all[j_index])

  n_ext_i <- as.integer(n_ext[i_index])
  n_ext_j <- as.integer(n_ext[j_index])

  n_bg_i <- as.integer(n_bg[i_index])
  n_bg_j <- as.integer(n_bg[j_index])

  n_all_intersect <- as.integer(inter_all)
  n_ext_intersect <- as.integer(inter_ext)
  n_bg_intersect  <- as.integer(inter_bg)

  n_all_unique_i  <- n_all_i - n_all_intersect
  n_all_unique_j  <- n_all_j - n_all_intersect

  prop_all_i_in_j <-
    ifelse(n_all_i > 0L, n_all_intersect / n_all_i, NA_real_)
  prop_all_j_in_i <-
    ifelse(n_all_j > 0L, n_all_intersect / n_all_j, NA_real_)

  # Jaccard vectors (duplicated from matrices, but handy here)
  jacc_all <- J_all[cbind(i_index, j_index)]
  jacc_ext <- J_ext[cbind(i_index, j_index)]
  jacc_bg  <- J_bg [cbind(i_index, j_index)]

  overlap_tbl <- data.table::data.table(
    i_index = i_index,
    j_index = j_index,

    # bucket sizes (all SNPs)
    n_all_i = n_all_i,
    n_all_j = n_all_j,
    n_all_intersect = n_all_intersect,
    n_all_unique_i  = n_all_unique_i,
    n_all_unique_j  = n_all_unique_j,

    # extremes & background intersection sizes
    n_ext_i = n_ext_i,
    n_ext_j = n_ext_j,
    n_bg_i  = n_bg_i,
    n_bg_j  = n_bg_j,

    n_bg_intersect  = n_bg_intersect,
    n_ext_intersect = n_ext_intersect,

    # directional overlap for all-SNP buckets
    prop_all_i_in_j = prop_all_i_in_j,
    prop_all_j_in_i = prop_all_j_in_i,

    # Jaccard (duplicated from matrices, but handy here)
    jacc_ext = jacc_ext,
    jacc_bg  = jacc_bg,
    jacc_all = jacc_all
  )

  overlap_tbl[, rule_i := factor(i_index, levels = seq_len(K), labels = rule_ids)]
  overlap_tbl[, rule_j := factor(j_index, levels = seq_len(K), labels = rule_ids)]
  data.table::setcolorder(overlap_tbl,
                          c("i_index", "j_index", "rule_i", "rule_j",
                            setdiff(
                              names(overlap_tbl),
                              c("i_index", "j_index", "rule_i", "rule_j")
                            )))

  list(
    summary     = summary_tbl[],
    jaccard_ext = J_ext,
    jaccard_bg  = J_bg,
    jaccard_all = J_all,
    overlap     = overlap_tbl[],
    rule_ids    = rule_ids
  )
}
