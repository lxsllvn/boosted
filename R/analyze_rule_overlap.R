#' Pairwise SNP overlap analysis between candidate rules
#'
#' For each pair of candidate rules, computes the intersection and union of
#' their SNP buckets (the sets of SNPs covered by matching
#' \code{(Tree, leaf_id)} pairs across the ensemble) and derives Jaccard
#' similarity coefficients and directional overlap proportions, separately for
#' all SNPs, extreme SNPs, and background SNPs. A sparse incidence matrix
#' approach replaces \eqn{O(K^2)}{} pairwise set-intersection loops, making
#' the computation tractable for large rule sets. The output supports
#' downstream deduplication via \code{\link{build_rule_dag}} and rule
#' selection via \code{\link{greedy_select_rules}} or
#' \code{\link{milp_select_rules}}.
#'
#' @inheritParams .boosted_params
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{overlap}}{A \code{data.table} with one row per ordered pair
#'     \code{(i, j)} with \code{i < j} and \code{n_all_intersect > 0}.
#'     In addition to the columns produced previously, now includes
#'     \code{prop_ext_i_in_j} (fraction of rule \code{i}'s extreme SNPs
#'     covered by rule \code{j}) and \code{prop_ext_j_in_i} (fraction of
#'     rule \code{j}'s extreme SNPs covered by rule \code{i}). Entirely-\code{NA}
#'     columns are dropped.}
#'   \item{\code{rule_ids}}{Character vector of length \code{K}: the sorted
#'     rule strings corresponding to rows/columns of the Jaccard matrices and
#'     to rows of \code{rule_summary}.}
#'   \item{\code{rule_summary}}{A \code{data.table} with one row per entry
#'     in \code{rule_ids}, containing: \code{rule_str}, \code{n_clauses},
#'     \code{n_ext}, \code{n_bg}, \code{n_all}, and \code{omega} (the
#'     rule-level log-likelihood ratio with Jeffreys smoothing at
#'     \code{alpha}).}
#'   \item{\code{jaccard_ext}, \code{jaccard_bg}, \code{jaccard_all}}{
#'     Numeric \eqn{K \times K}{} matrices of Jaccard coefficients computed
#'     over extreme, background, and all SNPs respectively. Diagonal entries
#'     are 1.}
#'   \item{\code{sets}}{A list with elements \code{A_sets} and \code{B_sets}
#'     (per-rule compact extreme and background SNP ID vectors, indexed
#'     1-based into \code{rule_ids} order) and scalars \code{A_n} and
#'     \code{B_n} (total extreme and background set sizes).}
#'   \item{\code{meta}}{A list recording \code{which}, \code{fold_n_all},
#'     \code{fold_n_extr}, and \code{fold_n_bg}.}
#' }
#' @export
analyze_rule_overlap <- function(boosted,
                                 candidate_rules = NULL,
                                 which = c("train", "test"),
                                 fold_indices = NULL,
                                 alpha = 0.5,
                                 progress_every = NULL) {
  # Signature & basic checks
  FUN <- "analyze_rule_overlap"
  message(sprintf("[%s] start: %s", FUN, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  if (!inherits(boosted, "boosted_harvest")) {
    stop(sprintf("[%s] boosted is not ready; run prepare_harvest() first.",
                 FUN))
  }

  # Determine if we're working with test or training data
  which <- match.arg(which)

  # Pull test/training indices from `boosted` as requested
  extr_idx  <- boosted[[sprintf("extr_idx_%s",  which)]]
  bg_idx    <- boosted[[sprintf("bg_idx_%s",    which)]]
  n_tot     <- boosted[[sprintf("n_yvar_%s",    which)]]
  N_extr    <- boosted[[sprintf("N_extr_%s",    which)]]
  N_bg      <- boosted[[sprintf("N_bg_%s",      which)]]

  # Pull (Tree, leaf_id) -> SNP map
  snps_all_by_leaf <- boosted[[sprintf("snps_all_by_leaf_%s",   which)]]
  # Pull rule -> (Tree, leaf_id) map
  leaf_rule_cache <- boosted$leaf_rule_cache

  # If fold_indices != null, ensure that they are a valid vector of integers
  if (is.null(fold_indices)) {
    all_idx <- seq_len(n_tot)
  } else {
    all_idx <- .check_idx(fold_indices, n_tot, FUN, "fold_indices")
    # Restrict background/extreme labels to fold_indices
    extr_idx  <- intersect(extr_idx, all_idx)
    bg_idx    <- intersect(bg_idx, all_idx)
    N_extr    <- length(extr_idx)
    N_bg      <- length(bg_idx)
    n_tot     <- length(all_idx)
  }

  # Membership vectors for fast counting within buckets
  in_fold    <- rep(FALSE, n_tot)
  is_extreme <- rep(FALSE, n_tot)
  is_bg      <- rep(FALSE, n_tot)

  in_fold[all_idx]     <- TRUE
  is_extreme[extr_idx] <- TRUE
  is_bg[bg_idx]        <- TRUE

  # Check if the candidate rule strings are found in `boosted$leaf_rule_cache`
  # and return the rows with (Tree, leaf_id) under each rule.
  rules <- .validate_candidate_rules(
    candidate_rules = candidate_rules,
    leaf_rule_cache = leaf_rule_cache,
    caller          = FUN
  )

  uniq_rules  <- rules$candidate_rules
  pairs_all   <- rules$pairs_all

  # Rule universe for overlap
  rule_ids <- sort(unique(pairs_all$rule_str))
  K <- length(rule_ids)

  if (K < 2L) {
    warning(sprintf("[%s] fewer than two rules selected; nothing to compare.", FUN))
    return(invisible(NULL))
  }

  # Pre-split pairs_all by rule_str once to avoid repeated forderv calls
  # inside the bucket-building loop
  pairs_by_rule <- split(pairs_all, by = "rule_str")

  # n_clauses per rule, in rule_ids order, for inclusion in overlap_tbl
  n_clauses_by_rule <- vapply(
    rule_ids,
    function(rs) pairs_by_rule[[rs]]$n_clauses[1L],
    integer(1L)
  )

  # Precompute buckets per rule (fold-restricted) and counts
  bucket_ext <- vector("list", K)
  bucket_bg  <- vector("list", K)
  bucket_all <- vector("list", K)
  n_all <- integer(K)
  n_ext <- integer(K)
  n_bg  <- integer(K)

  # Map global SNP index -> element id within the A/B universes (1..N_extr / 1..N_bg)
  map_ext <- integer(n_tot); map_ext[extr_idx] <- seq_along(extr_idx)
  map_bg  <- integer(n_tot); map_bg[bg_idx]    <- seq_along(bg_idx)

  A_sets <- vector("list", K)
  B_sets <- vector("list", K)

  # Build SNP buckets per rule string.
  for (k in seq_len(K)) {
    rs <- rule_ids[k]

    # All (Tree, leaf_id) pairs under this rule
    pairs <- pairs_by_rule[[rs]]

    # Look up SNP indices covered by these (Tree, leaf_id) pairs.
    # .snp_lookup() handles deduplication (no double-counting SNPs
    # that land in the same rule via multiple trees).
    b <- .snp_lookup(
      pairs            = pairs,
      snps_all_by_leaf = snps_all_by_leaf
    )$bucket_all

    # Restrict bucket to fold_indices
    if (length(b)) {
      b <- b[in_fold[b]]
    }
    bucket_all[[k]] <- b

    # Split buckets into extreme/background subsets (fold-restricted)
    # Needed for tcrossprod-based intersections
    if (length(b)) {
      a <- b[is_extreme[b]]
      g <- b[is_bg[b]]

      bucket_ext[[k]] <- a
      bucket_bg[[k]]  <- g

      # Convert global SNP indices to consecutive element IDs for downstream
      # selection functions
      A_sets[[k]] <- map_ext[a]
      B_sets[[k]] <- map_bg[g]
    } else {
      bucket_ext[[k]] <- integer(0)
      bucket_bg[[k]]  <- integer(0)
      A_sets[[k]] <- integer(0)
      B_sets[[k]] <- integer(0)
    }

    # Optional progress over rule buckets
    if (!is.null(progress_every) && progress_every > 0L &&
        (k %% progress_every == 0L)) {
      message(sprintf("[%s] built SNP buckets for %d / %d rules",
                      FUN, k, K))
    }
  }

  # Compute pairwise intersections in compiled code via sparse incidence
  # matrices and tcrossprod(). This replaces the O(K^2) pairwise
  # set-intersection loops.

  # Build a compact SNP universe over the selected rules (all buckets)
  snps_universe <- sort.int(unique(unlist(bucket_all, use.names = FALSE)))
  N_used <- length(snps_universe)

  # If all selected rules are empty, return nothing.
  if (N_used == 0L) {
    warning(sprintf("[%s] all selected rules are empty.", FUN))
    return(invisible(NULL))
  }

  # SNP indices in buckets_* are row positions in the train/test leaf
  # matrices, so they are not contiguous. We remap them to a compact 1..N
  # index first.
  buckets_all_m <- lapply(bucket_all, function(b) match(b, snps_universe))

  # Build sparse incidence matrix: rows = rules, cols = SNPs (remapped 1..N_used)
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

  # Build incidence matrix once (rules x SNPs)
  M_all <- .build_incidence(buckets_all_m)

  # Column subsets for extreme/background SNPs within the compact universe
  ext_cols <- which(is_extreme[snps_universe])
  bg_cols  <- which(is_bg[snps_universe])

  # Intersection count matrices (K x K)
  I_all <- as.matrix(Matrix::tcrossprod(M_all))
  I_ext <- as.matrix(Matrix::tcrossprod(M_all[, ext_cols, drop = FALSE]))
  I_bg  <- as.matrix(Matrix::tcrossprod(M_all[, bg_cols,  drop = FALSE]))

  # Sizes per rule
  n_all <- as.integer(diag(I_all))
  n_ext <- as.integer(diag(I_ext))
  n_bg  <- as.integer(diag(I_bg))

  # Per-rule summary table: one row per entry in rule_ids.
  # omega is the rule-level LLR with Jeffreys smoothing at alpha.
  rule_summary <- data.table::data.table(
    rule_str  = rule_ids,
    n_clauses = n_clauses_by_rule,
    n_ext     = n_ext,
    n_bg      = n_bg,
    n_all     = n_all,
    omega     = vapply(seq_len(K), function(k) {
      .rule_llr(
        n_extreme    = n_ext[k],
        n_bg         = n_bg[k],
        N_extr_total = N_extr,
        N_bg_total   = N_bg,
        alpha        = alpha
      )
    }, numeric(1L))
  )

  # Jaccard matrices over extreme / background / all SNPs
  # Use the same convention as before: if union is 0, Jaccard = 0.
  .jaccard_from_intersections <- function(I, nA) {
    U <- outer(nA, nA, "+") - I
    J <- matrix(0, nrow = K, ncol = K)
    ok <- (U > 0)
    J[ok] <- I[ok] / U[ok]
    diag(J) <- 1
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
  keep <- (inter_all > 0)
  if (!any(keep)) {
    overlap_tbl <- data.table::data.table()
  } else {
    i_index <- i_index[keep]
    j_index <- j_index[keep]

    inter_all <- inter_all[keep]
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

    prop_all_i_in_j <- ifelse(n_all_i > 0L, n_all_intersect / n_all_i, NA_real_)
    prop_all_j_in_i <- ifelse(n_all_j > 0L, n_all_intersect / n_all_j, NA_real_)

    # Directional containment on extreme SNPs only.
    # prop_ext_i_in_j: fraction of i's extreme SNPs that are also in j's bucket.
    # prop_ext_j_in_i: fraction of j's extreme SNPs that are also in i's bucket.
    # These are the primary edge weights for the dominance DAG in
    # build_rule_dag(): a high prop_ext_j_in_i means i nearly subsumes j's
    # extreme coverage, making j a candidate for removal.
    prop_ext_i_in_j <- ifelse(n_ext_i > 0L, n_ext_intersect / n_ext_i, NA_real_)
    prop_ext_j_in_i <- ifelse(n_ext_j > 0L, n_ext_intersect / n_ext_j, NA_real_)

    # Jaccard vectors (duplicated from matrices, but handy here)
    jacc_all <- J_all[cbind(i_index, j_index)]
    jacc_ext <- J_ext[cbind(i_index, j_index)]
    jacc_bg  <- J_bg [cbind(i_index, j_index)]

    overlap_all <- n_all_intersect / pmax.int(1L, pmin.int(n_all_i, n_all_j))
    overlap_ext <- ifelse(
      pmin.int(n_ext_i, n_ext_j) > 0L,
      n_ext_intersect / pmin.int(n_ext_i, n_ext_j),
      NA_real_
    )
    overlap_bg  <- ifelse(
      pmin.int(n_bg_i,  n_bg_j)  > 0L,
      n_bg_intersect / pmin.int(n_bg_i,  n_bg_j),
      NA_real_
    )

    overlap_tbl <- data.table::data.table(
      i_index   = i_index,
      j_index   = j_index,
      n_clauses_i = as.integer(n_clauses_by_rule[i_index]),
      n_clauses_j = as.integer(n_clauses_by_rule[j_index]),

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

      # directional containment on extreme SNPs (primary DAG edge weights)
      prop_ext_i_in_j = prop_ext_i_in_j,
      prop_ext_j_in_i = prop_ext_j_in_i,

      # Jaccard (duplicated from matrices, but handy here)
      jacc_ext = jacc_ext,
      jacc_bg  = jacc_bg,
      jacc_all = jacc_all,

      overlap_all = overlap_all,
      overlap_ext = overlap_ext,
      overlap_bg  = overlap_bg
    )
  }

  overlap_tbl <- .drop_all_na_cols(overlap_tbl)

  list(
    overlap      = overlap_tbl[],
    rule_ids     = rule_ids,
    rule_summary = rule_summary[],
    jaccard_ext  = J_ext,
    jaccard_bg   = J_bg,
    jaccard_all  = J_all,
    sets         = list(
      A_sets = A_sets,
      B_sets = B_sets,
      A_n    = length(extr_idx),
      B_n    = length(bg_idx)),
    meta         = list(
      which       = which,
      fold_n_all  = n_tot,
      fold_n_extr = N_extr,
      fold_n_bg   = N_bg
    )
  )
}
