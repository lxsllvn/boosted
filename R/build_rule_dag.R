#' Build a dominance DAG over candidate rules for deduplication
#'
#' Constructs a directed acyclic graph (DAG) over the candidate rules returned
#' by \code{\link{analyze_rule_overlap}}, with two sequential operations:
#' (1) near-equivalence collapse, which identifies groups of rules with nearly
#' identical extreme SNP coverage (mutual \code{prop_ext} above
#' \code{equiv_threshold}) and reduces each group to a single representative;
#' and (2) directed dominance edges, where an edge \eqn{i \to j}{} means rule
#' \eqn{i}{} covers most of rule \eqn{j}{}'s extreme SNPs, making \eqn{j}{}
#' potentially redundant given \eqn{i}{}. The \code{mode} argument controls
#' how strictly dominance is defined, with \code{"milp"} requiring an
#' additional background non-inferiority condition appropriate for exact
#' optimisation.
#'
#' The antichain of the resulting DAG — nodes with no incoming dominance edges
#' — is the deduplicated candidate set passed to
#' \code{\link{greedy_select_rules}} or \code{\link{milp_select_rules}}.
#'
#' @param overlap_result Named list returned by
#'   \code{\link{analyze_rule_overlap}}. Must contain \code{overlap},
#'   \code{rule_ids}, \code{rule_summary}, and \code{sets}.
#' @param mode Character string, one of \code{"greedy"} or \code{"milp"}.
#'   Controls dominance edge criteria. In \code{"greedy"} mode a moderate
#'   threshold is used and preferability (fewer clauses or higher \eqn{\omega})
#'   determines edge direction. In \code{"milp"} mode a strict threshold is
#'   used and an edge \eqn{i \to j}{} additionally requires that rule \eqn{i}{}
#'   has no worse background count than rule \eqn{j}{}, ensuring that MILP
#'   cannot improve on the remaining candidates by reinstating a pruned rule.
#' @param equiv_threshold Numeric scalar in \code{(0, 1]}. Mutual
#'   \code{prop_ext} threshold above which two rules are treated as
#'   near-equivalent and collapsed to a single representative. Default
#'   \code{0.95}.
#' @param dominance_threshold Numeric scalar in \code{(0, 1]}. One-sided
#'   \code{prop_ext_j_in_i} threshold above which rule \eqn{i}{} is considered
#'   to dominate rule \eqn{j}{}. \code{NULL} (default) uses \code{0.90} for
#'   \code{"greedy"} and \code{0.99} for \code{"milp"}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{dag}}{An \code{igraph} directed graph over the post-collapse
#'     representative nodes. Vertex attributes: \code{rule_str},
#'     \code{n_clauses}, \code{omega}, \code{n_ext}, \code{n_bg},
#'     \code{is_antichain}.}
#'   \item{\code{equiv_map}}{A \code{data.table} with one row per original
#'     rule, recording \code{rule_str}, \code{rep_rule_str} (the elected
#'     representative for its equivalence class), and \code{component}
#'     (equivalence class ID).}
#'   \item{\code{reps}}{Integer vector of 1-based positions in
#'     \code{rule_ids} corresponding to the elected representatives.}
#'   \item{\code{antichain_reps}}{Integer vector of 1-based positions in
#'     \code{rule_ids} for antichain members (representatives with no
#'     incoming dominance edges).}
#'   \item{\code{antichain_rule_ids}}{Character vector of rule strings for
#'     the antichain members.}
#'   \item{\code{antichain_summary}}{A \code{data.table} with one row per
#'     antichain member, containing the same columns as
#'     \code{overlap_result$rule_summary}.}
#'   \item{\code{antichain_sets}}{Named list with \code{A_sets},
#'     \code{B_sets} (compact extreme and background SNP ID vectors,
#'     one entry per antichain member in \code{antichain_reps} order),
#'     and scalars \code{A_n}, \code{B_n}. Ready for direct use by
#'     \code{\link{greedy_select_rules}} and \code{\link{milp_select_rules}}.}
#'   \item{\code{n_original}}{Integer. Number of candidate rules before any
#'     deduplication.}
#'   \item{\code{n_after_equiv}}{Integer. Number of representatives after
#'     near-equivalence collapse.}
#'   \item{\code{n_antichain}}{Integer. Number of antichain members.}
#'   \item{\code{mode}, \code{equiv_threshold}, \code{dominance_threshold}}{
#'     The parameter values used, for reproducibility.}
#' }
#' @export
build_rule_dag <- function(overlap_result,
                           mode = c("greedy", "milp"),
                           equiv_threshold     = 0.95,
                           dominance_threshold = NULL) {
  FUN  <- "build_rule_dag"
  mode <- match.arg(mode)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(sprintf(
      "[%s] igraph is required. Install with install.packages('igraph').", FUN
    ))
  }

  # Validate inputs
  required <- c("overlap", "rule_ids", "rule_summary", "sets")
  missing  <- setdiff(required, names(overlap_result))
  if (length(missing)) {
    stop(sprintf(
      "[%s] overlap_result missing elements: %s. Re-run analyze_rule_overlap().",
      FUN, paste(missing, collapse = ", ")
    ))
  }

  overlap      <- overlap_result$overlap
  rule_ids     <- overlap_result$rule_ids
  summary      <- overlap_result$rule_summary
  A_sets_full  <- overlap_result$sets$A_sets
  B_sets_full  <- overlap_result$sets$B_sets
  A_n          <- overlap_result$sets$A_n
  B_n          <- overlap_result$sets$B_n
  K            <- length(rule_ids)

  # Check for directional extreme columns added by the updated
  # analyze_rule_overlap(); error clearly if running against old output.
  if (nrow(overlap) > 0L) {
    req_cols <- c("prop_ext_i_in_j", "prop_ext_j_in_i")
    missing_cols <- setdiff(req_cols, names(overlap))
    if (length(missing_cols)) {
      stop(sprintf(
        "[%s] overlap table missing columns: %s. Re-run analyze_rule_overlap().",
        FUN, paste(missing_cols, collapse = ", ")
      ))
    }
  }

  # Set default dominance threshold per mode
  if (is.null(dominance_threshold)) {
    dominance_threshold <- if (mode == "greedy") 0.90 else 0.99
  }

  # -------------------------------------------------------------------
  # Step 1: Near-equivalence collapse (undirected)
  # Two rules are near-equivalent when each covers most of the other's
  # extreme SNPs. The relation is symmetric so we use an undirected graph
  # and find connected components.
  # -------------------------------------------------------------------
  if (nrow(overlap) > 0L) {
    eq_mask <- !is.na(overlap$prop_ext_i_in_j) &
               !is.na(overlap$prop_ext_j_in_i) &
               overlap$prop_ext_i_in_j >= equiv_threshold &
               overlap$prop_ext_j_in_i >= equiv_threshold

    eq_rows <- overlap[eq_mask]
  } else {
    eq_rows <- overlap[integer(0)]
  }

  # Build undirected equivalence graph; all K rules are vertices so that
  # isolated nodes (no equivalence edges) form singleton components.
  all_verts <- data.frame(name = seq_len(K), stringsAsFactors = FALSE)

  if (nrow(eq_rows) > 0L) {
    eq_edges <- data.frame(
      from = eq_rows$i_index,
      to   = eq_rows$j_index,
      stringsAsFactors = FALSE
    )
    g_equiv <- igraph::graph_from_data_frame(
      eq_edges,
      directed = FALSE,
      vertices = all_verts
    )
  } else {
    g_equiv <- igraph::graph_from_data_frame(
      data.frame(from = integer(0), to = integer(0), stringsAsFactors = FALSE),
      directed = FALSE,
      vertices = all_verts
    )
  }

  comps      <- igraph::components(g_equiv)
  membership <- comps$membership   # length K; 1-based component ID per rule

  # Elect one representative per component: prefer fewest clauses, then
  # highest omega, then lexicographically first rule string for determinism.
  rep_of   <- integer(K)           # rep_of[k] = original 1..K index of k's rep
  comp_ids <- sort(unique(membership))

  for (cid in comp_ids) {
    members <- which(membership == cid)
    if (length(members) == 1L) {
      rep_of[members] <- members
    } else {
      nc  <- summary$n_clauses[members]
      om  <- summary$omega[members]
      rs  <- summary$rule_str[members]
      best_local <- order(nc, -om, rs)[1L]
      rep_of[members] <- members[best_local]
    }
  }

  reps   <- sort(unique(rep_of))   # unique representative original indices
  n_reps <- length(reps)

  # Map original index -> position in reps (0 if not a representative)
  rep_pos          <- integer(K)
  rep_pos[reps]    <- seq_len(n_reps)

  # -------------------------------------------------------------------
  # Step 2: Directed dominance edges among representatives
  # -------------------------------------------------------------------
  is_rep <- logical(K)
  is_rep[reps] <- TRUE

  # Filter overlap to pairs where both rules are representatives from
  # different equivalence classes.
  if (nrow(overlap) > 0L) {
    dom_rows <- overlap[is_rep[overlap$i_index] & is_rep[overlap$j_index]]
  } else {
    dom_rows <- overlap[integer(0)]
  }

  edge_from_pos <- integer(0)
  edge_to_pos   <- integer(0)

  if (nrow(dom_rows) > 0L) {
    pi_j <- dom_rows$prop_ext_i_in_j   # fraction of i's extremes in j
    pj_i <- dom_rows$prop_ext_j_in_i   # fraction of j's extremes in i

    nc_i <- summary$n_clauses[dom_rows$i_index]
    nc_j <- summary$n_clauses[dom_rows$j_index]
    om_i <- summary$omega[dom_rows$i_index]
    om_j <- summary$omega[dom_rows$j_index]
    bg_i <- summary$n_bg[dom_rows$i_index]
    bg_j <- summary$n_bg[dom_rows$j_index]

    if (mode == "greedy") {
      # i is preferable to j: fewer clauses, or equal clauses + higher omega
      i_pref_j <- (nc_i < nc_j) | (nc_i == nc_j & om_i > om_j)
      j_pref_i <- (nc_j < nc_i) | (nc_j == nc_i & om_j > om_i)

      # Forward edge i->j: j's extremes nearly inside i, and i preferred
      fwd <- !is.na(pj_i) & pj_i >= dominance_threshold & i_pref_j
      # Reverse edge j->i: i's extremes nearly inside j, and j preferred
      rev <- !is.na(pi_j) & pi_j >= dominance_threshold & j_pref_i

    } else {
      # MILP mode: strict dominance requires both extreme containment AND
      # i does not introduce more background than j.
      fwd <- !is.na(pj_i) & pj_i >= dominance_threshold & bg_i <= bg_j
      rev <- !is.na(pi_j) & pi_j >= dominance_threshold & bg_j <= bg_i
    }

    # Translate overlap-table indices to DAG vertex positions (1..n_reps)
    from_fwd <- rep_pos[dom_rows$i_index[fwd]]
    to_fwd   <- rep_pos[dom_rows$j_index[fwd]]
    from_rev <- rep_pos[dom_rows$j_index[rev]]
    to_rev   <- rep_pos[dom_rows$i_index[rev]]

    edge_from_pos <- c(from_fwd, from_rev)
    edge_to_pos   <- c(to_fwd,   to_rev)

    # Remove self-loops (can arise when the same rule qualifies in both
    # directions at the boundary of equal preferability).
    ok <- edge_from_pos != edge_to_pos
    edge_from_pos <- edge_from_pos[ok]
    edge_to_pos   <- edge_to_pos[ok]
  }

  # Build DAG
  dag_verts <- data.frame(name = seq_len(n_reps), stringsAsFactors = FALSE)

  if (length(edge_from_pos) > 0L) {
    dag_edges <- data.frame(
      from = edge_from_pos,
      to   = edge_to_pos,
      stringsAsFactors = FALSE
    )
    dag <- igraph::graph_from_data_frame(
      dag_edges,
      directed = TRUE,
      vertices = dag_verts
    )

    # Guard against cycles. These should not arise given the preferability
    # conditions above, but binning removal can create near-ties in omega
    # or n_clauses that generate mutual edges.
    if (!igraph::is_dag(dag)) {
      warning(sprintf(
        paste("[%s] dominance graph contains cycles;",
              "breaking by condensing strongly connected components."),
        FUN
      ))
      scc      <- igraph::components(dag, mode = "strong")
      scc_mem  <- scc$membership   # length n_reps; SCC ID per vertex

      # For each SCC with > 1 member, elect one representative by the
      # same preferability criterion as the equivalence collapse.
      new_reps_orig <- vapply(sort(unique(scc_mem)), function(sid) {
        vpos  <- which(scc_mem == sid)             # positions in 1..n_reps
        korig <- reps[vpos]                        # original 1..K indices
        nc    <- summary$n_clauses[korig]
        om    <- summary$omega[korig]
        rs    <- summary$rule_str[korig]
        korig[order(nc, -om, rs)[1L]]
      }, integer(1L))

      reps   <- new_reps_orig
      n_reps <- length(reps)
      rep_pos <- integer(K)
      rep_pos[reps] <- seq_len(n_reps)

      # Remap edges, drop intra-SCC edges
      new_from <- rep_pos[reps[edge_from_pos]]
      new_to   <- rep_pos[reps[edge_to_pos]]
      ok2 <- !is.na(new_from) & !is.na(new_to) &
             new_from > 0L & new_to > 0L & new_from != new_to
      edge_from_pos <- new_from[ok2]
      edge_to_pos   <- new_to[ok2]

      dag_verts <- data.frame(name = seq_len(n_reps), stringsAsFactors = FALSE)
      if (length(edge_from_pos)) {
        dag <- igraph::graph_from_data_frame(
          data.frame(from = edge_from_pos, to = edge_to_pos,
                     stringsAsFactors = FALSE),
          directed = TRUE, vertices = dag_verts
        )
      } else {
        dag <- igraph::graph_from_data_frame(
          data.frame(from = integer(0), to = integer(0),
                     stringsAsFactors = FALSE),
          directed = TRUE, vertices = dag_verts
        )
      }
    }
  } else {
    dag <- igraph::graph_from_data_frame(
      data.frame(from = integer(0), to = integer(0), stringsAsFactors = FALSE),
      directed = TRUE,
      vertices = dag_verts
    )
  }

  # Antichain: representatives with no incoming dominance edges.
  # These are the rules that are not subsumed by any other surviving rule
  # and form the candidate set for downstream selection.
  in_deg         <- igraph::degree(dag, mode = "in")
  antichain_pos  <- which(in_deg == 0L)          # positions in 1..n_reps
  antichain_reps <- reps[antichain_pos]           # original 1..K indices

  # Attach per-vertex attributes for downstream use and visualization
  igraph::V(dag)$rule_str     <- rule_ids[reps]
  igraph::V(dag)$n_clauses    <- summary$n_clauses[reps]
  igraph::V(dag)$omega        <- summary$omega[reps]
  igraph::V(dag)$n_ext        <- summary$n_ext[reps]
  igraph::V(dag)$n_bg         <- summary$n_bg[reps]
  igraph::V(dag)$is_antichain <- in_deg == 0L

  # Equivalence map: original rule -> representative
  equiv_map <- data.table::data.table(
    rule_str     = rule_ids,
    rep_rule_str = rule_ids[rep_of],
    component    = as.integer(membership)
  )

  # Antichain summary (same columns as rule_summary, restricted to antichain)
  antichain_summary <- summary[antichain_reps, ]

  # Pre-indexed A/B sets for antichain rules, ready for selection functions
  antichain_sets <- list(
    A_sets = A_sets_full[antichain_reps],
    B_sets = B_sets_full[antichain_reps],
    A_n    = A_n,
    B_n    = B_n
  )

  message(sprintf(
    "[%s] K=%d rules -> %d after equiv collapse -> %d antichain members (mode=%s)",
    FUN, K, n_reps, length(antichain_reps), mode
  ))

  list(
    dag                  = dag,
    equiv_map            = equiv_map[],
    reps                 = reps,
    antichain_reps       = antichain_reps,
    antichain_rule_ids   = rule_ids[antichain_reps],
    antichain_summary    = antichain_summary[],
    antichain_sets       = antichain_sets,
    n_original           = K,
    n_after_equiv        = n_reps,
    n_antichain          = length(antichain_reps),
    mode                 = mode,
    equiv_threshold      = equiv_threshold,
    dominance_threshold  = dominance_threshold
  )
}
