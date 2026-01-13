# function to split genes up

#' Fit per-state GLMs on pseudobulked CCS data
#'
#' Builds pseudobulks for the requested cell groups and fits a GLM per gene,
#' returning coefficient matrices suitable for downstream DEG scoring.
#'
#' @param ccs A CCS object to pseudobulk.
#' @param cell_groups Character vector of cell groups to include.
#' @param genes_to_test Optional integer indices of genes to keep.
#' @param group_nodes_by Optional column name to group states by.
#' @param assembly_group Optional assembly group filter (currently unused).
#' @param gene_ids Optional character vector of gene IDs to keep.
#' @param min_samples_detected Minimum number of samples a gene must be detected in.
#' @param min_cells_per_pseudobulk Minimum cells per pseudobulk to retain.
#' @param nuisance_model_formula_str Nuisance formula string (without response).
#' @param abs_expr_thresh Absolute expression threshold for shrinkage.
#' @param cores Number of cores for model fitting.
#'
#' @return A list-like object from `collect_coefficients_for_shrinkage()` that
#'   includes coefficient and standard error matrices.
compute_glm <- function(ccs,
                        cell_groups,
                        genes_to_test = NULL,
                        group_nodes_by = NULL,
                        assembly_group = NULL,
                        gene_ids = NULL,
                        min_samples_detected = 10,
                        min_cells_per_pseudobulk = 10,
                        nuisance_model_formula_str = "0",
                        abs_expr_thresh = 0.1,
                        cores = 1) {
  
  if (is.null(group_nodes_by)) {
    pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun = "sum")
    state_term <- "cell_group"
  } else {
    pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, state_col = group_nodes_by, cell_agg_fun = "sum")
    state_term <- group_nodes_by
  }

  if (is.null(gene_ids) == FALSE) {
    pb_cds <- pb_cds[intersect(rownames(pb_cds), gene_ids), ]
  }

  pb_cds <- pb_cds[, colData(pb_cds)$cell_group %in% cell_groups]

  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  
  if (is.null(genes_to_test)) {
    expr_over_thresh <- monocle3::normalized_counts(pb_cds, "size_only", pseudocount = 0)
    genes_to_test <- which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  }
  
  pb_cds <- pb_cds[genes_to_test, ]
  
  pseudobulks_to_test <- which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)

  message("fitting regression models")
  pb_cds <- pb_cds[, pseudobulks_to_test]

  colData(pb_cds)$Size_Factor <- colData(pb_cds)$num_cells_in_group

  full_model_str <- paste0("~ 0 + cell_group")
  nuisance_model_formula_str <- stringr::str_replace_all(nuisance_model_formula_str, "~", "")

  if (nuisance_model_formula_str != "0") {
    full_model_str <- paste0(full_model_str, " + ", nuisance_model_formula_str)
  }

  pb_group_models <- monocle3::fit_models(pb_cds,
    model_formula_str = full_model_str,
    cores = cores
  ) %>% dplyr::select(gene_short_name, id, model, model_summary, status)

  message("      collecting coefficients")

  n <- dim(pb_cds)[2]
  pb_coeffs <- collect_coefficients_for_shrinkage(pb_cds, pb_group_models, abs_expr_thresh, term_to_keep = "cell_group")
  rm(pb_group_models)

  return(pb_coeffs)
}

#' Combine gene model coefficient matrices
#'
#' Merges coefficients and unscaled standard errors across multiple model
#' batches into a single tibble with matrix columns.
#'
#' @param gene_model_list List of objects containing `coefficients` and
#'   `stdev.unscaled` matrices.
#'
#' @return A tibble with `coefficients` and `stdev.unscaled` matrix columns.
combine_gene_models <- function(gene_model_list) {
  
  pb_coeffs <- tibble(
    coefficients = do.call(cbind, map(gene_model_list, ~ .x$coefficients)),
    stdev.unscaled = do.call(cbind, map(gene_model_list, ~ .x$stdev.unscaled)),
    # sigma = unlist(map(coefs_list, ~ .x$sigma)),
    # df.residual = unlist(map(coefs_list, ~ .x$df.residual)),
    # trend = unlist(map(coefs_list, ~ .x$trend)),
    # est_dispersion = unlist(map(coefs_list, ~ .x$est_dispersion)),
    # disp_fit = unlist(map(coefs_list, ~ .x$disp_fit))
  )
  
  return(pb_coeffs)
  
}

#' Classify gene expression patterns across states
#'
#' Scores and classifies genes per cell state using a state transition graph
#' and coefficient matrices from pseudobulk GLMs.
#'
#' @param state_graph State transition graph (igraph).
#' @param pb_coeffs Output from `compute_glm()`, including coefficient matrices.
#' @param pb_cds Pseudobulked CDS used for model fitting.
#' @param states_to_assess Optional vector of states to evaluate.
#' @param state_term Column/term name identifying states in the model.
#' @param log_fc_thresh Log fold-change threshold for class assignment.
#' @param abs_expr_thresh Absolute expression threshold for filtering.
#' @param cv_threshold Coefficient of variation threshold for filtering.
#' @param sig_thresh Significance threshold for class assignment.
#' @param cores Number of cores for scoring.
#'
#' @return A tibble of per-state gene class scores, nested by cell state.
classify_gene_patterns <- function(state_graph,
                                   pb_coeffs,
                                   pb_cds,
                                   states_to_assess = c(),
                                   state_term = "cell_group",
                                   log_fc_thresh = 1,
                                   abs_expr_thresh = 1e-3,
                                   cv_threshold = 100,
                                   sig_thresh = 0.05,
                                   cores = 1) {
  
  n <- dim(pb_cds)[2]
  if (length(states_to_assess) == 0) {
    states_to_assess <- unlist(igraph::V(state_graph)$name)
  } else {
    states_to_assess <- intersect(unlist(states_to_assess), unlist(igraph::V(state_graph)$name))
  }
  cell_states <- tidyr::tibble(cell_state = states_to_assess)

  cell_states <- cell_states %>%
    dplyr::mutate(gene_classes = purrr::map(
      .f = purrr::possibly(compare_genes_in_cell_state, NA_real_),
      .x = cell_state,
      state_graph = state_graph,
      estimate_matrix = pb_coeffs$coefficients,
      stderr_matrix = pb_coeffs$stdev.unscaled,
      state_term = state_term,
      log_fc_thresh = log_fc_thresh,
      abs_expr_thresh = abs_expr_thresh,
      cv_threshold = cv_threshold,
      n = n,
      sig_thresh = sig_thresh,
      cores = cores
    ))

  cell_states <- cell_states %>%
    filter(is.na(gene_classes) == FALSE) %>%
    dplyr::mutate(gene_class_scores = purrr::map2(
      .f = purrr::possibly(score_genes_for_expression_pattern, NA_real_),
      .x = cell_state,
      .y = gene_classes,
      state_graph,
      pb_coeffs$coefficients
    )) %>%
    dplyr::select(-gene_classes)

  cell_states <- cell_states %>%
    tidyr::unnest(gene_class_scores) %>%
    dplyr::left_join(
      rowData(pb_cds) %>% as.data.frame() %>% dplyr::select(id, gene_short_name),
      by = c("gene_id" = "id")
    ) %>%
    dplyr::group_by(cell_state) %>%
    tidyr::nest() %>%
    dplyr::rename(gene_class_scores = data)

  return(cell_states)
}
