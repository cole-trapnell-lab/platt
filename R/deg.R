#' This function takes a numeric vector `p` and returns a probability vector
#' where each element is divided by the sum of all elements in `p`. If the sum
#' of `p` is zero, the function returns a vector of zeros.
#'
#' @param p A numeric vector.
#'
#' @return A numeric vector of probabilities.
#'
#' @examples
#' makeprobsvec(c(1, 2, 3))
#' makeprobsvec(c(0, 0, 0))
#'
#' @noRd
makeprobsvec <- function(p) {
  phat <- p / sum(p)
  phat[is.na(phat)] <- 0
  phat
}

#' This function takes a relative abundance matrix and calculates the probability matrix.
#' It normalizes each column by dividing by the column sum, and replaces any NA values with 0.
#'
#' @param a A numeric matrix representing relative abundances.
#' @return A numeric matrix where each column sums to 1, representing probabilities.
#' @importFrom Matrix t
#' @noRd
makeprobs <- function(a) {
  colSums <- apply(a, 2, sum)
  b <- Matrix::t(Matrix::t(a) / colSums)
  b[is.na(b)] <- 0
  b
}


#' This function calculates the Shannon entropy for a given probability matrix.
#'
#' @param p A numeric matrix where each row represents a probability distribution.
#' @return A numeric vector containing the Shannon entropy for each row of the input matrix.
#' @import Matrix
#' @examples
#' p <- matrix(c(0.2, 0.3, 0.5, 0.1, 0.1, 0.8), nrow = 2, byrow = TRUE)
#' shannon_entropy(p)
#' @noRd
shannon_entropy <- function(p) {
  # if (Matrix::rowMin(p) < 0 || (p) <=0)
  #  return(Inf)
  p.norm <- p / Matrix::rowSums(p)
  lg_pnorm <- log2(p.norm) * p.norm
  lg_pnorm[p.norm == 0] <- 0
  SE <- -Matrix::rowSums(lg_pnorm)
  return(SE)
}


#' This function calculates the Jensen-Shannon distance between a given matrix of
#' probability distributions and a specified pattern. The Jensen-Shannon distance
#' is a method of measuring the similarity between two probability distributions.
#'
#' @param x A matrix where each row represents a probability distribution.
#' @param pattern A numeric vector representing the pattern to compare against.
#' @return A numeric vector representing the pattern match score for each row in the matrix.
#' @examples
#' x <- matrix(c(0.2, 0.3, 0.5, 0.1, 0.4, 0.5), nrow = 2, byrow = TRUE)
#' pattern <- c(0.3, 0.3, 0.4)
#' js_dist_to_pattern(x, pattern)
#' @noRd
js_dist_to_pattern <- function(x, pattern) {
  stopifnot(ncol(x) == length(pattern))
  avg_x_pattern <- sweep(x, 2, pattern, "+") / 2
  JSdiv <- shannon_entropy(avg_x_pattern) -
    (shannon_entropy(x) + shannon_entropy(matrix(pattern, nrow = 1))) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  pattern_match_score <- 1 - JSdist
  return(pattern_match_score)
}

#' Measure Upregulation Effect
#'
#' #' This function calculates the effect size of upregulation by comparing
#' the self estimate to the maximum of other estimates.
#'
#' @param self_estimate A numeric value representing the self estimate.
#' @param other_estimates A numeric matrix where each row represents different estimates.
#'
#' @return A numeric vector representing the effect size of upregulation.
#'
#' @importFrom matrixStats rowMaxs
#' @noRd
measure_upregulation_effect <- function(self_estimate, other_estimates) {
  effect_size <- self_estimate - matrixStats::rowMaxs(other_estimates)
  return(effect_size)
}

#' Measure Downregulation Effect
#'
#' This function calculates the downregulation effect size by comparing the
#' self estimate to the maximum of other estimates.
#'
#' @param self_estimate A numeric vector representing the self estimate.
#' @param other_estimates A numeric matrix where each row represents different
#' estimates to compare against the self estimate.
#'
#' @return A numeric vector representing the effect size, calculated as the
#' difference between the maximum of other estimates and the self estimate.
#'
#' @importFrom matrixStats rowMaxs
#' @noRd
measure_downregulation_effect <- function(self_estimate, other_estimates) {
  effect_size <- matrixStats::rowMaxs(other_estimates) - self_estimate
  return(effect_size)
}

#' Measure Maintenance Effect
#'
#' This function calculates the maintenance effect size by comparing a self estimate
#' with other estimates. The effect size is computed as the inverse of the maximum
#' absolute difference between the self estimate and other estimates.
#'
#' @param self_estimate A numeric vector representing the self estimate.
#' @param other_estimates A numeric matrix where each row represents an estimate to compare against the self estimate.
#'
#' @return A numeric vector representing the effect size for each row in the other estimates.
#'
#' @importFrom matrixStats rowMaxs
#'
#' @examples
#' self_estimate <- c(1, 2, 3)
#' other_estimates <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
#' measure_maintenance_effect(self_estimate, other_estimates)
#' @noRd
measure_maintenance_effect <- function(self_estimate, other_estimates) {
  effect_size <- 1 / matrixStats::rowMaxs(abs(other_estimates - as.vector(self_estimate)))
  return(effect_size)
}

#' Get Parent Nodes in a State Graph
#'
#' This function retrieves the parent nodes of a given cell state in a state graph.
#'
#' @param state_graph An igraph object representing the state graph.
#' @param cell_state The specific cell state for which to find parent nodes.
#'
#' @return A character vector of parent node names. If there are no parent nodes, an empty character vector is returned.
#'
#' @import igraph
#'
#' @examples
#' # Assuming `state_graph` is a pre-defined igraph object and `cell_state` is a valid node in the graph
#' parents <- get_parents(state_graph, cell_state)
#' print(parents)
#' @noRd
get_parents <- function(state_graph, cell_state) {
  parents <- igraph::neighbors(state_graph, cell_state, mode = "in")
  if (length(parents) > 0) {
    return(parents$name)
  } else {
    return(c())
  }
}

#' Get Children of a Cell State in a State Graph
#'
#' This function retrieves the children of a given cell state in a state graph.
#'
#' @param state_graph An igraph object representing the state graph.
#' @param cell_state A vertex in the state graph for which to find the children.
#' @return A character vector of the names of the children of the given cell state. If there are no children, an empty character vector is returned.
#' @import igraph
#' @noRd
get_children <- function(state_graph, cell_state) {
  children <- igraph::neighbors(state_graph, cell_state, mode = "out")
  if (length(children) > 0) {
    return(children$name)
  } else {
    return(c())
  }
}

#' Get the siblings of a state in a state transition graph
#'
#' This function retrieves the sibling states of a given cell state in a state transition graph.
#' Sibling states are defined as states that share the same parent state.
#'
#' @param state_graph An igraph object representing the state transition graph.
#' @param cell_state A character string representing the cell state for which to find siblings.
#'
#' @return A character vector of sibling states. If the cell state has no parents, an empty vector is returned.
#'
#' @import igraph
#' @noRd
get_siblings <- function(state_graph, cell_state) {
  parents <- get_parents(state_graph, cell_state)
  if (length(parents) > 0) {
    siblings <- igraph::neighbors(state_graph, parents, mode = "out")
    siblings <- setdiff(siblings$name, cell_state) # exclude self
    return(siblings)
  } else {
    return(c())
  }
}


#' Score Genes for Expression Pattern
#'
#' This function scores genes based on their expression patterns in a given cell state. It takes into account the expression estimates of the cell state, its parents, children, and siblings in a state graph.
#'
#' @param cell_state A character string representing the cell state to be analyzed.
#' @param gene_patterns A data frame containing gene patterns with gene IDs and their interpretations.
#' @param state_graph An igraph object representing the state graph.
#' @param estimate_matrix A matrix of expression estimates with genes as rows and cell states as columns.
#' @param state_term A character string representing the term used for cell groups. Default is "cell_group".
#' @param cores An integer specifying the number of cores to use for parallel processing. Default is 1.
#'
#' @return A data frame with gene patterns and their corresponding activity scores.
#'
#' @importFrom tidyr tibble unnest
#' @importFrom dplyr mutate case_when
#' @importFrom Matrix rowSums
#'
#' @noRd
score_genes_for_expression_pattern <- function(cell_state, gene_patterns, state_graph, estimate_matrix, state_term = "cell_group", cores = 1) {
  parents <- get_parents(state_graph, cell_state) # igraph::neighbors(state_graph, cell_state, mode="in")
  parents <- intersect(parents, colnames(estimate_matrix))

  children <- get_children(state_graph, cell_state) # igraph::neighbors(state_graph, cell_state, mode="out")
  children <- intersect(children, colnames(estimate_matrix))

  siblings <- get_siblings(state_graph, cell_state) # igraph::neighbors(state_graph, parents, mode="out")
  siblings <- intersect(siblings, colnames(estimate_matrix))

  # states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()

  expr_df <- tidyr::tibble(gene_id = row.names(estimate_matrix))


  self_estimates <- estimate_matrix[gene_patterns$gene_id, cell_state, drop = FALSE]
  parents_estimates <- estimate_matrix[gene_patterns$gene_id, c(parents), drop = FALSE]
  parents_and_sib_estimates <- estimate_matrix[gene_patterns$gene_id, c(parents, siblings), drop = FALSE]

  children_estimates <- estimate_matrix[gene_patterns$gene_id, c(children), drop = FALSE]

  self_and_parent <- exp(estimate_matrix[gene_patterns$gene_id, c(cell_state, parents), drop = FALSE])
  self_and_parent <- self_and_parent / Matrix::rowSums(self_and_parent) # normalize so that rows sum to 1

  self_parent_sibs <- exp(estimate_matrix[gene_patterns$gene_id, c(cell_state, parents, siblings), drop = FALSE])
  self_parent_sibs <- self_parent_sibs / Matrix::rowSums(self_parent_sibs) # normalize so that rows sum to 1


  # self_estimates = data_mat[gene_patterns$gene_id, "cell_state_shrunken_lfc", drop=FALSE] %>% as.matrix()
  # parents_estimates = data_mat[gene_patterns$gene_id, parent_cols, drop=FALSE] #%>% as.matrix()
  #
  # parents_and_sib_estimates = data_mat[gene_patterns$gene_id, c(parent_cols, sibling_cols), drop=FALSE] %>% as.matrix()
  #
  # self_and_parent = exp(data_mat[gene_patterns$gene_id, c("cell_state_shrunken_lfc", parent_cols), drop=FALSE]) %>% as.matrix()
  # self_and_parent = self_and_parent / Matrix::rowSums(self_and_parent) #normalize so that rows sum to 1
  #
  # self_parent_sibs = exp(data_mat[gene_patterns$gene_id, c("cell_state_shrunken_lfc", parent_cols, sibling_cols), drop=FALSE]) %>% as.matrix()
  # self_parent_sibs = self_parent_sibs / Matrix::rowSums(self_parent_sibs) #normalize so that rows sum to 1

  # num_parents <- length(parents)
  # num_siblings <- length(siblings)
  gene_patterns_and_scores <- gene_patterns %>%
    tidyr::unnest(interpretation) %>%
    # group_by(interpretation) %>%
    # mutate(pattern_match_score =
    #          case_when(interpretation %in% c("Absent") ~ 0,
    #                    interpretation %in% c("Maintained") ~ js_dist_to_pattern(self_and_parent,c(1, rep(1, num_parents))),
    #                    interpretation %in% c("Selectively maintained", "Specifically maintained") ~ js_dist_to_pattern(self_parent_sibs, c(1, rep(1, num_parents), rep(0, num_siblings))),
    #                    interpretation %in% c("Upregulated", "Activated") ~ js_dist_to_pattern(self_and_parent, c(1, rep(0, num_parents))),
    #                    interpretation %in% c("Selectively upregulated", "Specifically upregulated", "Selectively activated", "Specifically activated") ~ js_dist_to_pattern(self_parent_sibs, c(1, rep(0, num_parents), rep(0, num_siblings))),
    #                    interpretation %in% c("Downregulated", "Deactivated") ~ js_dist_to_pattern(self_and_parent, c(0, rep(1, num_parents))),
    #                    interpretation %in% c("Selectively downregulated", "Specifically downregulated", "Selectively deactivated", "Specifically deactivated") ~ js_dist_to_pattern(self_parent_sibs, c(0, rep(1, num_parents), rep(0, num_siblings))),
    #                    TRUE ~ 0)
    # ) %>%
    # group_by(interpretation) %>%
    mutate(
      pattern_activity_score =
        dplyr::case_when(
          interpretation %in% c("Absent") ~ 0,
          interpretation %in% c("Maintained") ~ measure_maintenance_effect(self_estimates, parents_estimates),
          interpretation %in% c("Selectively maintained", "Specifically maintained") ~ measure_maintenance_effect(self_estimates, parents_and_sib_estimates),
          interpretation %in% c("Upregulated", "Activated") ~ measure_upregulation_effect(self_estimates, parents_estimates),
          interpretation %in% c("Precursor-specific") ~ measure_upregulation_effect(self_estimates, children_estimates),
          interpretation %in% c("Precursor-depleted") ~ measure_downregulation_effect(self_estimates, children_estimates),
          interpretation %in% c("Selectively upregulated", "Specifically upregulated", "Selectively activated", "Specifically activated") ~ measure_upregulation_effect(self_estimates, parents_and_sib_estimates),
          interpretation %in% c("Downregulated", "Deactivated") ~ measure_downregulation_effect(self_estimates, parents_estimates),
          interpretation %in% c("Selectively downregulated", "Specifically downregulated", "Selectively deactivated", "Specifically deactivated") ~ measure_downregulation_effect(self_estimates, parents_and_sib_estimates),
          TRUE ~ 0
        )
    )
  return(gene_patterns_and_scores)
}

#' helper function
#' @noRd
unnest_degs <- function(ccs,
                        cell_states,
                        label_nodes_by) {
  cell_states <- cell_states %>%
    dplyr::filter(is.na(gene_classes) == FALSE) %>%
    tidyr::unnest(gene_class_scores) %>%
    dplyr::select(cell_state, gene_id, interpretation, pattern_match_score, pattern_activity_score)

  cell_states <- dplyr::left_join(cell_states,
    rowData(ccs@cds) %>% tidyr::as_tibble() %>% dplyr::select(id, gene_short_name),
    by = c("gene_id" = "id")
  )

  cell_states <- dplyr::left_join(cell_states,
    collect_psg_node_metadata(ccs, color_nodes_by = NULL, group_nodes_by = NULL, label_nodes_by = label_nodes_by),
    by = c("cell_state" = "id")
  ) %>%
    dplyr::rename(cell_type = label_nodes_by)

  return(cell_states)
}


#' Estimate Ambient RNA
#'
#' This function estimates the ambient RNA in a given dataset.
#'
#' @param ccs A cell count matrix or similar object.
#' @param abs_expr_thresh Absolute expression threshold. Default is 1e-3.
#' @param min_samples_detected Minimum number of samples in which a gene must be detected. Default is 2.
#' @param min_cells_per_pseudobulk Minimum number of cells per pseudobulk. Default is 3.
#' @param cores Number of cores to use for parallel processing. Default is 1.
#' @param gene_ids Optional vector of gene IDs to subset the data. Default is NULL.
#' @param group_nodes_by Optional column name to group nodes by. Default is NULL.
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing the ambient RNA coefficients and related statistics.
#' @export
#'
#' @examples
#' # Example usage:
#' ambient_rna <- estimate_ambient_rna(ccs)
estimate_ambient_rna <- function(ccs,
                                 abs_expr_thresh = 1e-10,
                                 min_samples_detected = 2,
                                 min_cells_per_pseudobulk = 3,
                                 cores = 1,
                                 gene_ids = NULL,
                                 group_nodes_by = NULL,
                                 ...) {
  if (is.null(group_nodes_by)) {
    pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun = "sum")
    state_term <- "cell_group"
  } else {
    pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, state_col = group_nodes_by, cell_agg_fun = "sum")
    state_term <- group_nodes_by
  }


  if (is.null(gene_ids) == FALSE) {
    pb_cds <- pb_cds[gene_ids, ]
  }


  # FIXME: assess whether we actually want to do these filtering steps:
  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  expr_over_thresh <- normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test <- which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds <- pb_cds[genes_to_test, ]

  colData(pb_cds)$Size_Factor <- colData(pb_cds)$num_cells_in_group

  # mean_expr = Matrix::rowMeans(normalized_counts(pb_cds, norm_method="size_only", pseudocount=0))
  # max_expr = DelayedMatrixStats::rowMaxs(normalized_counts(pb_cds, norm_method="size_only", pseudocount=0))
  # mean_v_max = mean_expr / max_expr
  # mean_v_max = mean_v_max[max_expr > 0]
  # median_mean_v_max = median(mean_v_max)
  # message(paste("ambient RNA estimate: ", median_mean_v_max))
  #
  # colData(pb_cds)$Size_Factor = colData(pb_cds)$Size_Factor / median_mean_v_max

  pseudobulks_to_test <- which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)

  message("fitting ambient regression models")
  pb_cds <- pb_cds[, pseudobulks_to_test]

  # Collect background estimate
  pb_ambient <- fit_models(pb_cds,
    model_formula_str = "~ 1",
    cores = cores
  )


  message("\tcollecting coefficients")
  ambient_coeffs <- collect_coefficients_for_shrinkage(pb_cds, pb_ambient, abs_expr_thresh, term_to_keep = "(Intercept)")

  rm(pb_ambient)
  gc()

  # Scale the "number of cells" up by the fraction of ambient RNA to arrive
  # at transcript count estimates coming from the soup
  # colData(pb_cds)$Size_Factor = colData(pb_cds)$Size_Factor / median_mean_v_max


  # FIXME:
  # automatic estimate of ambient fraction. Total heuristic, replace with below-the-knee estimate

  message("\testimating background")
  mean_expr <- Matrix::rowMeans(normalized_counts(pb_cds, norm_method = "size_only", pseudocount = 0))
  max_expr <- DelayedMatrixStats::rowMaxs(normalized_counts(pb_cds, norm_method = "size_only", pseudocount = 0))
  mean_v_max <- mean_expr / max_expr
  mean_v_max <- mean_v_max[max_expr > 0]
  median_mean_v_max <- median(mean_v_max)
  message(paste("ambient RNA estimate: ", median_mean_v_max))

  ambient_coeffs$coefficients <- ambient_coeffs$coefficients + log(median_mean_v_max)
  ambient_coeffs$median_mean_v_max <- median_mean_v_max
  ambient_coeffs$max_expr <- max_expr
  ambient_coeffs$mean_v_max <- mean_v_max

  # ambient_coeffs$stdev.unscaled = ambient_coeffs$stdev.unscaled

  return(ambient_coeffs)
}


#' Compare Genes Over Graph
#'
#' This function compares gene expression over a given state graph.
#'
#' @param ccs Cell count summary object.
#' @param state_graph A graph object representing the state graph.
#' @param gene_ids A vector of gene IDs to be assessed. Default is NULL.
#' @param group_nodes_by A column name to group nodes by. Default is NULL.
#' @param assembly_group A specific assembly group to assess. Default is NULL.
#' @param label_nodes_by A label for nodes. Default is "cell_state".
#' @param states_to_assess A list of states to assess. Default is an empty list.
#' @param nuisance_model_formula_str A string representing the nuisance model formula. Default is "0".
#' @param log_fc_thresh Log fold change threshold. Default is 1.
#' @param abs_expr_thresh Absolute expression threshold. Default is 1e-3.
#' @param sig_thresh Significance threshold. Default is 0.05.
#' @param min_samples_detected Minimum number of samples detected. Default is 2.
#' @param min_cells_per_pseudobulk Minimum number of cells per pseudobulk. Default is 3.
#' @param cores Number of cores to use for parallel processing. Default is 1.
#' @param cv_threshold Coefficient of variation threshold. Default is 10.
#' @param ... Additional arguments.
#'
#' @return A tibble containing gene class scores for each cell state.
#'
#' @examples
#' # Example usage:
#' # result <- compare_genes_over_graph(ccs, state_graph, gene_ids = c("gene1", "gene2"))
#'
#' @importFrom dplyr select mutate filter group_by rename left_join
#' @importFrom purrr map map2 possibly
#' @importFrom tidyr unnest nest
#' @importFrom igraph graph_from_data_frame V
#' @importFrom Matrix rowSums
#' @importFrom stringr str_replace_all
#' @importFrom tibble tibble
#'
#' @export
compare_genes_over_graph <- function(ccs,
                                     state_graph,
                                     gene_ids = NULL,
                                     group_nodes_by = NULL,
                                     assembly_group = NULL,
                                     label_nodes_by = "cell_state",
                                     states_to_assess = list(),
                                     nuisance_model_formula_str = "0",
                                     # ambient_coeffs = NULL,
                                     log_fc_thresh = 1,
                                     abs_expr_thresh = 1e-3,
                                     sig_thresh = 0.05,
                                     min_samples_detected = 2,
                                     min_cells_per_pseudobulk = 3,
                                     cores = 1,
                                     cv_threshold = 100,
                                     ...) {
  if (is.null(group_nodes_by)) {
    pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun = "sum")
    state_term <- "cell_group"
  } else {
    pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, state_col = group_nodes_by, cell_agg_fun = "sum")
    state_term <- group_nodes_by
  }

  # if we want to run it by assembly group
  if (is.null(assembly_group) == FALSE) {
    pb_cds <- hooke:::add_covariate(ccs, pb_cds, "assembly_group")
    pb_cds <- pb_cds[, colData(pb_cds)$assembly_group == assembly_group]

    if (is(state_graph, "igraph")) {
      state_graph <- igraph::as_data_frame(state_graph)
    }
    state_graph <- state_graph[state_graph$assembly_group == assembly_group, ]
  }


  if (!is(state_graph, "igraph")) {
    state_graph <- state_graph %>% igraph::graph_from_data_frame()
  }

  if (is.null(gene_ids) == FALSE) {
    pb_cds <- pb_cds[intersect(rownames(pb_cds), gene_ids), ]
  }

  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  expr_over_thresh <- monocle3::normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test <- which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds <- pb_cds[genes_to_test, ]

  pseudobulks_to_test <- which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)


  # num_cell_types = colData(pb_cds) %>% as.data.frame %>% group_by(cell_group) %>% summarise(n = sum(num_cells_in_group))
  # num_cells = colData(pb_cds) %>% as.data.frame %>%  group_by(cell_group) %>% summarise(n = sum(num_cells_in_group)) %>% pull(n)
  # names(num_cells) =  colData(pb_cds) %>% as.data.frame %>%  group_by(cell_group) %>% pull(cell_group)

  message("fitting regression models")
  pb_cds <- pb_cds[, pseudobulks_to_test]

  # subset to just the graph cells
  pb_cds <- pb_cds[, colData(pb_cds)$cell_group %in% igraph::V(state_graph)$name]

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

  if (length(states_to_assess) == 0) {
    states_to_assess <- intersect(as.character(unique(colData(pb_cds)[, state_term])), unlist(igraph::V(state_graph)$name))
  } else {
    states_to_assess <- intersect(unlist(states_to_assess), unlist(igraph::V(state_graph)$name))
  }
  cell_states <- tidyr::tibble(cell_state = states_to_assess)

  # browser()
  # I made some adjustments to this to deal with the case where there is no ambient expression, but now the join seems to be failing and there is no shrunken LFC values in the result
  # debugonce(compare_genes_in_cell_state)
  cell_states <- cell_states %>%
    dplyr::mutate(gene_classes = purrr::map(
      .f = purrr::possibly(
        compare_genes_in_cell_state, NA_real_
      ),
      .x = cell_state,
      state_graph = state_graph,
      # ambient_estimate_matrix = ambient_coeffs$coefficients,
      # ambient_stderr_matrix = ambient_coeffs$stdev.unscaled,
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
    dplyr::select(-gene_classes) # otherwise it is duplicated

  # cell_states = cell_states %>% tidyr::unnest(gene_class_scores) %>%
  #   left_join(rowData(pb_cds) %>% as.data.frame %>% dplyr::select(id, gene_short_name), by = c("gene_id" = "id")) %>%
  #   group_by(cell_state) %>% tidyr::nest() %>% dplyr::rename(gene_class_scores = data)

  cell_states <- cell_states %>%
    tidyr::unnest(gene_class_scores) %>%
    dplyr::left_join(
      rowData(pb_cds) %>% as.data.frame() %>% dplyr::select(id, gene_short_name),
      by = c("gene_id" = "id")
    ) %>%
    dplyr::group_by(cell_state) %>%
    tidyr::nest() %>%
    dplyr::rename(gene_class_scores = data)

  rm(pb_group_models)
  return(cell_states)
}


#' Collect Coefficients for Shrinkage
#'
#' This function collects coefficients for shrinkage from a given cell data set (cds) and model table (model_tbl).
#' It calculates the dispersion, mean expression, and other statistics for each model in the table.
#'
#' @param cds A cell data set object.
#' @param model_tbl A table containing models and their summaries.
#' @param abs_expr_thresh Absolute expression threshold for filtering.
#' @param term_to_keep The term to keep in the coefficient table.
#'
#' @return A tibble containing coefficients, standard deviations, sigma, degrees of freedom, trend, estimated dispersion, and fitted dispersion.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Mutates the model table to include dispersion values.
#'   \item Calculates mean expression for each model.
#'   \item Fits a generalized additive model (GAM) to the log of dispersion values.
#'   \item Updates the model summary with fitted dispersion values.
#'   \item Extracts raw coefficients and extra model statistics.
#'   \item Joins the raw coefficient table with extra model statistics.
#'   \item Filters and mutates the coefficient table based on the term to keep.
#'   \item Handles failed gene models by setting their estimates and standard errors to NA.
#'   \item Pivots the coefficient and standard error tables.
#'   \item Calculates sigma and unscaled standard deviations.
#' }
#'
#' @examples
#' \dontrun{
#' coefs_for_shrinkage <- collect_coefficients_for_shrinkage(cds, model_tbl, abs_expr_thresh, term_to_keep)
#' }
#'
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import mgcv
#' @import speedglm
#' @import tibble
#' @import stringr
#' @importFrom stats predict
#' @importFrom utils head
#' @export
collect_coefficients_for_shrinkage <- function(cds, model_tbl, abs_expr_thresh, term_to_keep) {
  model_tbl <- model_tbl %>%
    dplyr::mutate(dispersion = purrr::map2(
      .f = purrr::possibly(
        extract_dispersion_helper, NA_real_
      ), .x = model,
      .y = model_summary
    )) %>%
    tidyr::unnest(dispersion)

  mean_expr_helper <- function(model, new_data, type = "response") {
    res <- tryCatch(
      {
        # mean(stats::predict(model, newdata=new_data, type=type))
        # FIXME: this is a hack to avoid attaching speedglm for now
        mean(speedglm:::predict.speedglm(model, newdata = new_data, type = type))
      },
      error = function(e) {
        NA
      }
    )
    return(res)
  }

  model_tbl <- model_tbl %>%
    dplyr::mutate(mean_expr = purrr::map(model, mean_expr_helper, new_data = colData(cds) %>% as.data.frame())) %>%
    tidyr::unnest(mean_expr)
  model_tbl$disp_fit <- rep(1.0, nrow(model_tbl))
  tryCatch(
    {
      disp_fit <- mgcv::gam(log(dispersion) ~ s(log(mean_expr), bs = "cs"), data = model_tbl)
      model_tbl$disp_fit <- exp(predict(disp_fit))
    },
    error = function(e) {
      print(e)
      print(head(model_tbl[, c("mean_expr", "dispersion")]))
      model_tbl$disp_fit <- rep(1.0, nrow(model_tbl))
    }
  )

  # model_tbl$disp_fit = 1

  model_tbl <- update_summary(model_tbl, dispersion_type = "fitted")
  print(head(model_tbl))

  print("\tdispersions updated")

  raw_coefficient_table <- coefficient_table(model_tbl)

  # print (head(raw_coefficient_table))
  raw_coefficient_table %>%
    select(id, term, estimate) %>%
    print()


  extract_extra_model_stats <- function(model, newdata) {
    if (class(model)[1] == "speedglm") {
      extra_stats <- tibble(RSS = model$RSS, df.residual = model$df, mean_expr = mean(predict(model, newdata = newdata)))
      return(extra_stats)
    } else {
      extra_stats <- tibble(RSS = NA_real_, df.residual = NA_real_, mean_expr = NA_real_)
      return(extra_stats)
    }
  }

  print("\tcomputing extra model stats")

  extra_model_stats <- model_tbl %>%
    dplyr::mutate(extra_stats = purrr::map(.f = purrr::possibly(
      extract_extra_model_stats, NA_real_
    ), .x = model, newdata = colData(cds) %>% as.data.frame())) %>%
    dplyr::select(id, extra_stats) %>%
    tidyr::unnest(extra_stats)

  raw_coefficient_table <- left_join(raw_coefficient_table, extra_model_stats %>% dplyr::select(id, RSS, df.residual), by = "id") %>%
    # left_join(model_tbl %>% select(id, disp_fit, dispersion), by = "id") %>%
    # coefficient_table(pb_group_models) %>%
    # dplyr::select(gene_short_name, id, term, estimate, std_err, p_value, status) %>%
    filter(grepl(term_to_keep, term)) %>%
    mutate(term = stringr::str_replace_all(term, term_to_keep, "")) %>%
    mutate(term = stringr::str_replace(term, "\\(\\)", "Intercept"))

  fail_gene_ids <- model_tbl %>%
    filter(status == "FAIL") %>%
    select(id)
  term_to_keep_levels <- levels(as.factor(colData(cds)[, term_to_keep]))
  # term_table = tibble(term=stringr::str_c(term_to_keep, term_to_keep_levels))
  term_table <- tibble(term = term_to_keep_levels)

  fail_coeff_rows <- tidyr::crossing(fail_gene_ids, term_table) %>% mutate(estimate = NA_real_)
  raw_coefficient_table <- raw_coefficient_table %>%
    mutate(estimate = if_else(id %in% fail_gene_ids$id, NA_real_, estimate))
  # raw_coefficient_table = bind_rows(raw_coefficient_table, fail_coeff_rows)

  estimate_matrix <- raw_coefficient_table %>% dplyr::select(id, term, estimate)
  if (term_to_keep != "(Intercept)") {
    estimate_matrix <- estimate_matrix %>% mutate(term = factor(term, levels = term_to_keep_levels))
  }
  estimate_matrix <- estimate_matrix %>% mutate(id = factor(id, levels = model_tbl$id))

  print("\tpivoting coefficient table:")
  # print(head(estimate_matrix))
  estimate_matrix <- estimate_matrix %>% tidyr::pivot_wider(names_from = term, values_from = estimate, id_expand = TRUE, values_fill = 0)

  gene_ids <- estimate_matrix$id
  estimate_matrix$id <- NULL
  estimate_matrix <- as.matrix(estimate_matrix)
  row.names(estimate_matrix) <- gene_ids
  colnames(estimate_matrix) <- as.character(colnames(estimate_matrix))

  stderr_matrix <- raw_coefficient_table %>% dplyr::select(id, term, std_err)
  if (term_to_keep != "(Intercept)") {
    stderr_matrix <- stderr_matrix %>% mutate(term = factor(term, levels = term_to_keep_levels))
  }
  stderr_matrix <- stderr_matrix %>% mutate(id = factor(id, levels = model_tbl$id))

  print("\tpivoting coefficient stderr table")
  stderr_matrix <- stderr_matrix %>% tidyr::pivot_wider(names_from = term, values_from = std_err, id_expand = TRUE, values_fill = 0)

  gene_ids <- stderr_matrix$id
  stderr_matrix$id <- NULL
  stderr_matrix <- as.matrix(stderr_matrix)
  row.names(stderr_matrix) <- gene_ids
  colnames(stderr_matrix) <- as.character(colnames(stderr_matrix))

  # collect the ids of any genes that threw an exception in fit_models and
  # set their estimates and std_errors to NA

  if (length(fail_gene_ids$id) > 0) {
    print("\tzero-ing out failed gene models")
    estimate_matrix[fail_gene_ids$id, ] <- log(abs_expr_thresh)
    stderr_matrix[fail_gene_ids$id, ] <- Inf
  }


  sigma_df <- raw_coefficient_table %>%
    dplyr::select(id, RSS, df.residual, mean_expr, disp_fit, dispersion) %>%
    distinct() %>%
    as.data.frame()
  row.names(sigma_df) <- sigma_df$id
  sigma_df$id <- NULL
  sigma_df <- sigma_df[row.names(estimate_matrix), ]
  sigma_df$sigma <- sigma_df$RSS / sigma_df$df.residual

  std_dev.unscaled <- stderr_matrix^2 / sigma_df$sigma

  print("estimates")
  print(head(estimate_matrix))

  print("std errs")
  print(head(stderr_matrix))

  print("\treporting final stats")
  coefs_for_shrinkage <- tibble(
    coefficients = estimate_matrix,
    stdev.unscaled = stderr_matrix,
    sigma = sigma_df$sigma,
    df.residual = sigma_df$df.residual,
    trend = sigma_df$mean_expr,
    est_dispersion = sigma_df$dispersion,
    disp_fit = sigma_df$disp_fit
  )

  return(coefs_for_shrinkage)
}

#' Compare Gene Expression Within State Graph
#'
#' This function compares gene expression within a state graph for a given cell cluster set (ccs).
#'
#' @param ccs A cell cluster set object.
#' @param perturbation_col The column name in the cell data set that contains perturbation information. Default is "perturbation".
#' @param control_ids A vector of control IDs. Default is c("Control").
#' @param nuisance_model_formula_str A string representing the nuisance model formula. Default is "0".
#' @param ambient_coeffs Coefficients for ambient noise. Default is NULL.
#' @param cell_groups A vector of cell groups to be analyzed. Default is NULL.
#' @param assembly_group The assembly group to be analyzed. Default is NULL.
#' @param state_graph A state graph object. Default is NULL.
#' @param perturbations A vector of perturbations to be analyzed. Default is NULL.
#' @param gene_ids A vector of gene IDs to be analyzed. Default is NULL.
#' @param group_nodes_by A parameter to group nodes by. Default is NULL.
#' @param log_fc_thresh Log fold change threshold. Default is 1.
#' @param abs_expr_thresh Absolute expression threshold. Default is 1e-3.
#' @param sig_thresh Significance threshold. Default is 0.05.
#' @param min_samples_detected Minimum number of samples detected. Default is 2.
#' @param min_cells_per_pseudobulk Minimum number of cells per pseudobulk. Default is NULL.
#' @param cores Number of cores to use for parallel processing. Default is 1.
#' @param write_dir Directory to write output files. Default is NULL.
#' @param max_simultaneous_genes Maximum number of genes to analyze simultaneously. Default is NULL.
#' @param ... Additional arguments.
#'
#' @return A data frame containing the comparison results.
#'
#' @examples
#' # Example usage:
#' compare_genes_within_state_graph(ccs, perturbation_col = "perturbation", control_ids = c("Control"))
#'
#' @export
compare_genes_within_state_graph <- function(ccs,
                                             perturbation_col = "perturbation",
                                             control_ids = c("Control"),
                                             nuisance_model_formula_str = "0",
                                             ambient_coeffs = NULL,
                                             cell_groups = NULL,
                                             assembly_group = NULL,
                                             state_graph = NULL,
                                             perturbations = NULL,
                                             group_nodes_by = NULL,
                                             log_fc_thresh = 1,
                                             abs_expr_thresh = 1e-3,
                                             sig_thresh = 0.05,
                                             min_samples_detected = 2,
                                             min_cells_per_pseudobulk = NULL,
                                             cores = 1,
                                             write_dir = NULL,
                                             max_simultaneous_genes = NULL,
                                             cv_threshold = 100,
                                             ...) {
  # to do make sure that ccs and state graph match
  # check if all nodes in the state graph exist in the cds

  # if (is.null(state_graph) == FALSE) {
  #
  # }

  expts <- unique(colData(ccs)$expt)

  pb_cds <- hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun = "sum")
  colData(pb_cds)$Size_Factor <- colData(pb_cds)$num_cells_in_group

  # Once we pseudobulks, we are going to overwrite the provided control ids with "Control"
  # pb_cds = hooke:::add_covariate(ccs, pb_cds, perturbation_col)
  colData(pb_cds)[["perturbation"]] <- colData(pb_cds)[[perturbation_col]]
  colData(pb_cds)[["perturbation"]] <- ifelse(colData(pb_cds)$perturbation %in% control_ids,
    "Control", colData(pb_cds)$perturbation
  )

  # to do: if we want to run it by assembly group, must provide a state graph
  if (is.null(assembly_group) == FALSE & is.null(state_graph) == FALSE) {
    pb_cds <- hooke:::add_covariate(ccs, pb_cds, "assembly_group")
    pb_cds <- pb_cds[, colData(pb_cds)$assembly_group == assembly_group]

    if (is(state_graph, "igraph")) {
      state_graph <- igraph::as_data_frame(state_graph)
    }
    state_graph <- state_graph %>% filter(assembly_group == assembly_group)
  }

  # If the user provided a set of perturbations, subset the pseudobulk CDS to include just them and the controls
  if (!is.null(perturbations)) {
    perturbations <- setdiff(perturbations, "Control")
    pb_cds <- pb_cds[, colData(pb_cds)[[perturbation_col]] %in% c("Control", perturbations)]
  } else { # otherwise, grab the perturbation ids from the CDS and process them all
    perturbations <- unique(colData(pb_cds)[["perturbation"]])
    print(paste("pre-filter controls:", "Control"))
    print(paste("pre-filter", perturbations))
    perturbations <- setdiff(perturbations, "Control")
    print(paste("post-filter", perturbations))
  }

  # subset to genes that are expressed over a certain min value
  expr_over_thresh <- normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test <- which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds <- pb_cds[genes_to_test, ]



  # # Collect background estimate
  # pb_ambient = fit_models(pb_cds,
  #                         model_formula_str="~ 1",
  #                         cores=cores)
  # ambient_coeffs = collect_coefficients_for_shrinkage(pb_cds, pb_ambient, abs_expr_thresh, term_to_keep = "(Intercept)")


  if (is.null(cell_groups)) {
    cell_groups <- rownames(counts(ccs))
  }

  # if (is.null(ambient_coeffs)){
  #   ambient_estimate_matrix = NULL
  #   ambient_stderr_matrix =  NULL
  # } else {
  #   ambient_estimate_matrix = ambient_coeffs$coefficients
  #   ambient_stderr_matrix = ambient_coeffs$stdev.unscaled
  # }

  if (is.null(write_dir) == FALSE) {
    if (!file.exists(write_dir)) {
      dir.create(write_dir)
    }
  }

  df <- data.frame(cell_group = cell_groups) %>%
    mutate(genes_within_cell_group = purrr::map(
      .f = purrr:::possibly(compare_gene_expression_within_node, NA_character_),
      .x = cell_group,
      pb_cds = pb_cds,
      ccs = ccs,
      control_ids = c("Control"),
      perturbation_ids = perturbations,
      ambient_estimate_matrix = ambient_coeffs$coefficients,
      ambient_stderr_matrix = ambient_coeffs$stdev.unscaled,
      cores = cores,
      max_simultaneous_genes = max_simultaneous_genes,
      nuisance_model_formula_str = nuisance_model_formula_str,
      min_cells_per_pseudobulk = min_cells_per_pseudobulk,
      write_dir = write_dir, 
      cv_threshold = cv_threshold, 
      abs_expr_thresh = abs_expr_thresh
    ))

  return(df)
}

#' Compare Gene Expression Within Node
#'
#' This function compares gene expression within a specified cell group,
#' fitting models per cell group and computing contrasts for perturbations
#' against controls. It can also compare perturbations to ambient RNA levels
#' if provided.
#'
#' @param cell_group Character. The cell group to analyze.
#' @param ccs Data frame or matrix. Cell cycle scores.
#' @param pb_cds CellDataSet. The pseudobulk cell dataset.
#' @param control_ids Character vector. IDs of control samples.
#' @param perturbation_ids Character vector. IDs of perturbation samples.
#' @param ambient_estimate_matrix Matrix. Ambient RNA estimate matrix.
#' @param ambient_stderr_matrix Matrix. Ambient RNA standard error matrix.
#' @param nuisance_model_formula_str Character. Formula string for nuisance model.
#' @param state_term Character. The term in colData to use for state.
#' @param log_fc_thresh Numeric. Log fold change threshold.
#' @param abs_expr_thresh Numeric. Absolute expression threshold.
#' @param sig_thresh Numeric. Significance threshold.
#' @param min_cells_per_pseudobulk Integer. Minimum cells per pseudobulk.
#' @param exclude_results_below_ambient Logical. Whether to exclude results below ambient RNA levels.
#' @param expected_effect_mode_interval Numeric vector. Interval for expected effect mode.
#' @param write_dir Character. Directory to write output files.
#' @param cores Integer. Number of cores to use for parallel processing.
#' @param max_simultaneous_genes Integer. Maximum number of genes to process simultaneously.
#'
#' @return A data frame with perturbation effects, or NULL if write_dir is specified.
#'
#' @examples
#' compare_gene_expression_within_node(
#'   cell_group = "group1",
#'   ccs = ccs_data,
#'   pb_cds = pb_cds_data,
#'   control_ids = c("ctrl1", "ctrl2"),
#'   perturbation_ids = c("pert1", "pert2"),
#'   write_dir = "output"
#' )
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr select group_by summarize mutate filter left_join
#' @importFrom tibble tibble
#' @importFrom monocle3 aggregate_gene_expression
#' @importFrom Matrix colSums rowSums
#' @importFrom stringr str_replace_all
#' @importFrom purrr map possibly
#' @importFrom tidyr unnest nest
#' @importFrom utils write.csv
#' @export
compare_gene_expression_within_node <- function(cell_group,
                                                ccs,
                                                pb_cds,
                                                control_ids,
                                                perturbation_ids,
                                                ambient_estimate_matrix = NULL,
                                                ambient_stderr_matrix = NULL,
                                                nuisance_model_formula_str = "0",
                                                state_term = "cell_group",
                                                log_fc_thresh = 1,
                                                abs_expr_thresh = 1e-3,
                                                sig_thresh = 0.05,
                                                min_cells_per_pseudobulk = NULL,
                                                exclude_results_below_ambient = TRUE,
                                                expected_effect_mode_interval = c(-10, 10),
                                                write_dir = NULL,
                                                cores = 1,
                                                max_simultaneous_genes = NULL,
                                                cv_threshold = NULL) {
  # now fit models per cell group

  cg_pb_cds <- pb_cds[, colData(pb_cds)[[state_term]] == cell_group]

  assertthat::assert_that(is.na(cell_group) == FALSE)

  # Set the factor levels to include all provided controls and all provided
  # perturbation ids so we get coefficient entries for all of them, even if no cells are present in some groups
  colData(cg_pb_cds)$perturbation <- factor(colData(cg_pb_cds)$perturbation, levels = c(control_ids, perturbation_ids))

  # if (is.null(min_cells_per_pseudobulk)){
  perturb_sf_summary <- colData(cg_pb_cds) %>%
    as.data.frame() %>%
    dplyr::select(perturbation, Size_Factor) %>%
    dplyr::group_by(perturbation) %>%
    dplyr::summarize(mean_log_sf = mean(log(Size_Factor)))

  mean_ctrl_sf <- perturb_sf_summary %>%
    dplyr::filter(perturbation %in% control_ids) %>%
    dplyr::pull(mean_log_sf) %>%
    mean()

  perturb_sf_summary <- perturb_sf_summary %>% mutate(ctrl_log_sf = mean_ctrl_sf)

  pb_group_by_term <- tibble::tibble(
    pb_id = row.names(colData(cg_pb_cds)),
    pb_group = colData(cg_pb_cds)$perturbation
  )

  mean_expr_mat <- monocle3::aggregate_gene_expression(cg_pb_cds,
    cell_group_df = pb_group_by_term,
    norm_method = "size_only",
    pseudocount = 0,
    scale_agg_values = FALSE
  )

  # if (is.null(ambient_estimate_matrix) == FALSE){
  #   mean_expr_mat = log(mean_expr_mat)
  #   above_ambient = mean_expr_mat > ambient_estimate_matrix[row.names(mean_expr_mat),1]
  #   above_ambient_genes = Matrix::rowSums(above_ambient) > 0
  # }

  detected_genes <- Matrix::colSums(counts(cg_pb_cds) > 0)
  detected_genes <- as.data.frame(detected_genes)
  detected_genes$perturbation <- colData(cg_pb_cds)$perturbation
  detected_genes <- detected_genes %>%
    dplyr::group_by(perturbation) %>%
    dplyr::summarize(detected_genes = mean(detected_genes))

  mean_ctrl_detected_genes <- detected_genes %>%
    dplyr::filter(perturbation %in% control_ids) %>%
    dplyr::pull(detected_genes) %>%
    mean()

  detected_genes <- detected_genes %>% mutate(ctrl_detected_genes = mean_ctrl_detected_genes)

  perturb_sf_summary <- perturb_sf_summary %>% left_join(detected_genes)

  nz_genes <- Matrix::rowSums(counts(cg_pb_cds)) > 0
  cg_pb_cds <- cg_pb_cds[nz_genes, ]

  message(paste0("fitting ", nrow(cg_pb_cds), " regression models for ", cell_group))

  colData(cg_pb_cds)$Size_Factor <- colData(cg_pb_cds)$num_cells_in_group

  full_model_str <- paste0("~ 0 + perturbation")
  nuisance_model_formula_str <- stringr::str_replace_all(nuisance_model_formula_str, "~", "")

  if (nuisance_model_formula_str != "0") {
    full_model_str <- paste0(full_model_str, " + ", nuisance_model_formula_str)
  }

  if (is.null(max_simultaneous_genes)) {
    max_simultaneous_genes <- nrow(cg_pb_cds)
  }
  message(paste0("\tregressing with formula ", full_model_str))
  message(paste0("\tUsing max_simultaneous_genes = ", max_simultaneous_genes))

  # fit_model_blocks = ceiling(rcg_pb_cds) / max_simultaneous_genes)
  # print (ceiling(seq_along(row.names(cg_pb_cds))/max_simultaneous_genes))

  gene_blocks <- split(row.names(cg_pb_cds), ceiling(seq_along(row.names(cg_pb_cds)) / max_simultaneous_genes))
  # print (head(gene_blocks))
  print("levels for contrast:")
  print(levels(colData(cg_pb_cds)$perturbation))

  coeffs_for_blocks <- lapply(gene_blocks, function(gb) {
    gb <- unlist(gb)
    # print (gb)
    gb_cds <- cg_pb_cds[gb, ]
    # print (gb_cds)
    pb_group_models <- fit_models(gb_cds,
      model_formula_str = full_model_str,
      cores = cores
    ) %>% dplyr::select(gene_short_name, id, model, model_summary, status)
    
    # pb_coeffs = collect_coefficients_for_limma(cg_pb_cds, pb_group_models, abs_expr_thresh) #coefficient_table(pb_group_models) %>%

    message("\tcollecting coefficients")

    pb_coeffs <- collect_coefficients_for_shrinkage(gb_cds, pb_group_models, abs_expr_thresh, term_to_keep = "perturbation") # coefficient_table(pb_group_models) %>%

    rm(pb_group_models) # DO NOT REMOVE. This is important for keeping the memory footprint of this analysis light.
    gc()

    return(pb_coeffs)
  })

  pb_coeffs <- list()
  pb_coeffs$coefficients <- do.call(rbind, lapply(coeffs_for_blocks, function(b) {
    b$coefficients
  }))
  pb_coeffs$stdev.unscaled <- do.call(rbind, lapply(coeffs_for_blocks, function(b) {
    b$stdev.unscaled
  }))
  # pb_coeffs$sigma = do.call(rbind, lapply(coeffs_for_blocks, function(b){b$sigma}))
  pb_coeffs$sigma <- do.call(c, lapply(coeffs_for_blocks, function(b) {
    b$sigma
  }))
  pb_coeffs$df.residual <- do.call(c, lapply(coeffs_for_blocks, function(b) {
    b$df.residual
  }))
  pb_coeffs$trend <- do.call(c, lapply(coeffs_for_blocks, function(b) {
    b$trend
  }))
  pb_coeffs$est_dispersion <- do.call(c, lapply(coeffs_for_blocks, function(b) {
    b$est_dispersion
  }))
  pb_coeffs$disp_fit <- do.call(c, lapply(coeffs_for_blocks, function(b) {
    b$disp_fit
  }))

  assertthat::assert_that(nrow(pb_coeffs$coefficients) == nrow(cg_pb_cds))

  # estimates_for_controls = Matrix::rowSums(pb_coeffs$coefficients[,control_ids, drop=F])
  # stderrs_for_controls = sqrt(Matrix::rowSums(pb_coeffs$stdev.unscaled[,control_ids, drop=F]^2))

  # perturbation_ids = setdiff(colnames(estimate_matrix))
  cell_perturbations <- tidyr::tibble(term = unique(perturbation_ids))

  cell_perturbations <- left_join(cell_perturbations, perturb_sf_summary, by = c("term" = "perturbation"))

  message(paste("\tcomputing contrasts for", unique(perturbation_ids)))
  cell_perturbations <- cell_perturbations %>%
    mutate(perturb_effects = purrr:::map(
      .f = purrr::possibly(contrast_helper, NA_character_),
      .x = term,
      state_2 = control_ids,
      PEM = pb_coeffs$coefficients,
      PSEM = pb_coeffs$stdev.unscaled,
      n = dim(cg_pb_cds)[2], 
      prefix = "perturb_to_ctrl",
      ash.control = list(mode = expected_effect_mode_interval), 
      cv_threshold = cv_threshold, 
      abs_expr_thresh = abs_expr_thresh
    ))

  gene_map <- rowData(cg_pb_cds) %>%
    as.data.frame() %>%
    dplyr::select(gene_short_name, id) %>%
    dplyr::distinct()

  cell_perturbations <- cell_perturbations %>%
    dplyr::filter(!is.na(perturb_effects)) %>%
    tidyr::unnest(perturb_effects) %>%
    dplyr::left_join(gene_map, by = "id") %>%
    dplyr::group_by(term) %>%
    tidyr::nest("perturb_effects" = -term)


  message(paste("\tcontrasts complete for", unique(perturbation_ids)))

  # cell_perturbations = cell_perturbations %>%
  #   tidyr::unnest(perturb_effects) %>%
  #   left_join(rowData(pb_cds) %>% as.data.frame %>% select(id, gene_short_name), by = c("id" = "id"))  %>%
  #   group_by(term) %>% tidyr::nest() %>% dplyr::rename(perturb_effects = data)



  if (exclude_results_below_ambient & is.null(ambient_estimate_matrix) == FALSE) {
    message("\tfiltering contrast for ambient RNA")
    ambient_res <- cell_perturbations %>%
      tidyr::unnest(ambient_effects) %>%
      dplyr::filter((ctrl_to_ambient_p_value < 0.05 & ctrl_to_ambient_shrunken_lfc > 0) |
        (perturb_to_ambient_p_value < 0.05 & perturb_to_ambient_shrunken_lfc > 0)) %>%
      dplyr::select(term, id)

    # perturb_res = left_join(cell_perturbations %>%
    #                         dplyr::select(-ambient_effects) %>%
    #                         tidyr::unnest(perturb_effects),
    #                         ambient_res, by=c("term", "id"))

    perturb_res <- left_join(ambient_res,
      cell_perturbations %>%
        dplyr::select(-ambient_effects) %>%
        tidyr::unnest(perturb_effects),
      by = c("term", "id")
    )

    perturb_res <- perturb_res %>%
      dplyr::left_join(rowData(pb_cds) %>% as.data.frame() %>% dplyr::select(id, gene_short_name), by = c("id" = "id"))

    # cell_perturbations = perturb_res %>%
    #  group_by(term) %>% tidyr::nest()
  }

  if (is.null(write_dir) == FALSE) {
    message(paste("\twriting output to", write_dir))
    # message(head(cell_perturbations))

    cell_group_no_spaces <- gsub("[[:punct:]]", "", cell_group)
    cell_group_no_spaces <- gsub(" ", "_", cell_group_no_spaces)
    print(colnames(cell_perturbations))
    deg_out <- cell_perturbations %>% tidyr::unnest(perturb_effects)

    write.csv(deg_out %>% mutate(cell_group = cell_group),
      file = paste0(write_dir, "/", cell_group, "_within_node_degs.csv"),
      row.names = FALSE
    )
    return(NULL)
  } else {
    
    return(cell_perturbations)
  }
}



#' @noRd
extract_dispersion_helper <- function(model, model_summary) {
  disp_res <- if (class(model)[1] == "speedglm") {
    model_summary$dispersion
  } else if (class(model)[1] == "negbin") {
    model_summary$dispersion
  } else if (class(model) == "zeroinfl") {
    model_summary$dispersion
  } else {
    model_summary$dispersion
  }
}

#' @noRd
update_summary <- function(model_tbl, dispersion_type = c("max", "fitted", "estimated"), min_dispersion = 1) {
  dispersion_type <- match.arg(dispersion_type)
  model_tbl <- model_tbl %>% mutate(
    disp_for_test = dplyr::case_when(
      dispersion_type == "estimated" ~ dispersion,
      dispersion_type == "fitted" ~ disp_fit,
      dispersion_type == "max" ~ pmax(disp_fit, dispersion, min_dispersion),
      .default = dispersion
    ),
    model_summary = purrr::map2(model, .y = disp_for_test, .f = summary)
  )
  return(model_tbl)
}


#' Calculate GSEA Enrichment on State-Specific Genes
#'
#' This function calculates Gene Set Enrichment Analysis (GSEA) enrichment scores
#' for state-specific genes based on their pattern activity scores.
#'
#' @param gene_patterns_within_state_graph A data frame containing gene pattern activity scores.
#' @param gene_df A data frame containing gene information with columns `gene_short_name` and `gs_name`.
#' @param sig_thresh A numeric value specifying the significance threshold for adjusted p-values. Default is 0.1.
#'
#' @return A tibble containing the GSEA results filtered by the specified significance threshold.
#'
#' @examples
#' \dontrun{
#' gene_patterns_within_state_graph <- data.frame(
#'   pattern_activity_score = matrix(runif(100), ncol = 1),
#'   gene_short_name = sample(letters, 100, replace = TRUE)
#' )
#' gene_df <- data.frame(
#'   gene_short_name = sample(letters, 100, replace = TRUE),
#'   gs_name = sample(LETTERS[1:5], 100, replace = TRUE)
#' )
#' result <- calc_gsea_enrichment_on_state_specific_genes(gene_patterns_within_state_graph, gene_df)
#' }
#' @export
calc_gsea_enrichment_on_state_specific_genes <- function(gene_patterns_within_state_graph,
                                                         gene_df,
                                                         sig_thresh = 0.1) {
  gene_set_list <- split(x = gene_df$gene_short_name, f = gene_df$gs_name)
  gene_ranking <- gene_patterns_within_state_graph$pattern_activity_score[, 1]
  names(gene_ranking) <- gene_patterns_within_state_graph %>% pull(gene_short_name)
  gsea_res <- fgsea::fgsea(pathways = gene_set_list, stats = gene_ranking) %>% as_tibble()
  gsea_res <- gsea_res %>% filter(padj < sig_thresh)
  return(gsea_res)
}

#' Calculate Pathway Enrichment on State-Specific Genes
#'
#' This function calculates pathway enrichment for a given set of state-specific genes.
#'
#' @param gene_df A data frame containing gene information, with at least one column named `gene_short_name` that contains gene symbols.
#' @param msigdbr_t2g A data frame containing the TERM2GENE mapping from the MSigDB database.
#' @param sig_thresh A numeric value specifying the significance threshold for enrichment results. Default is 0.1.
#' @param ... Additional arguments passed to the `clusterProfiler::enricher` function.
#'
#' @return A tibble containing the enrichment results.
#' @importFrom dplyr %>%
#' @export
calc_pathway_enrichment_on_state_specific_genes <- function(gene_df, msigdbr_t2g, sig_thresh = 0.1, ...) {
  gene_symbols_vector <- gene_df$gene_short_name
  enrich_res <- clusterProfiler::enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...) %>% as_tibble()
  return(enrich_res)
}


#' Function to find differentially expressed genes (DEGs) within a cell type
#'
#' This function identifies DEGs within a specified cell type by fitting models to pseudobulked cell data.
#'
#' @param ccm A list containing the output from `fit_genotype_ccm`.
#' @param cores An integer specifying the number of cores to use for parallel processing. Default is 1.
#' @param cell_agg_fun A character string specifying the function to use for cell aggregation. Options are "sum" or "mean". Default is "sum".
#' @param ... Additional filtering criteria for subsetting cell types.
#'
#' @return A list of fitted models for each gene.
#'
#' @examples
#' # Example usage:
#' # ccm <- fit_genotype_ccm(...)
#' # degs <- fit_genotype_deg(ccm, cores = 2, cell_agg_fun = "mean", cell_type = "T cells")
#'
#' @export
fit_genotype_deg <- function(ccm,
                             cores = 1,
                             cell_agg_fun = c("sum", "mean"),
                             ...) {
  cell_agg_fun <- match.arg(cell_agg_fun)

  # filter by cell types
  sub_ccs <- hooke:::subset_ccs(ccm@ccs, ...)
  sub_pb_cds <- hooke::pseudobulk_ccs_for_states(sub_ccs, cell_agg_fun = cell_agg_fun)
  sub_pb_cds <- add_covariate(
    sub_ccs,
    sub_pb_cds,
    ccm@info$perturbation_col
  )
  # filter by perturbation
  sub_pb_cds <- sub_pb_cds[, colData(sub_pb_cds)[[ccm@info$perturbation_col]] %in% c(ccm@info$genotype, ccm@info$ctrl_ids)]
  colData(sub_pb_cds)$knockout <- ifelse(colData(sub_pb_cds)[[ccm@info$perturbation_col]] %in% ccm@info$ctrl_ids, F, T)

  gene_fits <- fit_models(sub_pb_cds,
    model_formula_str = "~ knockout",
    weights = colData(sub_pb_cds)$num_cells_in_group,
    cores = cores
  )

  return(gene_fits)
}



# state 1 - state 2
# is state 1 higher
#
#' Contrast Helper Function
#'
#' This function calculates the contrast between two states using provided effect matrices and standard error matrices.
#'
#' @param state_1 A vector of column indices or names for the first state in the PEM and PSEM matrices.
#' @param state_2 A vector of column indices or names for the second state in the PEM and PSEM matrices.
#' @param PEM A matrix of posterior effect means.
#' @param PSEM A matrix of posterior standard errors of the means.
#' @param PEM_2 An optional second matrix of posterior effect means. Defaults to PEM.
#' @param PSEM_2 An optional second matrix of posterior standard errors of the means. Defaults to PSEM.
#' @param prefix An optional prefix to add to the column names of the resulting data frame.
#' @param ash.control A list of control parameters for the `ashr::ash` function, including `mixcompdist` and `mode`.
#' @param cv_threshold A numeric value specifying the coefficient of variation threshold. Defaults to 10.
#'
#' @return A tibble containing the contrast results, including raw and shrunken log fold changes, standard errors, p-values, and other metrics.
#'
#' @importFrom Matrix rowSums
#' @importFrom ashr ash
#' @importFrom tibble tibble
#' @importFrom stats pnorm
#'
#' @examples
#' # Example usage:
#' contrast_res <- contrast_helper(
#'   state_1 = c("state1_col1", "state1_col2"),
#'   state_2 = c("state2_col1", "state2_col2"),
#'   PEM = pem_matrix,
#'   PSEM = psem_matrix
#' )
contrast_helper <- function(state_1,
                            state_2,
                            PEM,
                            PSEM,
                            PEM_2 = PEM,
                            PSEM_2 = PSEM,
                            prefix = NULL,
                            n = NULL, 
                            ash.control = NULL,
                            cv_threshold = NULL, 
                            abs_expr_thresh = 1e-3) {
  
  ash.mixcompdist <- "uniform"
  coefficient_mode <- 0

  if (is.null(ash.control) == FALSE) {
    if (is.null(ash.control[["mixcompdist"]]) == FALSE) {
      ash.mixcompdist <- ash.control[["mixcompdist"]]
    }
    if (is.null(ash.control[["mode"]]) == FALSE) {
      coefficient_mode <- ash.control[["mode"]]
    }
  }

  ids <- intersect(rownames(PEM), rownames(PEM_2))

  state_1_effects <- Matrix::rowSums(PEM[, state_1, drop = F])
  state_1_effects_se <- sqrt(Matrix::rowSums(PSEM[, state_1, drop = F]^2))

  state_2_effects <- Matrix::rowSums(PEM_2[, state_2, drop = F])
  state_2_effects_se <- sqrt(Matrix::rowSums(PSEM_2[, state_2, drop = F]^2))

  log_mean_expression <- log(exp(state_1_effects[ids]) + exp(state_2_effects[ids]))

  effect_est <- state_1_effects[ids] - state_2_effects[ids]
  
  # if the CV is big, take the smaller standard error
  if (is.null(cv_threshold) == FALSE) {
    
    cv_state_1 <- state_1_effects_se[ids] * sqrt(n) / state_1_effects[ids]
    cv_state_2 <- state_2_effects_se[ids] * sqrt(n) / state_2_effects[ids]
    
    # get the max of each id
    cv_max <- pmax(abs(cv_state_1), abs(cv_state_2))
    
    reliable_state_1_effects = abs(cv_state_1) < cv_threshold
    reliable_state_2_effects = abs(cv_state_2) < cv_threshold
    
    # what is the smallest nonzero expression value in each condition
    min_state_1_effect = min(state_1_effects[reliable_state_1_effects])
    min_state_1_effect_se = mean(state_1_effects_se[which(state_1_effects == min_state_1_effect)])
    
    min_state_2_effect = min(state_2_effects[reliable_state_2_effects])
    min_state_2_effect_se = mean(state_2_effects_se[which(state_2_effects == min_state_2_effect)])
    
    # set any gene that has a lower value to that determined threshold value
    state_1_effects_thresholded = pmax(state_1_effects, min_state_1_effect)
    state_2_effects_thresholded = pmax(state_2_effects, min_state_2_effect)
    
    effect_est = state_1_effects_thresholded[ids] - state_2_effects_thresholded[ids]
    
    # does it make sense to use 
    # things that have been cv thresholded have been capped, take the smaller std error
    se_est <- ifelse(cv_max > cv_threshold,
                     pmin(state_1_effects_se[ids], state_2_effects_se[ids]),
                     sqrt(state_1_effects_se[ids]^2 + state_2_effects_se[ids]^2))
    
    
  } else {
    se_est <- sqrt(state_1_effects_se[ids]^2 + state_2_effects_se[ids]^2)
  }

  up_down_large_effect_skew <- NA
  effect_skewness_classic <- NA
  if (length(coefficient_mode) >= 2) {
    extreme_effect_upper_thresh <- max(coefficient_mode[1:2])
    extreme_effect_lower_thresh <- min(coefficient_mode[1:2])

    num_large_up_effects <- sum(effect_est > extreme_effect_upper_thresh)
    num_large_down_effects <- sum(effect_est < extreme_effect_upper_thresh)
    up_down_large_effect_skew <- log(num_large_up_effects / num_large_down_effects)
    # if (abs(up_down_large_effect_skew) > log(2)){
    #  ash.mixcompdist = "halfuniform"
    # }

    numerical_mode <- function(x) {
      median(x)
    }

    effect_mode_range <- effect_est[effect_est > extreme_effect_lower_thresh & effect_est < extreme_effect_upper_thresh]
    coefficient_mode <- numerical_mode(effect_mode_range)

    effect_est <- effect_est - coefficient_mode
  }

  shrunken_res <- tryCatch(
    {
      ashr::ash(effect_est, se_est, method = "fdr", mixcompdist = ash.mixcompdist, mode = 0)
    },
    error = function(e) {
      print(e)
      res <- list(
        result = list(
          PosteriorMean = rep(0, length(ids)),
          PosteriorSD = rep(Inf, length(ids)),
          lfsr = rep(1, length(ids))
        )
      )
      return(res)
    }
  )

 contrast_res <- tidyr::tibble(
    id = ids, # row.names(PEM),
    raw_lfc = effect_est,
    raw_lfc_se = se_est,
    raw_p_value = pnorm(abs(effect_est), sd = se_est, lower.tail = FALSE),
    shrunken_lfc = shrunken_res$result$PosteriorMean,
    shrunken_lfc_se = shrunken_res$result$PosteriorSD,
    p_value = shrunken_res$result$lfsr,
    # ash_mixcompdist=ash.mixcompdist,
    effect_skew = up_down_large_effect_skew,
    log_mean_expression,
    coefficient_mode
  )

  if (!is.null(prefix)) {
    colnames(contrast_res)[2:7] <- paste(prefix, colnames(contrast_res)[2:7], sep = "_")
  }
  return(contrast_res)
}


#' Compare Genes in Cell State
#'
#' This function compares gene expression levels in a given cell state with its parent, sibling, and child states.
#'
#' @param cell_state A character string representing the cell state to be analyzed.
#' @param state_graph An igraph object representing the state graph.
#' @param estimate_matrix A matrix of estimated gene expression values.
#' @param stderr_matrix A matrix of standard errors for the estimated gene expression values.
#' @param state_term A character string representing the state term. Default is "cell_group".
#' @param log_fc_thresh A numeric value representing the log fold change threshold. Default is 1.
#' @param abs_expr_thresh A numeric value representing the absolute expression threshold. Default is 1e-3.
#' @param sig_thresh A numeric value representing the significance threshold. Default is 0.05.
#' @param cores An integer representing the number of cores to use for parallel processing. Default is 1.
#' @param expected_effect_mode_interval A numeric vector representing the expected effect mode interval. Default is c(-10, 10).
#' @param cv_threshold A numeric value representing the coefficient of variation threshold. Default is 10.
#'
#' @return A tibble containing gene expression comparisons and interpretations.
#'
#' @details
#' The function performs the following steps:
#' 1. Identifies parent, sibling, and child states of the given cell state.
#' 2. Computes gene expression levels and significance for the cell state and its related states.
#' 3. Compares gene expression levels between the cell state and its parent, sibling, and child states.
#' 4. Interprets the gene expression patterns based on the comparisons.
#'
#' @examples
#' \dontrun{
#' compare_genes_in_cell_state(
#'   cell_state = "state1",
#'   state_graph = my_state_graph,
#'   estimate_matrix = my_estimate_matrix,
#'   stderr_matrix = my_stderr_matrix
#' )
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import igraph
#' @import Matrix
#' @importFrom stats pnorm p.adjust
#' @export
compare_genes_in_cell_state <- function(cell_state,
                                        state_graph,
                                        estimate_matrix,
                                        stderr_matrix,
                                        n,
                                        state_term = "cell_group",
                                        log_fc_thresh = 1,
                                        abs_expr_thresh = 1e-3,
                                        sig_thresh = 0.05,
                                        cores = 1,
                                        expected_effect_mode_interval = c(-10, 10),
                                        cv_threshold = NULL) {
  parents <- get_parents(state_graph, cell_state) # igraph::neighbors(state_graph, cell_state, mode="in")
  parents <- intersect(parents, colnames(estimate_matrix))

  children <- get_children(state_graph, cell_state) # igraph::neighbors(state_graph, cell_state, mode="out")
  children <- intersect(children, colnames(estimate_matrix))

  siblings <- get_siblings(state_graph, cell_state) # igraph::neighbors(state_graph, parents, mode="out")
  siblings <- intersect(siblings, colnames(estimate_matrix))

  states_in_contrast <- c(cell_state, parents, children, siblings) %>% unique()

  message("      examining coefficients ", cell_state)

  expr_df <- tibble(gene_id = row.names(estimate_matrix))

  # if (is.null(ambient_estimate_matrix)) {

  # save the expr value
  expr_df$expr_self_est <- estimate_matrix[, cell_state]
  expr_df$expr_self_sd <- stderr_matrix[, cell_state]

  expr_df$expr_self <- pnorm(estimate_matrix[, cell_state] - log(abs_expr_thresh),
    sd = stderr_matrix[, cell_state], lower.tail = FALSE
  )

  expr_df$expr_self <- p.adjust(expr_df$expr_self, method = "BH") < sig_thresh
  genes_to_test <- expr_df$gene_id

  expr_df$expressed_in_parents <- NA
  expr_df$expressed_in_siblings <- NA
  expr_df$higher_than_parents <- NA
  expr_df$lower_than_parents <- NA
  expr_df$higher_than_all_siblings <- NA
  expr_df$lower_than_all_siblings <- NA
  expr_df$higher_than_siblings <- NA
  expr_df$lower_than_siblings <- NA
  expr_df$expressed_in_children <- NA
  expr_df$higher_than_all_children <- NA
  expr_df$lower_than_all_children <- NA
  expr_df$higher_than_children <- NA
  expr_df$lower_than_children <- NA


  if (length(parents) > 0) {
    # if (is.null(ambient_estimate_matrix)) {
    expressed_in_parents_mat <- pnorm(estimate_matrix[, parents, drop = F] - log(abs_expr_thresh),
      sd = stderr_matrix[, parents, drop = F], lower.tail = FALSE
    )
    expressed_in_parents_mat <- apply(expressed_in_parents_mat, 2, p.adjust, method = "BH")
    expressed_in_parents_mat <- expressed_in_parents_mat < sig_thresh
    expr_df$expressed_in_parents <- Matrix::rowSums(expressed_in_parents_mat) > 0

    cell_state_to_parents <- contrast_helper(cell_state, parents,
      PEM = estimate_matrix[genes_to_test, ],
      PSEM = stderr_matrix[genes_to_test, ],
      prefix = "cell_state_to_parents",
      ash.control = list(mode = expected_effect_mode_interval),
      cv_threshold = cv_threshold, 
      n = n, 
      abs_expr_thresh = abs_expr_thresh
    )


    parents_to_cell_state <- contrast_helper(parents, cell_state,
      PEM = estimate_matrix[genes_to_test, ],
      PSEM = stderr_matrix[genes_to_test, ],
      prefix = "parents_to_cell_state",
      ash.control = list(mode = expected_effect_mode_interval),
      cv_threshold = cv_threshold, 
      n = n, 
      abs_expr_thresh = abs_expr_thresh
    )

    expr_df <- left_join(expr_df, cell_state_to_parents, by = c("gene_id" = "id"))
    expr_df <- left_join(
      expr_df,
      # FIXME this is to deal with the case where there is no ambient expression provided, columns are duplicated
      parents_to_cell_state %>% select(-c(effect_skew, log_mean_expression, coefficient_mode)),
      by = c("gene_id" = "id")
    )

    expr_df$higher_than_parents <- p.adjust(expr_df$cell_state_to_parents_p_value, method = "BH") < sig_thresh &
      expr_df$cell_state_to_parents_shrunken_lfc > log_fc_thresh

    expr_df$lower_than_parents <- p.adjust(expr_df$parents_to_cell_state_p_value, method = "BH") < sig_thresh &
      expr_df$parents_to_cell_state_shrunken_lfc > log_fc_thresh
  } else {
    expr_df$expressed_in_parents <- NA
    expr_df$expressed_in_siblings <- NA
    expr_df$higher_than_parents <- NA
    expr_df$lower_than_parents <- NA
    expr_df$higher_than_all_siblings <- NA
    expr_df$lower_than_all_siblings <- NA
    expr_df$higher_than_siblings <- NA
    expr_df$lower_than_siblings <- NA
    expr_df$expressed_in_children <- NA
    expr_df$higher_than_all_children <- NA
    expr_df$lower_than_all_children <- NA
    expr_df$higher_than_children <- NA
    expr_df$lower_than_children <- NA
  }

  if (length(siblings) > 0) {
    # if (is.null(ambient_estimate_matrix)) {
    expressed_in_siblings_mat <- pnorm(estimate_matrix[, siblings, drop = F] - log(abs_expr_thresh),
      sd = stderr_matrix[, siblings, drop = F], lower.tail = FALSE
    )
    expressed_in_siblings_mat <- apply(expressed_in_siblings_mat, 2, p.adjust, method = "BH")
    expressed_in_siblings_mat <- expressed_in_siblings_mat < sig_thresh
    expr_df$expressed_in_siblings <- Matrix::rowSums(expressed_in_siblings_mat) > 0

    higher_than_siblings_mat <- lapply(siblings, function(sibling) {
      cell_state_to_sibling <- contrast_helper(cell_state, sibling,
        PEM = estimate_matrix[genes_to_test, ],
        PSEM = stderr_matrix[genes_to_test, ],
        prefix = "cell_state_to_sibling",
        ash.control = list(mode = expected_effect_mode_interval),
        cv_threshold = cv_threshold, 
        n = n, 
        abs_expr_thresh = abs_expr_thresh
      )

      p.adjust(cell_state_to_sibling$cell_state_to_sibling_p_value) < sig_thresh &
        cell_state_to_sibling$cell_state_to_sibling_shrunken_lfc > log_fc_thresh
    })
    higher_than_siblings_mat <- do.call(cbind, higher_than_siblings_mat)
    names(higher_than_siblings_mat) <- siblings
    row.names(higher_than_siblings_mat) <- genes_to_test

    lower_than_siblings_mat <- lapply(siblings, function(sibling) {
      sibling_to_cell_state <- contrast_helper(sibling, cell_state,
        PEM = estimate_matrix[genes_to_test, ],
        PSEM = stderr_matrix[genes_to_test, ],
        prefix = "sibling_to_cell_state",
        ash.control = list(mode = expected_effect_mode_interval),
        cv_threshold = cv_threshold, 
        n = n, 
        abs_expr_thresh = abs_expr_thresh
      )

      p.adjust(sibling_to_cell_state$sibling_to_cell_state_p_value) < sig_thresh &
        sibling_to_cell_state$sibling_to_cell_state_shrunken_lfc > log_fc_thresh
    })
    lower_than_siblings_mat <- do.call(cbind, lower_than_siblings_mat)
    names(lower_than_siblings_mat) <- siblings
    row.names(lower_than_siblings_mat) <- genes_to_test

    expr_df$higher_than_siblings <- left_join(expr_df,
      data.frame("higher_than_siblings" = Matrix::rowSums(higher_than_siblings_mat) > 0) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(higher_than_siblings.y)
    expr_df$lower_than_siblings <- left_join(expr_df,
      data.frame("lower_than_siblings" = Matrix::rowSums(lower_than_siblings_mat) > 0) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(lower_than_siblings.y)


    # change this to pairwise
    expr_df$higher_than_all_siblings <- left_join(expr_df,
      data.frame("higher_than_all_siblings" = Matrix::rowSums(higher_than_siblings_mat) == ncol(higher_than_siblings_mat)) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(higher_than_all_siblings.y)

    expr_df$lower_than_all_siblings <- left_join(expr_df,
      data.frame("lower_than_all_siblings" = Matrix::rowSums(lower_than_siblings_mat) == ncol(lower_than_siblings_mat)) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(lower_than_all_siblings.y)
  } else {
    expr_df$expressed_in_siblings <- NA
    expr_df$higher_than_all_siblings <- NA
    expr_df$lower_than_all_siblings <- NA
    expr_df$higher_than_siblings <- NA
    expr_df$lower_than_siblings <- NA
  }

  if (length(children) > 0) {
    # if (is.null(ambient_estimate_matrix)) {
    expressed_in_children_mat <- pnorm(estimate_matrix[, children, drop = F] - log(abs_expr_thresh),
      sd = stderr_matrix[, children, drop = F], lower.tail = FALSE
    )
    expressed_in_children_mat <- apply(expressed_in_children_mat, 2, p.adjust, method = "BH")
    expressed_in_children_mat <- expressed_in_children_mat < sig_thresh
    expr_df$expressed_in_children <- Matrix::rowSums(expressed_in_children_mat) > 0

    higher_than_children_mat <- lapply(children, function(child) {
      cell_state_to_child <- contrast_helper(cell_state, child,
        PEM = estimate_matrix[genes_to_test, ],
        PSEM = stderr_matrix[genes_to_test, ],
        prefix = "cell_state_to_child",
        ash.control = list(mode = expected_effect_mode_interval),
        cv_threshold = cv_threshold, 
        n = n, 
        abs_expr_thresh = abs_expr_thresh
      )

      p.adjust(cell_state_to_child$cell_state_to_child_p_value) < sig_thresh &
        cell_state_to_child$cell_state_to_child_shrunken_lfc > log_fc_thresh
    })
    higher_than_children_mat <- do.call(cbind, higher_than_children_mat)
    names(higher_than_children_mat) <- children
    row.names(higher_than_children_mat) <- genes_to_test

    lower_than_children_mat <- lapply(children, function(child) {
      child_to_cell_state <- contrast_helper(child, cell_state,
        PEM = estimate_matrix[genes_to_test, ],
        PSEM = stderr_matrix[genes_to_test, ],
        prefix = "child_to_cell_state",
        ash.control = list(mode = expected_effect_mode_interval),
        cv_threshold = cv_threshold, 
        n = n, 
        abs_expr_thresh = abs_expr_thresh
      )

      p.adjust(child_to_cell_state$child_to_cell_state_p_value) < sig_thresh &
        child_to_cell_state$child_to_cell_state_shrunken_lfc > log_fc_thresh
    })
    lower_than_children_mat <- do.call(cbind, lower_than_children_mat)
    names(lower_than_children_mat) <- children
    row.names(lower_than_children_mat) <- genes_to_test

    expr_df$higher_than_children <- left_join(expr_df,
      data.frame("higher_than_children" = Matrix::rowSums(higher_than_children_mat) > 0) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(higher_than_children.y)
    expr_df$lower_than_children <- left_join(expr_df,
      data.frame("lower_than_children" = Matrix::rowSums(lower_than_children_mat) > 0) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(lower_than_children.y)


    expr_df$higher_than_all_children <- left_join(expr_df,
      data.frame("higher_than_all_children" = Matrix::rowSums(higher_than_children_mat) == ncol(higher_than_children_mat)) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(higher_than_all_children.y)

    expr_df$lower_than_all_children <- left_join(expr_df,
      data.frame("lower_than_all_children" = Matrix::rowSums(lower_than_children_mat) == ncol(lower_than_children_mat)) %>% tibble::rownames_to_column("gene_id"),
      by = "gene_id"
    ) %>%
      pull(lower_than_all_children.y)
  } else {
    expr_df$expressed_in_children <- NA
    expr_df$higher_than_all_children <- NA
    expr_df$lower_than_all_children <- NA
    expr_df$higher_than_children <- NA
    expr_df$lower_than_children <- NA
  }

  expr_df <- expr_df %>% tidyr::nest(data = !gene_id)

  message("      interpreting patterns")
  interpret_expression_pattern <- function(pat_df) {
    if (pat_df$expr_self) {
      if (is.na(pat_df$expressed_in_parents)) {
        # no parents, therefore no siblings
        # return ("Maintained")
        if (is.na(pat_df$expressed_in_children)) {
          return("Maintained")
        } else {
          # no parent, but there are children
          if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children) {
            # Higher than parent, and higher than children
            return("Precursor-specific")
          } else if (pat_df$higher_than_children) {
            # no parent, higher than children
            return("Precursor-specific")
          } else if (pat_df$lower_than_children) {
            # no parent, but lower than children
            return("Precursor-depleted")
          } else { # no parent same as children
            return("Maintained")
          }
        }
      } else if (pat_df$expressed_in_parents) {
        # Expressed in self and parent
        if (is.na(pat_df$expressed_in_siblings)) {
          # Expressed in self and parent and there are no siblings
          if (pat_df$higher_than_parents) {
            if (is.na(pat_df$expressed_in_children)) {
              return("Upregulated")
            } else {
              # there are children
              if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children) {
                # Higher than parent, and higher than siblings
                return("Transiently upregulated")
              } else if (pat_df$lower_than_all_children) {
                # lower than children
                return("Increasingly upregulated")
              } else { # same as children
                return("Upregulated")
              }
            }
          } else if (pat_df$lower_than_parents) {
            if (is.na(pat_df$expressed_in_children)) {
              return("Downregulated")
            } else {
              # there are children
              if (pat_df$lower_than_all_children) {
                # Lower than parent, and lower than children
                return("Decreasingly downregulated")
              } else if (pat_df$lower_than_all_children) {
                # lower than children
                return("Transiently downregulated")
              } else { # same as children
                return("Downregulated")
              }
            }
          } else {
            if (is.na(pat_df$expressed_in_children)) {
              return("Maintained")
            } else {
              # same as parent, and there are children
              if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children) {
                # Higher than parent, and higher than children
                return("Precursor-specific")
              } else if (pat_df$lower_than_all_children) {
                # no parent, but lower than children
                return("Precursor-depleted")
              } else { # no parent same as children
                return("Maintained")
              }
            }
          }
        } else {
          # Expressed in self and parent and there are siblings
          if (pat_df$higher_than_parents) {
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings) {
              # Higher than parent, and higher than siblings
              return("Specifically upregulated")
            } else if (pat_df$higher_than_siblings) {
              # Higher than parent, and higher than siblings
              return("Selectively upregulated")
            } else if (pat_df$lower_than_siblings) {
              # Higher than parent, but lower than siblings
              return("Upregulated")
            } else { # higher than parent, same as siblings
              return("Upregulated")
            }
          } else if (pat_df$lower_than_parents) {
            if (pat_df$expressed_in_siblings == TRUE & pat_df$lower_than_all_siblings) {
              # Lower than parent, and higher than siblings
              return("Specifically downregulated")
            } else if (pat_df$expressed_in_siblings == TRUE & pat_df$lower_than_siblings) {
              # Lower than parent, and higher than some siblings
              return("Selectively downregulated")
            } else { # lower than parent, same as or higher than siblings
              return("Downregulated")
            }
          } else { # same as parent
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings) {
              # Same as parent, and higher than all siblings
              return("Specifically maintained")
            } else if (pat_df$higher_than_siblings) {
              # Same as parent, and higher than some siblings
              return("Selectively maintained")
            } else if (pat_df$lower_than_all_siblings) {
              # Same as parent, but lower than siblings
              return("Maintained")
            } else { # same as parent, same as siblings
              return("Maintained")
            }
          }
        }
      } else {
        # expressed in self but not in parent
        if (is.na(pat_df$expressed_in_siblings)) {
          # Expressed in self, not in parent and there are no siblings
          if (pat_df$higher_than_parents) {
            return("Activated")
          } else if (pat_df$lower_than_parents) {
            return("Downregulated")
          } # shouldn't happen
          else {
            return("Activated")
          } # might happen if its above threshold but not significantly above parent (and parent is below thresh)
        } else {
          # expressed in self, not in parent, and there are siblings
          if (pat_df$higher_than_parents) { # Expressed in self and higher than parent and there are siblings
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings) {
              # Higher than parent, and higher than all siblings
              return("Specifically activated")
            } else if (pat_df$higher_than_siblings) {
              # Higher than parent, and higher than some siblings
              return("Selectively activated")
            } else if (pat_df$lower_than_all_siblings) {
              # Higher than parent (which is off), but lower than all siblings
              return("Activated")
            }
            if (pat_df$lower_than_siblings) {
              # Higher than parent (which is off), but lower than some siblings
              return("Activated")
            } else { # Higher than parent (which is off), same as siblings
              return("Activated")
            }
          } else if (pat_df$lower_than_parents) {
            # gene is expressed, lower in the parent (which is off)
            if (pat_df$higher_than_all_siblings) {
              # Lower than parent, and higher than all siblings
              return("Absent") # shouldn't happen
            } else if (pat_df$higher_than_siblings) {
              # Lower than parent, and higher than some siblings
              return("Absent") # shouldn't happen
            } else if (pat_df$lower_than_all_siblings) {
              # Lower than parent and  lower than all siblings
              return("Absent")
            } else if (pat_df$lower_than_siblings) {
              # Lower than parent and  lower than some siblings
              return("Absent")
            } else { # Lower than parent, same as siblings
              return("Absent")
            }
          } else { # same as parent (which is off)
            if (pat_df$higher_than_all_siblings) {
              # Same as parent, and higher than all siblings
              return("Absent")
            } else if (pat_df$higher_than_siblings) {
              # Same as parent, and higher than some siblings
              return("Absent")
            } else if (pat_df$lower_than_all_siblings) {
              # Same as parent, but lower than all siblings
              return("Absent")
            } else if (pat_df$lower_than_siblings) {
              # Same as parent, but lower than some siblings
              return("Absent")
            } else { # same as parent, same as siblings
              return("Absent")
            }
          }
        }
      }
      return("Expressed")
    } else {
      # Not expressed in self
      if (is.na(pat_df$expressed_in_parents)) {
        # no parents, therefore no siblings
        return("Absent")
      } else if (pat_df$expressed_in_parents) {
        # Not expressed in self, but expressed in parents
        if (is.na(pat_df$expressed_in_siblings)) {
          # Not expressed in self, expressed parent and there are no siblings
          if (pat_df$lower_than_parents) {
            return("Deactivated")
          } else {
            return("Absent")
          } # shouldn't happen
        } else {
          # Not expressed in self, expressed in parent and there are siblings
          if (pat_df$lower_than_parents) {
            # Lower than parent
            if (pat_df$lower_than_all_siblings) {
              # Lower than parent and  lower than siblings
              return("Specifically deactivated")
            } else if (pat_df$lower_than_siblings) {
              # Lower than parent and  lower than siblings
              return("Selectively deactivated")
            }
            return("Deactivated")
          } else {
            # Not expressed in self, not lower than parent
            return("Absent")
          }
        }
      } else {
        # Not expressed in self or parents
        return("Absent")
      }
      return("Absent")
    }
    return("Absent")
    # match_row = match(data.frame(t(pat_df)), data.frame(t(interp_table)))
    # interpetation[match_row]
  }
  # debug(interpret_expression_pattern)
  #
  expr_df <- expr_df %>% mutate(interpretation = purrr::map(.f = purrr::possibly(
    interpret_expression_pattern, NA_character_
  ), .x = data))
  expr_df <- expr_df %>% tidyr::unnest(interpretation)
  message("      completed ", cell_state)
  return(expr_df)
}



#' Calculate Differentially Expressed Genes (DVEGs)
#'
#' This function calculates differentially expressed genes (DVEGs) by comparing
#' perturbation data with reference data. It identifies genes that are
#' significantly underexpressed or overexpressed based on a specified p-value threshold.
#'
#' @param perturb_degs A data frame containing perturbation differentially expressed genes (DEGs).
#' @param ref_degs A data frame containing reference differentially expressed genes (DEGs).
#' @param sig_p_val_thresh A numeric value specifying the significance p-value threshold for filtering DVEGs. Default is 1.
#'
#' @return A data frame containing significantly differentially expressed genes (DVEGs) with their dysregulation type.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item Mutates the `ref_degs` data frame to simplify the interpretation of gene expression.
#'   \item Joins the `ref_degs` and `perturb_degs` data frames based on gene ID and cell state.
#'   \item Selects and renames relevant columns for further analysis.
#'   \item Creates a display name for each gene based on its interpretation and cell state.
#'   \item Determines the dysregulation type (underexpressed or overexpressed) based on log fold change (LFC) and interpretation.
#'   \item Filters the results to include only significantly dysregulated genes based on the specified p-value threshold.
#' }
#'
#' @examples
#' \dontrun{
#' perturb_degs <- data.frame(id = c("gene1", "gene2"), cell_group = c("state1", "state2"), perturb_to_ctrl_shrunken_lfc = c(-1.5, 2.0), perturb_to_ctrl_p_value = c(0.01, 0.05))
#' ref_degs <- data.frame(gene_id = c("gene1", "gene2"), cell_state = c("state1", "state2"), interpretation = c("Upregulated", "Downregulated"), gene_short_name = c("G1", "G2"))
#' sig_dvegs <- calculate_dvegs(perturb_degs, ref_degs, sig_p_val_thresh = 0.05)
#' print(sig_dvegs)
#' }
#'
#' @import dplyr
#' @import stringr
#' @export
calculate_dvegs <- function(perturb_degs,
                            ref_degs,
                            sig_p_val_thresh = 1) {
  if ("gene_class_scores" %in% colnames(ref_degs)) {
    ref_degs <- ref_degs %>% tidyr::unnest(gene_class_scores)
  }

  if ("genes_within_cell_group" %in% colnames(perturb_degs)) {
    perturb_degs <- perturb_degs %>% tidyr::unnest(genes_within_cell_group)
  }

  if ("perturb_effects" %in% colnames(perturb_degs)) {
    perturb_degs <- perturb_degs %>% tidyr::unnest(perturb_effects)
  }


  ref_degs <- ref_degs %>%
    mutate(interpretation_simple = case_when(
      NA ~ NA,
      interpretation %in% c(
        "Upregulated",
        "Activated",
        "Selectively upregulated",
        "Specifically upregulated",
        "Increasingly upregulated",
        "Transiently upregulated",
        "Precursor-depleted"
      ) ~ "Up",
      interpretation %in% c(
        "Downregulated",
        "Deactivated",
        "Selectively downregulated",
        "Specifically downregulated",
        "Decreasingly downregulated",
        "Transiently downregulated",
        "Precursor-specific"
      ) ~ "Down",
      interpretation %in% c(
        "Maintained",
        "Specifically maintained",
        "Selectively maintained"
      ) ~ "Maintained",
      TRUE ~ interpretation
    ))

  wt_vs_perturb_degs <- ref_degs %>%
    select(cell_state, gene_id, gene_short_name, interpretation, interpretation_simple) %>%
    left_join(perturb_degs, by = c("gene_id" = "id", "cell_state" = "cell_group", "gene_short_name"))

  wt_vs_perturb_degs <- wt_vs_perturb_degs %>%
    select(
      cell_state,
      term,
      gene_id,
      gene_short_name,
      interpretation,
      interpretation_simple,
      perturb_to_ctrl_shrunken_lfc,
      perturb_to_ctrl_p_value
    ) %>%
    distinct()

  wt_vs_perturb_degs <- wt_vs_perturb_degs %>%
    mutate(display_name = stringr::str_c(gene_short_name, ":\n", "Normally",
      stringr::str_to_lower(interpretation), "in", cell_state,
      sep = " "
    ))

  wt_vs_perturb_degs <- wt_vs_perturb_degs %>%
    mutate(dysreg_type = case_when(
      interpretation_simple %in% c("Up", "Maintained") & perturb_to_ctrl_shrunken_lfc < 0 ~ "Underexpressed",
      interpretation_simple %in% c("Down", "Maintained") & perturb_to_ctrl_shrunken_lfc > 0 ~ "Overexpressed",
      TRUE ~ "Other DEG"
    ))

  sig_dvegs <- wt_vs_perturb_degs %>% filter(perturb_to_ctrl_p_value < sig_p_val_thresh & dysreg_type %in% c("Underexpressed", "Overexpressed"))
  sig_dvegs$dysreg_type <- factor(sig_dvegs$dysreg_type, levels = c("Underexpressed", "Overexpressed"))

  return(sig_dvegs)
}




#' Calculate Overrepresented Genes in Pathways
#'
#' This function identifies overrepresented genes in pathways for a given cell state
#' based on differential expression analysis. It uses the `fgsea` package to perform
#' overrepresentation analysis (ORA) and collapses pathways to main pathways for clarity.
#'
#' @param degs_for_cell_state A data frame containing differential expression results
#'   for a specific cell state. It must include columns `interpretation_simple`
#'   (with values "Up", "Down", or "Maintained") and `gene_short_name`.
#' @param pathway_type A character string specifying the type of pathway analysis to perform.
#'   (Currently unused in the function but included for potential future use.)
#' @param pathway_df A data frame containing pathway information. It must include columns
#'   `gene_short_name` (gene names) and `gs_name` (pathway names).
#' @param gene_universe A character vector of all possible genes to consider in the analysis.
#'   If `NULL`, the universe is derived from the input data (`degs_for_cell_state`).
#' @param q_val_thresh A numeric value specifying the q-value threshold for significance.
#'   Default is 0.05.
#'
#' @return A data frame containing the results of the overrepresentation analysis (ORA),
#'   including the pathways and their significance. The results are filtered to include
#'   only main pathways for "Up" and "Down" genes.
#'
#' @details
#' The function splits the input pathways into gene sets and performs ORA separately
#' for "Up" and "Down" genes. It uses the `fgsea::fora` function for ORA and
#' `fgsea::collapsePathwaysORA` to collapse pathways to main pathways. The results
#' are combined and returned as a single data frame.
#'
#' @import dplyr
#' @import fgsea
#' @export
calc_overrep_genes <- function(degs_for_cell_state,
                               pathway_type,
                               pathway_df,
                               gene_universe = NULL,
                               q_val_thresh = 0.05) {
  up_genes <- degs_for_cell_state %>%
    filter(interpretation_simple == "Up") %>%
    pull(gene_short_name) %>%
    unique()
  down_genes <- degs_for_cell_state %>%
    filter(interpretation_simple == "Down") %>%
    pull(gene_short_name) %>%
    unique()
  maint_genes <- degs_for_cell_state %>%
    filter(interpretation_simple == "Maintained") %>%
    pull(gene_short_name) %>%
    unique()
  if (is.null(gene_universe)) {
    gene_universe <- unique(c(up_genes, down_genes, maint_genes))
  } else {
    gene_universe <- unique(c(up_genes, down_genes, gene_universe))
  }


  gene_set_list <- split(x = pathway_df$gene_short_name, f = pathway_df$gs_name)
  up_fora_res <- NULL
  if (length(up_genes) > 1) {
    # print (head(gene_set_list))
    # print (head(up_genes))
    # print (head(gene_universe))
    up_fora_res <- fgsea::fora(gene_set_list, genes = up_genes, universe = gene_universe)
    up_fora_res$interpretation_simple <- "Up"
    up_fora_res <- up_fora_res[up_fora_res$padj < q_val_thresh,]
    if (nrow(up_fora_res) > 0) {
      up_fora_res_collapsed <- fgsea::collapsePathwaysORA(up_fora_res[order(up_fora_res$pval),],
        gene_set_list,
        genes = up_genes,
        universe = gene_universe
      )
      if (length(up_fora_res_collapsed$mainPathways) > 0) {
        up_fora_res <- up_fora_res[pathway %in% up_fora_res_collapsed$mainPathways,]
      }
    }
  }

  down_fora_res <- NULL
  if (length(down_genes) > 0) {
    # print (head(gene_set_list))
    # print (head(down_genes))
    # print (head(gene_universe))
    down_fora_res <- fgsea::fora(gene_set_list, genes = down_genes, universe = gene_universe)
    down_fora_res$interpretation_simple <- "Down"
    down_fora_res <- down_fora_res[down_fora_res$padj < q_val_thresh,]
    if (nrow(down_fora_res) > 0) {
      down_fora_res_collapsed <- fgsea::collapsePathwaysORA(down_fora_res[order(down_fora_res$pval),],
        gene_set_list,
        genes = down_genes,
        universe = gene_universe
      )
      if (length(down_fora_res_collapsed$mainPathways) > 0) {
        down_fora_res <- down_fora_res[pathway %in% down_fora_res_collapsed$mainPathways,]
      }
    }
  }

  fora_res <- bind_rows(up_fora_res, down_fora_res)
  fora_res <- fora_res %>% relocate(interpretation_simple)

  return(fora_res)
}


#' Calculate Graph Overrepresented Genes
#'
#' This function calculates overrepresented pathways for genes in a given dataset,
#' grouped by cell state. It processes the input data, simplifies gene interpretation,
#' and computes overrepresented pathways for each cell state.
#'
#' @param ref_degs A data frame containing reference differentially expressed genes (DEGs).
#'   It should include columns such as `cell_state`, `gene_short_name`, `gene_id`,
#'   and `interpretation`. If the data frame contains a `gene_class_scores` column,
#'   it will be unnested.
#' @param pathway_type A character vector specifying the type of pathways to analyze.
#'   Options are `"GO:BP"` (Biological Process) or `"GO:MF"` (Molecular Function).
#'   Defaults to `c("GO:BP", "GO:MF")`.
#' @param species A character string specifying the species for pathway analysis.
#'   Defaults to `"Danio rerio"`.
#'
#' @return A nested data frame where each row corresponds to a cell state, and the
#'   `overrep_pathways` column contains the overrepresented pathways for that cell state.
#'
#' @details The function simplifies the `interpretation` column into categories
#'   (`"Up"`, `"Down"`, `"Maintained"`, or retains the original value if none match).
#'   It then groups the data by `cell_state` and calculates overrepresented pathways
#'   using the `calc_overrep_genes` function.
#'
#' @examples
#' # Example usage:
#' ref_degs <- data.frame(
#'   cell_state = c("State1", "State2"),
#'   gene_short_name = c("GeneA", "GeneB"),
#'   gene_id = c("ID1", "ID2"),
#'   interpretation = c("Upregulated", "Downregulated")
#' )
#' result <- calc_graph_overrep_genes(ref_degs, pathway_type = "GO:BP", species = "Danio rerio")
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#'
#' @export
calc_graph_overrep_genes <- function(ref_degs,
                                     pathway_type = c("GO:BP", "GO:MF"),
                                     species = "Danio rerio") {
  pathway_type <- match.arg(pathway_type)
  
  gene_set =  msigdbr::msigdbr(species = species, subcategory = pathway_type)
  gene_set_df = gene_set %>% dplyr::distinct(gs_name, gene_short_name=gene_symbol) %>% as.data.frame()
  

  if ("gene_class_scores" %in% colnames(ref_degs)) {
    ref_degs <- ref_degs %>% tidyr::unnest(gene_class_scores)
  }

  deg_nested_df <- ref_degs %>%
    select(cell_state, gene_short_name, gene_id, interpretation) %>%
    mutate(interpretation_simple = case_when(
      NA ~ NA,
      interpretation %in% c(
        "Upregulated",
        "Activated",
        "Selectively upregulated",
        "Specifically upregulated",
        "Increasingly upregulated",
        "Transiently upregulated",
        "Precursor-depleted"
      ) ~ "Up",
      interpretation %in% c(
        "Downregulated",
        "Deactivated",
        "Selectively downregulated",
        "Specifically downregulated",
        "Decreasingly downregulated",
        "Transiently downregulated",
        "Precursor-specific"
      ) ~ "Down",
      interpretation %in% c(
        "Maintained",
        "Specifically maintained",
        "Selectively maintained"
      ) ~ "Maintained",
      TRUE ~ interpretation
    )) %>%
    group_by(cell_state) %>%
    tidyr::nest()

  graph_degs_overrep_pathways <- deg_nested_df %>%
    mutate(overrep_pathways = purrr::map(
      .f = calc_overrep_genes,
      .x = data,
      pathway_type = pathway_type,
      pathway_df = gene_set_df
    ))

  return(graph_degs_overrep_pathways)
}
