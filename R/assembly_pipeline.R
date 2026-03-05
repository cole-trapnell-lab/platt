#' Assemble a wild-type state graph
#'
#' Fit the wild-type model and assemble a state graph for a partition of the data.
#'
#' @param cds A SingleCellExperiment/CellDataSet with expression and metadata.
#' @param sample_group Column name in colData specifying sample grouping.
#' @param cell_group Column name in colData specifying cell grouping.
#' @param partition_name Optional partition label appended to graph nodes.
#' @param main_model_formula_str Model formula string for the main effect.
#' @param start_time Start time for model fitting.
#' @param stop_time Stop time for model fitting.
#' @param interval_col Column name for time intervals.
#' @param nuisance_model_formula_str Nuisance model formula string.
#' @param ctrl_ids Optional vector of control IDs.
#' @param sparsity_factor Sparsity factor for model fitting.
#' @param perturbation_col Column name for perturbation labels.
#' @param batch_col Column name for batch labels.
#' @param verbose Logical, whether to log progress.
#' @param keep_ccs Logical, whether to keep cell_count_set outputs.
#' @param num_threads Number of threads to use.
#' @param backend Optimization backend to use.
#' @param q_val FDR threshold.
#' @param vhat_method Method to estimate vhat.
#' @param num_bootstraps Number of bootstraps for vhat.
#' @param newdata Optional data frame of new timepoints for prediction.
#' @param edge_allowlist Optional allowlist of edges.
#' @param edge_denylist Optional denylist of edges.
#' @param links_between_components Strategy for linking graph components.
#' @param component_col Column name for component labels.
#' @param embryo_size_factors Optional size factors for embryo data.
#' @param log_abund_detection_thresh Log abundance detection threshold.
#' @param interval_step Step size for interval grid.
#' @param min_interval Minimum interval length.
#' @param max_interval Maximum interval length.
#' @param min_pathfinding_lfc Minimum log-fold-change for pathfinding.
#' @param num_time_breaks Number of time breaks for fitting.
#' @param batches_excluded_from_assembly Vector of batches to exclude.
#' @param force_allowlist Logical, whether to force the edge allowlist.
#' @param min_penalty Minimum penalty for model fitting.
#' @param max_penalty Maximum penalty for model fitting.
#' @param break_cycles Logical, whether to break cycles in the graph.
#'
#' @return An igraph object for the assembled wild-type graph, or NA on failure.
#'
#' @export
run_wildtype_assembly <- function(cds,
                                  sample_group,
                                  cell_group,
                                  partition_name = NULL,
                                  main_model_formula_str = NULL,
                                  start_time = 18,
                                  stop_time = 72,
                                  interval_col = "timepoint",
                                  nuisance_model_formula_str = "~1",
                                  ctrl_ids = NULL,
                                  sparsity_factor = 0.01,
                                  perturbation_col = "perturbation",
                                  batch_col = "expt",
                                  verbose = FALSE,
                                  keep_ccs = TRUE,
                                  num_threads = 1,
                                  backend = "nlopt",
                                  q_val = 0.1,
                                  vhat_method = "bootstrap",
                                  num_bootstraps = 10,
                                  newdata = tibble(),
                                  edge_allowlist = NULL,
                                  edge_denylist = NULL,
                                  links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                                  component_col = "partition",
                                  embryo_size_factors = NULL,
                                  log_abund_detection_thresh = -5,
                                  interval_step = 2,
                                  min_interval = 4,
                                  max_interval = 24,
                                  min_pathfinding_lfc = 0,
                                  num_time_breaks = 4,
                                  batches_excluded_from_assembly = c(),
                                  force_allowlist = FALSE,
                                  min_penalty = 0.01,
                                  max_penalty = 1e+06,
                                  break_cycles = TRUE) {
    colData(cds)$subassembly_group <- stringr::str_c(partition_name, colData(cds)[, cell_group], sep = "-")
    colData(cds)[["cell_state"]] <- as.character(colData(cds)[[cell_group]])

    selected_colData <- colData(cds) %>%
        tibble::as_tibble() %>%
        dplyr::select(cell, !!sym(sample_group), !!sym(cell_group), subassembly_group, cell_state)

    selected_colData$cds_row_id <- colData(cds) %>%
        as.data.frame() %>%
        row.names()

    partition_results <- selected_colData %>%
        tidyr::nest(data = c(cds_row_id, cell, !!sym(sample_group), !!sym(cell_group), subassembly_group))

    # if there is only one cell state, return NA
    if (length(unique(selected_colData[["cell_state"]])) <= 1) {
        wt_graph <- list(NA)
        return(wt_graph)
    }

    cds <- cds[, colData(cds)[[batch_col]] %in% batches_excluded_from_assembly == FALSE]

    tryCatch(
        {
            message("Starting wild-type fit...")
            wt_ccm <- suppressWarnings(fit_wt_model(cds,
                sample_group = sample_group,
                cell_group = cell_group,
                main_model_formula_str = main_model_formula_str,
                start_time = start_time,
                stop_time = stop_time,
                interval_col = interval_col,
                nuisance_model_formula_str = nuisance_model_formula_str,
                ctrl_ids = ctrl_ids,
                sparsity_factor = sparsity_factor,
                perturbation_col = perturbation_col,
                batch_col = batch_col,
                keep_ccs = keep_ccs,
                verbose = verbose,
                num_threads = num_threads,
                backend = backend,
                vhat_method = vhat_method,
                edge_allowlist = edge_allowlist,
                edge_denylist = edge_denylist,
                num_bootstraps = num_bootstraps,
                embryo_size_factors = embryo_size_factors,
                num_time_breaks = num_time_breaks,
                min_penalty = min_penalty,
                max_penalty = max_penalty
            ))

            if (is.null(wt_ccm) || is.na(wt_ccm)) {
                partition_results$wt_graph <- list(NA)
                partition_results$mt_graph <- list(NA)
                partition_results$perturbation_effects <- list(NA)
                partition_results$wt_state_graph_plot <- list(NA)
                partition_results$mt_state_graph_plot <- list(NA)
                return(partition_results)
            }

            message("Assembling wild-type graph...")
            wt_graph <- assemble_wt_graph(cds,
                wt_ccm,
                sample_group = sample_group,
                cell_group = cell_group,
                main_model_formula_str = main_model_formula_str,
                start_time = start_time,
                stop_time = stop_time,
                interval_col = interval_col,
                interval_step = interval_step,
                min_interval = min_interval,
                max_interval = max_interval,
                min_pathfinding_lfc = min_pathfinding_lfc,
                newdata = newdata,
                log_abund_detection_thresh = log_abund_detection_thresh,
                links_between_components = links_between_components,
                ctrl_ids = ctrl_ids,
                edge_allowlist = edge_allowlist,
                edge_denylist = edge_denylist,
                sparsity_factor = sparsity_factor,
                perturbation_col = perturbation_col,
                component_col = component_col,
                verbose = verbose,
                q_val = q_val,
                force_allowlist = force_allowlist,
                break_cycles = break_cycles
            )

            if (is.null(wt_graph) == FALSE) {
                if (cell_group == "cell_state") {
                    igraph::V(wt_graph)$name <- stringr::str_c(partition_name, igraph::V(wt_graph)$name, sep = "-")
                } else {
                    igraph::V(wt_graph)$name <- igraph::V(wt_graph)$name
                }
                if (is.null(partition_name)) {
                    partition_name <- ""
                }
                igraph::E(wt_graph)$assembly_group <- partition_name
                wt_graph <- list(wt_graph)
                return(wt_graph[[1]])
            } else {
                wt_graph <- list(NA)
                return(wt_graph)
            }
        },
        error = function(e) {
            print(e)
            wt_graph <- list(NA)
            return(wt_graph)
        }
    )
}


#' Assemble a mutant state graph
#'
#' Fit mutant models and assemble a state graph based on perturbation effects.
#'
#' @param cds A SingleCellExperiment/CellDataSet with expression and metadata.
#' @param sample_group Column name in colData specifying sample grouping.
#' @param cell_group Column name in colData specifying cell grouping.
#' @param wt_graph Wild-type graph to anchor mutant assembly.
#' @param partition_name Optional partition label.
#' @param main_model_formula_str Model formula string for the main effect.
#' @param start_time Start time for model fitting.
#' @param stop_time Stop time for model fitting.
#' @param interval_col Column name for time intervals.
#' @param nuisance_model_formula_str Nuisance model formula string.
#' @param ctrl_ids Optional vector of control IDs.
#' @param mt_ids Optional vector of mutant IDs.
#' @param sparsity_factor Sparsity factor for model fitting.
#' @param perturbation_col Column name for perturbation labels.
#' @param batch_col Column name for batch labels.
#' @param max_num_cells Optional cap on number of cells.
#' @param verbose Logical, whether to log progress.
#' @param keep_ccs Logical, whether to keep cell_count_set outputs.
#' @param num_threads Number of threads to use.
#' @param backend Optimization backend to use.
#' @param q_val FDR threshold.
#' @param interval_step Step size for interval grid.
#' @param vhat_method Method to estimate vhat.
#' @param num_bootstraps Number of bootstraps for vhat.
#' @param newdata Optional data frame of new timepoints for prediction.
#' @param edge_allowlist Optional allowlist of edges.
#' @param min_lfc Minimum log-fold-change for perturbation effects.
#' @param links_between_components Strategy for linking graph components.
#' @param log_abund_detection_thresh Log abundance detection threshold.
#' @param discordant_config Optional list controlling discordant pruning with
#'   keys `power_threshold`, `prune_mode`, `K_paths`, and `lambda_edge`.
#' @param batches_excluded_from_assembly Vector of batches to exclude.
#' @param component_col Column name for component labels.
#' @param embryo_size_factors Optional size factors for embryo data.
#'
#' @return An igraph object for the assembled mutant graph, or the WT graph on failure.
run_perturbation_assembly <- function(cds,
                                      sample_group,
                                      cell_group,
                                      wt_graph,
                                      partition_name = NULL,
                                      main_model_formula_str = NULL,
                                      start_time = 18,
                                      stop_time = 72,
                                      interval_col = "timepoint",
                                      nuisance_model_formula_str = "~1",
                                      ctrl_ids = NULL,
                                      mt_ids = NULL,
                                      sparsity_factor = 0.01,
                                      perturbation_col = "perturbation",
                                      batch_col = "expt",
                                      max_num_cells = NULL,
                                      verbose = FALSE,
                                      keep_ccs = TRUE,
                                      num_threads = 1,
                                      backend = "nlopt",
                                      q_val = 0.1,
                                      interval_step = 2,
                                      vhat_method = "bootstrap",
                                      num_bootstraps = 10,
                                      newdata = tibble(),
                                      edge_allowlist = NULL,
                                      min_lfc = 0,
                                      links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                                      log_abund_detection_thresh = -5,
                                      discordant_config = NULL,
                                      discordant_pruning_mode = c("none", "greedy"),
                                      discordant_pruning_k = 1,
                                      discordant_pruning_power_threshold = 0,
                                      discordant_pruning_cost_attr = "total_path_score_supporting",
                                      discordant_pruning_lambda_edge = 0,
                                      batches_excluded_from_assembly = c(),
                                      component_col = "partition",
                                      embryo_size_factors = NULL) {
    message("Starting mutant fits...")
    perturb_models_tbl <- suppressWarnings(fit_mt_models(cds,
        sample_group = sample_group,
        cell_group = cell_group,
        main_model_formula_str = main_model_formula_str,
        nuisance_model_formula_str = nuisance_model_formula_str,
        start_time = start_time,
        stop_time = stop_time,
        interval_col = interval_col,
        ctrl_ids = ctrl_ids,
        mt_ids = mt_ids,
        sparsity_factor = sparsity_factor,
        perturbation_col = perturbation_col,
        verbose = verbose,
        num_threads = num_threads,
        backend = backend,
        keep_ccs = keep_ccs,
        batch_col = batch_col,
        # edge_allowlist = edge_allowlist,
        vhat_method = vhat_method,
        num_bootstraps = num_bootstraps,
        embryo_size_factors = embryo_size_factors
    ))

    perturb_models_tbl <- perturb_models_tbl %>% filter(!is.na(perturb_ccm))

    if (is.null(perturb_models_tbl)) {
        partition_results$wt_graph <- list(NA)
        partition_results$mt_graph <- list(NA)
        partition_results$perturbation_effects <- list(NA)
        partition_results$wt_state_graph_plot <- list(NA)
        partition_results$mt_state_graph_plot <- list(NA)
        return(partition_results)
        # stop("Error: fit_mt_models() failed")
    }

    perturb_models_tbl <- assess_perturbation_effects(perturb_models_tbl,
        q_val = q_val,
        start_time = start_time,
        stop_time = stop_time,
        # perturbation_col = perturbation_col,
        interval_col = interval_col,
        log_abund_detection_thresh = log_abund_detection_thresh,
        min_lfc = min_lfc,
        verbose = verbose,
        newdata = newdata
    )

    # this makes a prediction for every measured timepoint
    # this is for useful to save for plotting later
    perturb_models_tbl <- perturb_models_tbl %>%
        mutate(perturbation_table = purrr::map(
            .f = purrr::possibly(get_perturbation_effects),
            .x = perturb_ccm,
            interval_col = interval_col,
            newdata = newdata
        ))

    ccs <- new_cell_count_set(cds, sample_group = sample_group, cell_group = cell_group)
    message("Assembling mutant graphs...")
    mt_graph <- assemble_mt_graph(ccs,
        wt_graph,
        perturb_models_tbl,
        newdata = newdata,
        start_time = start_time,
        stop_time = stop_time,
        interval_col = interval_col,
        interval_step = interval_step,
        q_val = q_val,
        log_abund_detection_thresh = log_abund_detection_thresh,
        links_between_components = links_between_components,
        discordant_config = discordant_config,
        discordant_pruning_mode = discordant_pruning_mode,
        discordant_pruning_k = discordant_pruning_k,
        discordant_pruning_power_threshold = discordant_pruning_power_threshold,
        discordant_pruning_cost_attr = discordant_pruning_cost_attr,
        discordant_pruning_lambda_edge = discordant_pruning_lambda_edge,
        component_col = component_col,
        verbose = verbose
    )

    if (is.na(mt_graph)) {
        print("returning wt_graph")
        return(wt_graph)
    }
    return(mt_graph)
}

# Backward-compatible wrappers
#' @export
wt_assembly <- function(...) {
    .Deprecated("run_wildtype_assembly")
    run_wildtype_assembly(...)
}

#' @export
mt_assembly <- function(...) {
    .Deprecated("run_perturbation_assembly")
    run_perturbation_assembly(...)
}
