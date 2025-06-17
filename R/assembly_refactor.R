

wt_assembly <- function(cds,
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
                              links_between_components = c("ctp", "none", "strongest-pcor", "strong-pcor"),
                              component_col = "partition",
                              embryo_size_factors = NULL,
                              batches_excluded_from_assembly = c()) {
  
    colData(cds)$subassembly_group <- stringr::str_c(partition_name, colData(cds)[, cell_group], sep = "-")
    colData(cds)[["cell_state"]] <- as.character(colData(cds)[[cell_group]])
    
    selected_colData <- colData(cds) %>%
        tibble::as_tibble() %>%
        dplyr::select(cell, !!sym(sample_group), cluster, !!sym(cell_group), subassembly_group, cell_state)

    selected_colData$cds_row_id <- colData(cds) %>%
        as.data.frame() %>%
        row.names()

    partition_results <- selected_colData %>% 
      tidyr::nest(data = c(cds_row_id, cell, !!sym(sample_group), cluster, !!sym(cell_group), subassembly_group))
    
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
                num_bootstraps = num_bootstraps,
                embryo_size_factors = embryo_size_factors
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
                newdata = newdata,
                links_between_components = links_between_components,
                ctrl_ids = ctrl_ids,
                edge_allowlist = edge_allowlist,
                sparsity_factor = sparsity_factor,
                perturbation_col = perturbation_col,
                component_col = component_col,
                verbose = verbose
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
            } else {
                wt_graph <- list(NA)
            }
        },
        error = function(e) {
            print(e)
            wt_graph <- list(NA)
            return(wt_graph)
        }
    )

    return(wt_graph[[1]])
}



mt_assembly <- function(cds,
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
                              perturbation_col = "gene_target",
                              batch_col = "expt",
                              max_num_cells = NULL,
                              verbose = FALSE,
                              keep_ccs = TRUE,
                              num_threads = 1,
                              backend = "nlopt",
                              q_val = 0.1,
                              vhat_method = "bootstrap",
                              num_bootstraps = 10,
                              newdata = tibble(),
                              edge_allowlist = NULL,
                              min_lfc = 0,
                              links_between_components = c("ctp", "none", "strongest-pcor", "strong-pcor"),
                              log_abund_detection_thresh = -5,
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

    perturb_models_tbl <- assess_perturbation_effects(wt_ccm,
        perturb_models_tbl,
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


    message("Assembling mutant graphs...")
    mt_graph <- assemble_mt_graph(wt_ccm,
        perturb_models_tbl,
        newdata = newdata,
        start_time = start_time,
        stop_time = stop_time,
        interval_col = interval_col,
        links_between_components = links_between_components,
        # perturbation_col = perturbation_col,
        # edge_allowlist = edge_allowlist,
        component_col = component_col,
        verbose = verbose
    )

    return(mt_graph)
}
