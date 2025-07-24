#' @export
get_time_window <- function(genotype, ccs, interval_col, perturbation_col = "perturbation") {
  subset_ccs <- ccs[, replace_na(colData(ccs)[[perturbation_col]] %in% genotype, F)]
  colData(subset_ccs)$knockout <- colData(subset_ccs)[[perturbation_col]] %in% genotype
  knockout_time_start <- min(colData(subset_ccs)[[interval_col]][colData(subset_ccs)$knockout])
  knockout_time_stop <- max(colData(subset_ccs)[[interval_col]][colData(subset_ccs)$knockout])
  return(tibble(start_time = knockout_time_start, stop_time = knockout_time_stop))
}

#' @export
get_perturbation_effects <- function(ccm, interval_col = "timepoint", newdata = tibble(), adjust_q_values = FALSE) {
  timepoints <- colData(ccm@ccs)[[interval_col]] %>% unique()
  df <- data.frame(timepoint = timepoints)

  if (nrow(newdata) > 0) {
    df <- cross_join(df, newdata)
  }

  df <- df %>%
    group_split(row_number(), .keep = FALSE) %>%
    purrr::map_df(tidyr::nest) %>%
    mutate(genotype_eff = purrr::map(
      .f = make_contrast,
      .x = data,
      ccm = ccm,
      adjust_q_values = adjust_q_values
    )) %>%
    unnest(c(data, genotype_eff))
  return(df)
}

#' Fit Genotype Cell Count Model (CCM)
#'
#' This function fits a cell count model (CCM) for a given genotype using a
#' specified dataset of cell counts (CCS). It allows for flexible modeling of
#' time intervals, perturbations, and batch effects, and supports various
#' backend options for optimization.
#'
#' @param genotype Character. The genotype to fit the model for.
#' @param ccs A cell count dataset (CCS) object containing the data to model.
#' @param prior_state_transition_graph Optional. A prior state transition graph
#'   to use as a constraint in the model.
#' @param interval_col Character. The column name in `ccs` that represents the
#'   time intervals. Default is "timepoint".
#' @param perturbation_col Character. The column name in `ccs` that represents
#'   the perturbation (e.g., perturbation). Default is "perturbation".
#' @param batch_col Character. The column name in `ccs` that represents the
#'   batch information. Default is "expt".
#' @param ctrl_ids Character vector. A list of control IDs to include in the
#'   model. Default includes several predefined control IDs.
#' @param contrast_time_start Numeric. The start time for the contrast interval.
#'   If NULL, the minimum time in the data is used. Default is NULL.
#' @param contrast_time_stop Numeric. The stop time for the contrast interval.
#'   If NULL, the maximum time in the data is used. Default is NULL.
#' @param num_time_breaks Numeric. The number of time breakpoints to use for
#'   spline modeling. Default is NULL, which sets it to 3.
#' @param independent_spline_for_ko Logical. Whether to use an independent
#'   spline for knockout modeling. Default is TRUE.
#' @param edge_allowlist Optional. A list of allowed edges for the model. Default
#'   is NULL.
#' @param edge_denylist Optional. A list of denied edges for the model. Default
#'   is NULL.
#' @param penalize_by_distance Logical. Whether to penalize edges by distance.
#'   Default is TRUE.
#' @param keep_ccs Logical. Whether to keep the CCS object in the model. Default
#'   is TRUE.
#' @param vhat_method Character. The method to compute vhat. Default is
#'   "bootstrap".
#' @param num_threads Numeric. The number of threads to use for computation.
#'   Default is 1.
#' @param backend Character. The backend to use for optimization. Default is
#'   "nlopt".
#' @param sparsity_factor Numeric. The sparsity factor for model selection.
#'   Default is 0.01.
#' @param num_bootstraps Numeric. The number of bootstraps to use for vhat
#'   computation. Default is 10.
#' @param ftol_rel Numeric. The relative tolerance for optimization. Default is
#'   1e-06.
#'
#' @return A fitted genotype CCM object with additional metadata and model
#'   information.
#'
#' @details
#' The function subsets the input CCS data based on the specified genotype,
#' control IDs, and time intervals. It constructs a model formula based on the
#' number of time points and other parameters, and fits the model using the
#' specified backend. Batch effects are corrected if multiple batches are
#' present. The function also validates the rank of the model matrix to ensure
#' it is full rank.
#'
#' @examples
#' # Example usage:
#' genotype_ccm <- fit_genotype_ccm(
#'   genotype = "example_genotype",
#'   ccs = example_ccs,
#'   interval_col = "timepoint",
#'   perturbation_col = "gene_target",
#'   batch_col = "expt"
#' )
#'
#' @export
fit_genotype_ccm <- function(genotype,
                             ccs,
                             prior_state_transition_graph = NULL,
                             interval_col = "timepoint",
                             perturbation_col = "perturbation",
                             batch_col = NULL,
                             ctrl_ids = c("ctrl-uninj", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                             contrast_time_start = NULL,
                             contrast_time_stop = NULL,
                             num_time_breaks = NULL,
                             independent_spline_for_ko = TRUE,
                             edge_allowlist = NULL,
                             edge_denylist = NULL,
                             penalize_by_distance = TRUE,
                             keep_ccs = TRUE,
                             vhat_method = "bootstrap",
                             num_threads = 1,
                             backend = "nlopt",
                             sparsity_factor = 0.01,
                             num_bootstraps = 10,
                             ftol_rel = 1e-06) {
  message(paste("Fitting knockout model for", genotype))
  # subset_ccs = ccs[,colData(ccs)$gene_target == genotype | colData(ccs)$gene_target %in% ctrl_ids]

  if (is.null(batch_col)) {
    colData(ccs)[["batch"]] <- "DUMMY"
    batch_col <- "batch"
  }

  if (!is.null(ccs@cds@metadata$umap_space)) {
    print(paste0("You are running this in ", ccs@cds@metadata$umap_space))
  }

  subset_ccs <- ccs[, replace_na(colData(ccs)[[perturbation_col]] == genotype, F)]
  # expts = unique(colData(subset_ccs)[[batch_col]])

  if (is.null(contrast_time_start)) {
    knockout_time_start <- min(colData(subset_ccs)[[interval_col]])
  } else {
    knockout_time_start <- contrast_time_start
  }

  if (is.null(contrast_time_stop)) {
    knockout_time_stop <- max(colData(subset_ccs)[[interval_col]])
  } else {
    knockout_time_stop <- contrast_time_stop
  }
  subset_ccs <- subset_ccs[, replace_na(colData(subset_ccs)[[interval_col]] <= knockout_time_stop, F)]
  subset_ccs <- subset_ccs[, replace_na(colData(subset_ccs)[[interval_col]] >= knockout_time_start, F)]

  num_knockout_timepoints <- length(unique(colData(subset_ccs)[[interval_col]]))

  message(paste("\ttime range:", knockout_time_start, "to", knockout_time_stop))
  # subset_ccs = ccs[,( replace_na(colData(ccs)[[perturbation_col]] == genotype, F) | colData(ccs)[[perturbation_col]] %in% ctrl_ids) & colData(ccs)[[batch_col]] %in% expts]
  subset_ccs <- ccs[, (replace_na(colData(ccs)[[perturbation_col]] == genotype, F) | colData(ccs)[[perturbation_col]] %in% ctrl_ids)]
  expts <- unique(colData(subset_ccs)[[batch_col]])

  colData(subset_ccs)$knockout <- colData(subset_ccs)[[perturbation_col]] == genotype
  subset_ccs <- subset_ccs[, (colData(subset_ccs)[[interval_col]] >= knockout_time_start & colData(subset_ccs)[[interval_col]] <= knockout_time_stop)]
  time_breakpoints <- c()

  if (is.null(num_time_breaks)) {
    num_time_breaks <- 3
  }

  # Set up the knockout model formula. Considers the number of timepoints in the knockout
  if (num_knockout_timepoints > 2 & knockout_time_stop > knockout_time_start) {
    time_breakpoints <- seq(knockout_time_start, knockout_time_stop, length.out = num_time_breaks)
    time_breakpoints <- time_breakpoints[2:(length(time_breakpoints) - 1)] # exclude the first and last entry as these will become boundary knots

    time_term <- paste("~ ns(", interval_col, ", knots=", paste("c(", paste(time_breakpoints, collapse = ","), ")", sep = ""), ")")

    # If there are only two timepoints for a knockout, we can't have interaction terms between
    # the time spline and the knockout indicator variable. The model won't be
    # full rank
    if (independent_spline_for_ko & num_knockout_timepoints >= 3) {
      knockout_terms <- paste("~ ns(", interval_col, ", knots=", paste("c(", paste(time_breakpoints, collapse = ","), ")", sep = ""), "):knockout + knockout")
    } else {
      knockout_terms <- "~ knockout"
    }

    main_model_formula_str <- knockout_terms # paste(time_term, knockout_terms , sep="+")
    # print (main_model_formula_str)
    nuisance_model_formula_str <- time_term
  } else {
    time_term <- paste("~ ", interval_col)
    main_model_formula_str <- "~ knockout"
    if (num_knockout_timepoints > 1) {
      nuisance_model_formula_str <- time_term
    } else {
      nuisance_model_formula_str <- "~ 1"
    }
  }

  # make this any column
  if (length(unique(colData(subset_ccs)[[batch_col]])) > 1) {
    colData(subset_ccs)[[batch_col]] <- as.factor(colData(subset_ccs)[[batch_col]])
    # FIXME: This is identical code to what is immediately after the final construction of the model matrix, there may be value in writing a generic checker function
    # FIXME: This is also bespoke and may not be the best strategy to use in the general case.
    full_model_matrix <- Matrix::sparse.model.matrix(
      as.formula(paste(
        "~",
        stringr::str_replace_all(main_model_formula_str, "~", ""),
        "+", stringr::str_replace_all(nuisance_model_formula_str, "~", ""),
        "+", batch_col
      )),
      data = colData(subset_ccs)
    )
    if (Matrix::rankMatrix(full_model_matrix) < ncol(full_model_matrix)) {
      print(paste("Error: cannot correct for batches because the full model matrix is not full rank [model rank = ", Matrix::rankMatrix(full_model_matrix), "ncol =", ncol(full_model_matrix), "]"))
      print(colnames(full_model_matrix))
    } else {
      # main_model_formula_str = paste(main_model_formula_str, "+expt")
      nuisance_model_formula_str <- paste(nuisance_model_formula_str, "+", batch_col)
    }
  }

  main_model_formula_str_xxx <- stringr::str_replace_all(main_model_formula_str, "~", "")
  nuisance_model_formula_str_xxx <- stringr::str_replace_all(nuisance_model_formula_str, "~", "")

  full_model_formula_str <- paste("~", nuisance_model_formula_str_xxx, "+", main_model_formula_str_xxx)


  message(paste("\tformula:", full_model_formula_str))

  full_model_matrix <- Matrix::sparse.model.matrix(as.formula(full_model_formula_str), data = colData(subset_ccs))
  if (Matrix::rankMatrix(full_model_matrix) < ncol(full_model_matrix)) {
    print(paste("Error: full model matrix is not full rank [model rank = ", Matrix::rankMatrix(full_model_matrix), "ncol =", ncol(full_model_matrix), "]"))
    print(colnames(full_model_matrix))
    return(NA)
  }

  if (is.null(prior_state_transition_graph) == FALSE) {
    wt_prior_allowlist <- prior_state_transition_graph %>% igraph::as_data_frame()
  } else {
    wt_prior_allowlist <- NULL
  }

  genotype_ccm <- suppressMessages(suppressWarnings(new_cell_count_model(subset_ccs,
    main_model_formula_str = main_model_formula_str,
    # main_model_formula_str = "~ splines::ns(timepoint, knots=c(24, 30, 36)) + knockout",
    # main_model_formula_str = "~ as.factor(timepoint) + knockout",

    nuisance_model_formula_str = nuisance_model_formula_str,
    allowlist = edge_allowlist,
    denylist = edge_denylist,
    vhat_method = vhat_method,
    penalize_by_distance = penalize_by_distance,
    num_threads = num_threads,
    keep_ccs = keep_ccs,
    num_bootstraps = num_bootstraps,
    backend = backend,
    ftol_rel = ftol_rel,
    verbose = FALSE
  )))
  # FIXME: we should maybe be pulling out the sparsity_factor used for the WT prior model and using that here rather
  # than hardcoding?
  genotype_ccm <- select_model(genotype_ccm, sparsity_factor = sparsity_factor)

  # save information
  genotype_ccm@info$genotype <- genotype
  genotype_ccm@info$perturbation_col <- perturbation_col
  genotype_ccm@info$ctrl_ids <- ctrl_ids
  genotype_ccm@info$batch_col <- batch_col
  genotype_ccm@info$interval_col <- interval_col
  genotype_ccm@info$contrast_time_start <- contrast_time_start
  genotype_ccm@info$contrast_time_stop <- contrast_time_stop

  return(genotype_ccm)
}

#' wrapper function to easily plot output of fit_genotype_ccm
#' @param ccm a cell_count_model object
#' @export
make_contrast <- function(ccm, newdata = tibble(), adjust_q_values = FALSE) {
  if (nrow(newdata) > 0) {
    newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
    newdata_mt <- cross_join(tibble(knockout = TRUE), newdata)
  } else {
    newdata_wt <- tibble(knockout = FALSE)
    newdata_mt <- tibble(knockout = TRUE)
  }

  wt_cond <- estimate_abundances(ccm, newdata = newdata_wt)
  mt_cond <- estimate_abundances(ccm, newdata = newdata_mt)
  tbl <- compare_abundances(ccm, wt_cond, mt_cond, adjust_q_values = adjust_q_values)
  return(tbl)
}

#' Assemble Partition
#'
#' This function assembles a partition for a given cell dataset (CDS) by fitting wild-type (WT)
#' and mutant (MT) models, constructing state transition graphs, and assessing perturbation effects.
#'
#' @param cds A cell dataset object.
#' @param sample_group A string specifying the sample group column in the CDS.
#' @param cell_group A string specifying the cell group column in the CDS.
#' @param partition_name A string specifying the name of the partition (default: NULL).
#' @param main_model_formula_str A string specifying the main model formula (default: NULL).
#' @param start_time Numeric value specifying the start time for the analysis (default: 18).
#' @param stop_time Numeric value specifying the stop time for the analysis (default: 72).
#' @param interval_col A string specifying the column for time intervals (default: "timepoint").
#' @param nuisance_model_formula_str A string specifying the nuisance model formula (default: "~1").
#' @param ctrl_ids A vector of control IDs (default: NULL).
#' @param mt_ids A vector of mutant IDs (default: NULL).
#' @param sparsity_factor Numeric value specifying the sparsity factor (default: 0.01).
#' @param perturbation_col A string specifying the perturbation column (default: "gene_target").
#' @param batch_col A string specifying the batch column (default: "expt").
#' @param max_num_cells Numeric value specifying the maximum number of cells to include (default: NULL).
#' @param verbose Logical indicating whether to print verbose messages (default: FALSE).
#' @param keep_ccs Logical indicating whether to keep connected components (default: TRUE).
#' @param num_threads Numeric value specifying the number of threads to use (default: 1).
#' @param backend A string specifying the backend to use for optimization (default: "nlopt").
#' @param q_val Numeric value specifying the q-value threshold for perturbation effects (default: 0.1).
#' @param vhat_method A string specifying the method for variance estimation (default: "bootstrap").
#' @param num_bootstraps Numeric value specifying the number of bootstraps for variance estimation (default: 10).
#' @param newdata A tibble containing new data for predictions (default: tibble()).
#' @param edge_allowlist A list of edges to allow in the graph (default: NULL).
#' @param min_lfc Numeric value specifying the minimum log fold change for perturbation effects (default: 0).
#' @param links_between_components A character vector specifying the type of links between components
#'        (default: c("ctp", "none", "strongest-pcor", "strong-pcor")).
#' @param log_abund_detection_thresh Numeric value specifying the log abundance detection threshold (default: -5).
#' @param batches_excluded_from_assembly A vector of batch names to exclude from assembly (default: c()).
#' @param component_col A string specifying the column for components (default: "partition").
#' @param embryo_size_factors A vector of size factors for embryos (default: NULL).
#'
#' @return A tibble containing the results of the partition assembly, including WT and MT graphs,
#'         perturbation effects, and state graph plots.
#'
#' @details
#' The function performs the following steps:
#' - Prepares the CDS by adding subassembly group and cell state information.
#' - Fits a wild-type model and assembles a WT state transition graph.
#' - Fits mutant models and assembles MT state transition graphs.
#' - Assesses perturbation effects and constructs annotated graphs.
#' - Handles errors gracefully and returns NA for failed steps.
#'
#' @examples
#' # Example usage:
#' results <- assemble_partition(
#'   cds = my_cds,
#'   sample_group = "sample",
#'   cell_group = "cell_type",
#'   partition_name = "partition_1",
#'   main_model_formula_str = "~timepoint",
#'   start_time = 18,
#'   stop_time = 72
#' )
#'
#' @export
assemble_partition <- function(cds,
                               sample_group,
                               cell_group,
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
  colData(cds)$subassembly_group <- stringr::str_c(partition_name, colData(cds)[, cell_group], sep = "-")
  colData(cds)[["cell_state"]] <- as.character(colData(cds)[[cell_group]])
  # selected_colData = selected_colData %>% mutate(cell_state = paste0(partition_name, cell_state))
  selected_colData <- colData(cds) %>%
    as_tibble() %>%
    dplyr::select(cell, !!sym(sample_group), cluster, !!sym(cell_group), subassembly_group)

  selected_colData$cds_row_id <- colData(cds) %>%
    as.data.frame() %>%
    row.names()

  partition_results <- selected_colData %>% tidyr::nest(data = c(cds_row_id, cell, !!sym(sample_group), cluster, !!sym(cell_group), subassembly_group))

  # partition_results$cell_plot_state = list(plot_cells(cds, color_cells_by="cell_state"))
  # partition_results$cell_plot_time = list(plot_cells(cds, color_cells_by=interval_col))
  # partition_results$cell_plot_type = list(plot_cells(cds, color_cells_by="cell_type_sub"))

  if (length(unique(selected_colData[["cell_state"]])) <= 1) {
    partition_results$wt_graph <- list(NA)
    partition_results$mt_graph <- list(NA)
    partition_results$perturbation_effects <- list(NA)
    partition_results$wt_state_graph_plot <- list(NA)
    partition_results$mt_state_graph_plot <- list(NA)
    return(partition_results)
  }

  # Exclude user-specified experiments prior to assembly but after iterative UMAP & clustering:
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
        # stop("Error: fit_wt_model() failed")
        return(partition_results)
        # return(partition_results)
      }

      # FIXME: probably need to pass additional args here sometimes:
      wt_extant_cell_type_df <- get_extant_cell_types(wt_ccm, start_time, stop_time, interval_col = interval_col, newdata = newdata)

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
        num_time_breaks = num_time_breaks,
        # nuisance_model_formula_str = "~expt",
        links_between_components = links_between_components,
        ctrl_ids = ctrl_ids,
        edge_allowlist = edge_allowlist,
        sparsity_factor = sparsity_factor,
        perturbation_col = perturbation_col,
        component_col = component_col,
        verbose = verbose
      )

      if (is.null(wt_graph) == FALSE) {
        # partition_results$wt_ccm = list(wt_ccm)
        igraph::E(wt_graph)$assembly_group <- partition_name
        if (cell_group == "cell_state") {
          igraph::V(wt_graph)$name <- stringr::str_c(partition_name, igraph::V(wt_graph)$name, sep = "-")
        } else {
          igraph::V(wt_graph)$name <- igraph::V(wt_graph)$name
        }
        partition_results$wt_graph <- list(wt_graph)

        # partition_results$wt_state_graph_plot = list(plot_state_graph_annotations(wt_ccm, wt_graph,
        #                                                                           color_nodes_by = "timepoint",
        #                                                                           group_nodes_by="cell_type_sub",
        #                                                                           edge_weights = "support",
        #                                                                           hide_unlinked_nodes = FALSE))
      } else {
        partition_results$wt_graph <- list(NA)
        partition_results$wt_state_graph_plot <- list(NA)
      }

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


      message("Assembling mutant graphs...")
      mt_graph <- assemble_mt_graph(ref_ccs,
                                    wt_graph, 
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

      # partition_results$wt_ccm = list(wt_ccm)
      if (is.null(mt_graph) == FALSE) {
        igraph::E(mt_graph)$assembly_group <- partition_name

        if (cell_group == "cell_state") {
          igraph::V(mt_graph)$name <- stringr::str_c(partition_name, igraph::V(mt_graph)$name, sep = "-")
        } else {
          igraph::V(mt_graph)$name <- igraph::V(mt_graph)$name
        }
        merge_wt_graph_edges <- igraph::as_data_frame(wt_graph)
        merge_mt_graph_edges <- igraph::as_data_frame(mt_graph)

        merge_wt_graph_edges <- merge_wt_graph_edges %>% filter(to %in% merge_mt_graph_edges$to == FALSE) # exclude WT edges to nodes that have at least one mutant-supported parent

        merge_wt_graph_nodes <- igraph::as_data_frame(wt_graph, what = "vertices")
        merge_mt_graph_nodes <- igraph::as_data_frame(mt_graph, what = "vertices")

        if (cell_group == "cell_state") {
          merge_annotated_graph_nodes <- data.frame(name = stringr::str_c(partition_name, row.names(wt_ccm@ccs), sep = "-"))
        } else {
          merge_annotated_graph_nodes <- data.frame(name = row.names(wt_ccm@ccs))
        }

        merge_annotated_graph_nodes <- left_join(merge_annotated_graph_nodes, merge_mt_graph_nodes)

        mt_only <- setdiff(merge_mt_graph_edges %>% select(from, to), merge_wt_graph_edges %>%
          select(from, to)) %>% as.data.frame()


        # merge_wt_graph_edges doesn't have the assembly group name
        merge_graph_annotated <- left_join(merge_wt_graph_edges, merge_mt_graph_edges, by = c("from", "to", "assembly_group"))
        merge_graph_annotated <- merge_graph_annotated %>% select(-support)

        # BREAK
        merge_graph_annotated <- rbind(
          merge_graph_annotated,
          merge_mt_graph_edges %>% inner_join(mt_only)
        )
        merge_graph_annotated <- igraph::graph_from_data_frame(merge_graph_annotated, vertices = merge_annotated_graph_nodes)
        merge_graph_annotated <- platt:::break_cycles_in_state_transition_graph(merge_graph_annotated, "total_perturb_path_score_supporting")

        mt_graph <- merge_graph_annotated

        partition_results$mt_graph <- list(mt_graph)
        perturbation_effects <- perturb_models_tbl %>%
          dplyr::select(perturb_name, perturb_summary_tbl) %>%
          tidyr::unnest(perturb_summary_tbl)

        perturbation_effects$cell_group <- stringr::str_c(partition_name, perturbation_effects$cell_group, sep = "-")
        perturbation_effects <- perturbation_effects %>% tidyr::nest(perturb_summary_tbl = !perturb_name)
        partition_results$perturbation_effects <- list(perturbation_effects)

        perturbation_table <- perturb_models_tbl %>%
          dplyr::select(perturb_name, perturbation_table) %>%
          tidyr::unnest(perturbation_table)
        partition_results$perturbation_table <- list(perturbation_table)

        # this is failing because cell state is not matching
        # partition_results$mt_state_graph_plot = list(plot_state_graph_annotations(wt_ccm@ccs,
        #                                                                         mt_graph,
        #                                                                         label_nodes_by="global_cell_state",
        #                                                                         #color_nodes_by = "timepoint",
        #                                                                         group_nodes_by="cell_type_sub",
        #                                                                         label_edges_by="support_label",
        #                                                                         edge_weights = "num_perturbs_supporting",
        #                                                                         hide_unlinked_nodes = TRUE))
      } else {
        partition_results$mt_graph <- list(NA)
        partition_results$mt_state_graph_plot <- list(NA)
        partition_results$perturbation_effects <- list(NA)
        partition_results$perturbation_table <- list(NA)
      }
      partition_results$mt_graph_denylist <- list(NA)
      partition_results$mt_graph_denylist_plot <- list(NA)
    },
    error = function(e) {
      print(e)
      partition_results$wt_graph <- list(NA)
      partition_results$mt_graph <- list(NA)
      partition_results$perturbation_effects <- list(NA)
      partition_results$wt_state_graph_plot <- list(NA)
      partition_results$mt_state_graph_plot <- list(NA)
      return(partition_results)
    }
  )

  return(partition_results)
}


#' Fit Wild Type Model
#'
#' This function fits a wild type (WT) model to a cell dataset (CDS) by estimating
#' cell count dynamics over time and accounting for nuisance variables. It supports
#' various customization options for model fitting, including user-defined formulas,
#' size factors, and control IDs.
#'
#' @param cds A cell dataset (CDS) object containing cell data.
#' @param sample_group A string specifying the column in `colData(cds)` that defines sample groups.
#' @param cell_group A string specifying the column in `colData(cds)` that defines cell groups.
#' @param main_model_formula_str A string specifying the main model formula. If `NULL`, it will be generated automatically.
#' @param num_time_breaks An integer specifying the number of time breaks for the main model formula. Default is 4.
#' @param nuisance_model_formula_str A string specifying the nuisance model formula. Default is `"~1"`.
#' @param ctrl_ids A vector of control IDs. If `NULL`, control IDs are inferred from `perturbation_col`.
#' @param sparsity_factor A numeric value for sparsity factor used in model selection. Default is 1.
#' @param vhat_method A string specifying the method for estimating variance (`"bootstrap"` by default).
#' @param interval_col A string specifying the column in `colData(cds)` that defines time intervals. Default is `"timepoint"`.
#' @param perturbation_col A string specifying the column in `colData(cds)` that defines perturbation groups. Default is `"knockout"`.
#' @param batch_col A string specifying the column in `colData(cds)` that defines batch groups. Default is `"expt"`.
#' @param start_time A numeric value specifying the start time for the model. If `NULL`, it is inferred from the data.
#' @param stop_time A numeric value specifying the stop time for the model. If `NULL`, it is inferred from the data.
#' @param interval_step A numeric value specifying the step size for time intervals. Default is 2.
#' @param log_abund_detection_thresh A numeric threshold for log abundance detection. Default is -5.
#' @param keep_ccs A logical value indicating whether to retain the cell count set (CCS). Default is `TRUE`.
#' @param q_val A numeric value specifying the q-value threshold for significance. Default is 0.1.
#' @param edge_allowlist A vector of edges to allow in the model. Default is `NULL`.
#' @param edge_denylist A vector of edges to deny in the model. Default is `NULL`.
#' @param base_penalty A numeric value specifying the base penalty for model selection. Default is 1.
#' @param keep_cds A logical value indicating whether to retain the CDS. Default is `TRUE`.
#' @param verbose A logical value indicating whether to print verbose messages. Default is `FALSE`.
#' @param num_threads An integer specifying the number of threads to use. Default is 1.
#' @param backend A string specifying the backend for model fitting. Default is `"nlopt"`.
#' @param penalize_by_distance A logical value indicating whether to penalize by distance. Default is `TRUE`.
#' @param embryo_size_factors A named vector of size factors for embryos. Default is `NULL`.
#' @param batches_excluded_from_assembly A vector of batch IDs to exclude from assembly. Default is an empty vector.
#' @param ... Additional arguments passed to the underlying model fitting functions.
#'
#' @return A fitted wild type cell count model object, or `NULL` if no control cells are available.
#'
#' @details
#' The function first subsets the CDS to include only control cells based on the `perturbation_col`
#' and `ctrl_ids`. It then constructs a cell count set (CCS) and fits a model using the specified
#' formulas and parameters. If no control cells are available or only a single cell group is present,
#' the function returns `NULL`.
#'
#' @examples
#' # Example usage:
#' wt_model <- fit_wt_model(
#'   cds = my_cds,
#'   sample_group = "sample",
#'   cell_group = "cell_type",
#'   interval_col = "timepoint",
#'   perturbation_col = "knockout",
#'   batch_col = "batch"
#' )
#' @export
fit_wt_model <- function(cds,
                         sample_group,
                         cell_group,
                         main_model_formula_str = NULL,
                         num_time_breaks = 4,
                         nuisance_model_formula_str = "~1",
                         ctrl_ids = NULL,
                         sparsity_factor = 1,
                         vhat_method = "bootstrap",
                         interval_col = "timepoint",
                         perturbation_col = "knockout",
                         batch_col = "expt",
                         start_time = NULL,
                         stop_time = NULL,
                         log_abund_detection_thresh = -5,
                         keep_ccs = TRUE,
                         edge_allowlist = NULL,
                         edge_denylist = NULL,
                         base_penalty = 1,
                         keep_cds = TRUE,
                         verbose = FALSE,
                         num_threads = 1,
                         backend = "nlopt",
                         penalize_by_distance = TRUE,
                         embryo_size_factors = NULL,
                         batches_excluded_from_assembly = c(),
                         include_time_in_nuisance = FALSE,
                         ...) {
  if (is.null(ctrl_ids)) {
    ctrl_ids <- unique(colData(cds)[[perturbation_col]])
    ctrl_ids <- ctrl_ids[grepl("wt|ctrl", ctrl_ids)]
  }

  wt_cds <- cds[, colData(cds)[[perturbation_col]] %in% ctrl_ids]

  if (is.null(batch_col) == FALSE) {
    wt_cds <- wt_cds[, colData(wt_cds)[[batch_col]] %in% batches_excluded_from_assembly == FALSE]
  }


  if (ncol(wt_cds) == 0) {
    message("No control cells. Skipping...")
    return(NULL)
  }

  timepoints <- as.numeric(unique(colData(wt_cds)[[interval_col]]))
  timepoints <- timepoints[!is.na(timepoints)]

  if (is.null(start_time)) {
    start_time <- min(timepoints)
  }
  if (is.null(stop_time)) {
    stop_time <- max(timepoints)
  }

  wt_ccs <- new_cell_count_set(wt_cds,
    sample_group = sample_group,
    cell_group = cell_group,
    keep_cds = keep_cds,
    norm_method = "size_factors"
  )

  if (is.null(embryo_size_factors) == FALSE) {
    message("Using user-supplied size factors")
    colData(wt_ccs)$Size_Factor <- embryo_size_factors[colnames(wt_ccs)]
  }

  num_cell_groups <- nrow(wt_ccs)
  if (num_cell_groups <= 1) {
    stop("Only a single cell group. Skipping...")
  }

  # make this any column
  if (is.null(batch_col) == FALSE) {
    if (length(unique(colData(wt_ccs)[[batch_col]])) > 1) {
      # main_model_formula_str = paste(main_model_formula_str, "+ expt")
      nuisance_model_formula_str <- paste(nuisance_model_formula_str, "+", batch_col)
      colData(wt_ccs)[[batch_col]] <- as.factor(colData(wt_ccs)[[batch_col]])
    }
  }

  if (is.null(main_model_formula_str)) {
    main_model_formula_str <- build_interval_formula(wt_ccs,
      interval_var = interval_col,
      interval_start = start_time,
      interval_stop = stop_time,
      num_breaks = num_time_breaks
    )

    main_model_formula_str_xxx <- stringr::str_replace_all(main_model_formula_str, "~", "")
    nuisance_model_formula_str_xxx <- stringr::str_replace_all(nuisance_model_formula_str, "~", "")
    full_model_formula_str <- paste("~", nuisance_model_formula_str_xxx, "+", main_model_formula_str_xxx)

    # make these formulas the same
    if (include_time_in_nuisance) {
      nuisance_model_formula_str <- full_model_formula_str
    }

    message(paste("Fitting wild type model with main effects:", full_model_formula_str))
    message(paste("Nuisance effects:", nuisance_model_formula_str))
  }

  # undebug(new_cell_count_model)
  wt_ccm <- new_cell_count_model(wt_ccs,
    main_model_formula_str = full_model_formula_str,
    nuisance_model_formula_str = nuisance_model_formula_str,
    vhat_method = vhat_method,
    allowlist = edge_allowlist,
    denylist = edge_denylist,
    base_penalty = base_penalty,
    num_threads = num_threads,
    keep_ccs = keep_ccs,
    backend = backend,
    verbose = verbose,
    penalize_by_distance = penalize_by_distance,
    # covariance_type="spherical",
    ...
  )

  wt_ccm <- select_model(wt_ccm, criterion = "EBIC", sparsity_factor = sparsity_factor)

  return(wt_ccm)
}


#' Assemble Wild-Type State Transition Graph
#'
#' This function constructs a state transition graph for wild-type (WT) cells
#' based on a provided cell count model (CCM) and other parameters. It processes
#' the input data, filters for control IDs, and generates a graph representing
#' state transitions over time.
#'
#' @param cds A CellDataSet object containing single-cell data.
#' @param wt_ccm A cell count model object for wild-type cells.
#' @param sample_group A grouping variable for samples (not used directly in this function).
#' @param cell_group A grouping variable for cells (not used directly in this function).
#' @param newdata A tibble containing new data for predictions (default: empty tibble).
#' @param main_model_formula_str A string specifying the main model formula (default: NULL).
#' @param num_time_breaks Number of breaks for discretizing continuous variables (default: 4).
#' @param nuisance_model_formula_str A string specifying the nuisance model formula (default: "~1").
#' @param ctrl_ids A vector of control IDs (default: NULL, automatically inferred).
#' @param mt_ids A vector of mitochondrial IDs (default: NULL).
#' @param sparsity_factor A numeric factor for sparsity adjustment (default: 1).
#' @param vhat_method Method for estimating variance ("bootstrap" by default).
#' @param interval_col Column name representing time intervals (default: "timepoint").
#' @param perturbation_col Column name representing perturbation groups (default: "knockout").
#' @param start_time Start time for the state transition graph (default: NULL, inferred from data).
#' @param stop_time Stop time for the state transition graph (default: NULL, inferred from data).
#' @param interval_step Step size for time intervals (default: 2).
#' @param links_between_components Method for linking components in the graph
#'        (default: c("ctp", "none", "strongest-pcor", "strong-pcor")).
#' @param log_abund_detection_thresh Log abundance detection threshold (default: -5).
#' @param q_val Q-value threshold for significance (default: 0.1).
#' @param break_cycles Logical indicating whether to break cycles in the graph (default: TRUE).
#' @param edge_allowlist A list of edges to allow in the graph (default: NULL).
#' @param edge_denylist A list of edges to deny in the graph (default: NULL).
#' @param component_col Column name representing components (default: "partition").
#' @param verbose Logical indicating whether to print verbose output (default: FALSE).
#'
#' @return A state transition graph object representing the transitions between
#'         states in the wild-type cell population.
#'
#' @details
#' The function filters the input dataset to include only control IDs and valid
#' time intervals. It then constructs a state transition graph using the provided
#' cell count model and parameters. If `break_cycles` is TRUE, cycles in the graph
#' are removed to ensure acyclic transitions.
#'
#' @examples
#' # Example usage:
#' wt_graph <- assemble_wt_graph(cds, wt_ccm, sample_group, cell_group)
#'
#' @export
assemble_wt_graph <- function(cds,
                              wt_ccm,
                              sample_group,
                              cell_group,
                              newdata = tibble(),
                              main_model_formula_str = NULL,
                              nuisance_model_formula_str = "~1",
                              ctrl_ids = NULL,
                              mt_ids = NULL,
                              sparsity_factor = 1,
                              vhat_method = "bootstrap",
                              interval_col = "timepoint",
                              perturbation_col = "knockout",
                              start_time = NULL,
                              stop_time = NULL,
                              interval_step = 2,
                              links_between_components = c("ctp", "none", "strongest-pcor", "strong-pcor"),
                              log_abund_detection_thresh = -5,
                              q_val = 0.1,
                              min_interval = 4,
                              max_interval = 24,
                              min_pathfinding_lfc = 0,
                              break_cycles = TRUE,
                              edge_allowlist = NULL,
                              edge_denylist = NULL,
                              component_col = "partition",
                              verbose = FALSE) {
  if (is.null(ctrl_ids)) {
    ctrl_ids <- unique(colData(cds)[[perturbation_col]])
    ctrl_ids <- ctrl_ids[grepl("wt|ctrl", ctrl_ids)]
  }

  wt_cds <- cds[, colData(cds)[[perturbation_col]] %in% ctrl_ids]
  wt_cds <- wt_cds[, !is.na(colData(wt_cds)[[interval_col]])]

  timepoints <- unique(colData(wt_cds)[[interval_col]])
  timepoints <- timepoints[!is.na(timepoints)]

  if (is.null(start_time)) {
    start_time <- min(timepoints)
  }
  if (is.null(stop_time)) {
    stop_time <- max(timepoints)
  }

  if (is.null(wt_ccm) || is.na(wt_ccm)) {
    stop("No cell count model. Skipping.")
  }

  if (nrow(wt_ccm@ccs) <= 1) {
    stop("Model has only a single cell type. Skipping.")
  }

  wt_state_transition_graph <- assemble_timeseries_transitions(wt_ccm,
    start_time = start_time,
    stop_time = stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    min_interval = min_interval,
    max_interval = max_interval,
    min_pathfinding_lfc = min_pathfinding_lfc,
    log_abund_detection_thresh = log_abund_detection_thresh,
    q_val = q_val,
    links_between_components = links_between_components,
    edge_allowlist = edge_allowlist,
    edge_denylist = edge_denylist,
    components = component_col,
    newdata = newdata
  )
  if (break_cycles) {
    print("breaking cycles in control timeseries graph...")
    wt_state_transition_graph <- platt:::break_cycles_in_state_transition_graph(wt_state_transition_graph, "support")
  }

  return(wt_state_transition_graph)
}



# FIXME: allow using WT graph as a prior?
#' Fit Models for Multi-Timepoint Perturbation Analysis
#'
#' This function fits models for analyzing multi-timepoint perturbation data
#' using a cell dataset (cds). It allows for the specification of various
#' parameters to customize the analysis, including the main model formula,
#' nuisance model formula, and perturbation-specific settings.
#'
#' @param cds A cell dataset object (cds) containing the data to be analyzed.
#' @param sample_group A character string specifying the sample grouping variable.
#' @param cell_group A character string specifying the cell grouping variable.
#' @param main_model_formula_str A character string specifying the main model formula. Default is NULL.
#' @param num_time_breaks An integer specifying the number of time breaks for the analysis. Default is 3.
#' @param nuisance_model_formula_str A character string specifying the nuisance model formula. Default is "~1".
#' @param ctrl_ids A vector of control identifiers. Default is NULL.
#' @param mt_ids A vector of multi-timepoint identifiers. Default is NULL.
#' @param sparsity_factor A numeric value for sparsity adjustment. Default is 1.
#' @param vhat_method A character string specifying the method for variance estimation. Default is "bootstrap".
#' @param interval_col A character string specifying the column name for time intervals. Default is "timepoint".
#' @param perturbation_col A character string specifying the column name for perturbation identifiers. Default is "knockout".
#' @param batch_col A character string specifying the column name for batch identifiers. Default is "expt".
#' @param newdata A tibble containing new data for prediction. Default is an empty tibble.
#' @param start_time A numeric value specifying the start time for the analysis. Default is NULL.
#' @param stop_time A numeric value specifying the stop time for the analysis. Default is NULL.
#' @param interval_step An integer specifying the step size for time intervals. Default is 2.
#' @param log_abund_detection_thresh A numeric threshold for log abundance detection. Default is -5.
#' @param q_val A numeric value specifying the q-value threshold for significance. Default is 0.1.
#' @param edge_allowlist A vector of edges to allow in the analysis. Default is NULL.
#' @param edge_denylist A vector of edges to deny in the analysis. Default is NULL.
#' @param keep_cds A logical value indicating whether to keep the cds object in the output. Default is TRUE.
#' @param keep_ccs A logical value indicating whether to keep the cell count set (ccs) in the output. Default is TRUE.
#' @param verbose A logical value indicating whether to print verbose messages. Default is FALSE.
#' @param num_threads An integer specifying the number of threads to use. Default is 1.
#' @param backend A character string specifying the backend to use for optimization. Default is "nlopt".
#' @param penalize_by_distance A logical value indicating whether to penalize by distance in the analysis. Default is TRUE.
#' @param independent_spline_for_ko A logical value indicating whether to use independent splines for knockouts. Default is TRUE.
#' @param num_bootstraps An integer specifying the number of bootstraps for variance estimation. Default is 10.
#' @param embryo_size_factors A named vector of size factors for embryos. Default is NULL.
#' @param batches_excluded_from_assembly A vector of batch identifiers to exclude from the analysis. Default is an empty vector.
#'
#' @return A tibble containing the fitted perturbation models and associated metadata.
#'
#' @details This function performs multi-timepoint perturbation analysis by
#' fitting models to cell count data. It supports various customization options
#' for handling batch effects, time intervals, and perturbation-specific settings.
#' The function also allows for the inclusion of user-supplied size factors and
#' the exclusion of specific batches from the analysis.
#'
#' @examples
#' # Example usage:
#' # fit_mt_models(cds, sample_group = "sample", cell_group = "cell_type")
#'
#' @export
fit_mt_models <- function(cds,
                          sample_group,
                          cell_group,
                          main_model_formula_str = NULL,
                          num_time_breaks = 3,
                          nuisance_model_formula_str = "~1",
                          ctrl_ids = NULL,
                          mt_ids = NULL,
                          sparsity_factor = 1,
                          vhat_method = "bootstrap",
                          interval_col = "timepoint",
                          perturbation_col = "knockout",
                          batch_col = "expt",
                          newdata = tibble(),
                          start_time = NULL,
                          stop_time = NULL,
                          interval_step = 2,
                          log_abund_detection_thresh = -5,
                          q_val = 0.1,
                          edge_allowlist = NULL,
                          edge_denylist = NULL,
                          keep_cds = TRUE,
                          keep_ccs = TRUE,
                          verbose = FALSE,
                          num_threads = 1,
                          backend = "nlopt",
                          penalize_by_distance = TRUE,
                          independent_spline_for_ko = TRUE,
                          num_bootstraps = 10,
                          embryo_size_factors = NULL,
                          batches_excluded_from_assembly = c()) {
  if (!is.null(mt_ids)) {
    cds <- cds[, replace_na(colData(cds)[[perturbation_col]] %in% c(ctrl_ids, mt_ids), F)]
  }

  cds <- cds[, colData(cds)[[batch_col]] %in% batches_excluded_from_assembly == FALSE]


  ccs <- new_cell_count_set(cds,
    sample_group = sample_group,
    cell_group = cell_group,
    keep_cds = keep_cds,
    norm_method = "size_factors"
  )

  timepoints <- as.numeric(unique(colData(cds)[[interval_col]]))
  timepoints <- timepoints[!is.na(timepoints)]

  if (is.null(start_time)) {
    start_time <- min(timepoints)
  }
  if (is.null(stop_time)) {
    stop_time <- max(timepoints)
  }

  ccs@cds_coldata[["perturb_name"]] <- ccs@cds_coldata[[perturbation_col]]

  if (is.null(embryo_size_factors) == FALSE) {
    message("Using user-supplied size factors")
    colData(ccs)$Size_Factor <- embryo_size_factors[colnames(ccs)]
  }

  num_cell_groups <- nrow(ccs)
  if (num_cell_groups <= 1) {
    stop("Only a single cell group. Skipping...")
  }

  perturb_df <- ccs@cds_coldata %>%
    as_tibble() %>%
    dplyr::select(perturb_name) %>%
    filter(!perturb_name %in% ctrl_ids) %>%
    distinct()
  perturb_models_tbl <- perturb_df %>%
    dplyr::mutate(perturb_time_window = purrr::map(
      .f = purrr::possibly(get_time_window, NA_real_),
      .x = perturb_name,
      ccs = ccs,
      interval_col = interval_col,
      perturbation_col = perturbation_col
    )) %>%
    dplyr::mutate(perturb_ccm = purrr::map(
      .f = purrr::possibly(fit_genotype_ccm, NA_real_),
      .x = perturb_name,
      ccs,
      # prior_state_transition_graph = wt_state_transition_graph,
      ctrl_ids = ctrl_ids,
      interval_col = interval_col,
      perturbation_col = perturbation_col,
      num_time_breaks = num_time_breaks,
      batch_col = batch_col,
      # assembly_time_start=start_time,
      # assembly_time_stop=stop_time,
      keep_ccs = keep_ccs,
      edge_allowlist = edge_allowlist,
      edge_denylist = edge_denylist,
      penalize_by_distance = penalize_by_distance,
      independent_spline_for_ko = independent_spline_for_ko,
      num_threads = num_threads,
      vhat_method = vhat_method,
      backend = backend,
      num_bootstraps = num_bootstraps
    ))

  return(perturb_models_tbl)
}

#' assembles a graph using the perturbation data
#' @export
assemble_mt_graph <- function(ref_ccs,
                              wt_graph,
                              perturb_models_tbl,
                              interval_col = "timepoint",
                              # perturbation_col = "knockout",
                              start_time = NULL,
                              stop_time = NULL,
                              interval_step = 2,
                              links_between_components = "none",
                              log_abund_detection_thresh = -5,
                              q_val = 0.1,
                              newdata = tibble(),
                              break_cycles = TRUE,
                              component_col = "partition",
                              edge_allowlist = NULL,
                              edge_denylist = NULL,
                              verbose = FALSE) {
  # if (is.null(wt_ccm) || is.na(wt_ccm)) {
  #   stop("No control timeseries cell count model. Skipping.")
  # }

  wt_cds <- ref_ccs@cds

  timepoints <- unique(colData(wt_cds)[[interval_col]])
  timepoints <- timepoints[!is.na(timepoints)]

  if (is.null(start_time)) {
    start_time <- min(timepoints)
  }
  if (is.null(stop_time)) {
    stop_time <- max(timepoints)
  }

  if (nrow(ref_ccs) <= 1) {
    stop("Control timeseries model has only a single cell type. Skipping.")
  }

  perturb_models_tbl <- perturb_models_tbl %>%
    filter(is.na(perturb_ccm) == FALSE)

  if (nrow(perturb_models_tbl) == 0) {
    stop("No valid perturbation models.")
  }

  mutant_supergraph <- assemble_transition_graph_from_perturbations(
    ref_ccs,
    wt_graph,
    perturb_models_tbl,
    start_time = start_time,
    stop_time = stop_time,
    # perturbation_col = perturbation_col,
    interval_col = interval_col,
    interval_step = interval_step,
    log_abund_detection_thresh = log_abund_detection_thresh,
    q_val = q_val,
    newdata = newdata,
    links_between_components = links_between_components,
    edge_allowlist = edge_allowlist,
    edge_denylist = edge_denylist,
    components = component_col,
    verbose = verbose
  )
  if (break_cycles & inherits(mutant_supergraph, "igraph")) {
    print("breaking cycles in perturbation graph...")
    mutant_supergraph <- platt:::break_cycles_in_state_transition_graph(mutant_supergraph, "total_perturb_path_score_supporting")
  }

  return(mutant_supergraph)
}


default_resolution_fun <- function(num_cells, min_res = 5e-6, max_res = 1e-5, max_num_cells = NULL) {
  if (is.null(max_num_cells)) {
    max_num_cells <- num_cells
  }
  resolution <- approxfun(c(0, log10(max_num_cells)), c(min_res, max_res))(log10(num_cells))
  reflected_resolution <- (max_res - resolution) + min_res
  return(reflected_resolution)
}


#' @export
subcluster_cds <- function(cds,
                           recursive_subcluster = FALSE,
                           partition_name = NULL,
                           num_dim = NULL,
                           max_components = 3,
                           resolution_fun = NULL,
                           max_num_cells = NULL,
                           min_res = 5e-6,
                           max_res = 1e-5,
                           cluster_k = 20,
                           num_threads = 1) {
  message("Clustering all cells")

  if (is.null(max_num_cells)) {
    max_num_cells <- ncol(cds)
  }

  if (is.null(resolution_fun)) {
    resolution_fun <- function(num_cells) {
      resolution <- approxfun(c(0, log10(max_num_cells)), c(min_res, max_res))(log10(num_cells))
      reflected_resolution <- (max_res - resolution) + min_res
      return(reflected_resolution)
    }
    # resolution_fun(min_res)
  }

  partition_resolution <- resolution_fun(ncol(cds))
  message(paste("Clustering", ncol(cds), "cells at resolution =", partition_resolution))
  cds <- monocle3::cluster_cells(cds, resolution = partition_resolution, k = cluster_k)
  colData(cds)$cluster <- monocle3::clusters(cds)
  colData(cds)$res <- partition_resolution
  partitions <- unique(monocle3::partitions(cds))

  # if (length(partitions) > 1) {
  #   message(paste0("This could be split up further into more partitions. Num partitions = ", length(partitions)))
  # }

  # if we want to recursively subcluster
  if (length(partitions) > 1 & recursive_subcluster) {
    partition_res <- lapply(partitions, function(partition) {
      if (is.null(partition_name)) {
        next_partition_name <- partition
      } else {
        next_partition_name <- paste0(partition_name, "_", partition)
      }
      print(next_partition_name)
      message(paste("Constructing sub-UMAP for partition", next_partition_name))
      cds <- cds[, monocle3::partitions(cds) == partition]

      RhpcBLASctl::blas_set_num_threads(num_threads)
      RhpcBLASctl::omp_set_num_threads(num_threads)

      cds <- suppressMessages(suppressWarnings(preprocess_cds(cds))) %>%
        align_cds(residual_model_formula_str = "~log.n.umi") %>%
        suppressMessages(suppressWarnings(reduce_dimension(
          max_components = max_components,
          preprocess_method = "Aligned",
          umap.fast_sgd = TRUE,
          cores = num_threads
        )))

      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)

      cds <- subcluster_cds(cds,
        partition_name = next_partition_name,
        recursive_subcluster = recursive_subcluster,
        num_dim = num_dim,
        max_components = max_components,
        resolution_fun = resolution_fun,
        max_num_cells = max_num_cells,
        min_res = min_res,
        max_res = max_res,
        cluster_k = cluster_k
      )

      # otherwise sometimes the matrix columns dim don't match when
      # trying to keep the reduced dims in the combine
      reducedDims(cds)$PCA <- NULL
      reducedDims(cds)$Aligned <- NULL
      cds
    })
    # undebug(combine_cds)

    cds <- combine_cds(partition_res, keep_reduced_dims = T)
  } else {
    # save the umap coordinates
    partion_umap_coords <- reducedDims(cds)[["UMAP"]]
    num_components <- dim(partion_umap_coords)[[2]]
    for (i in 1:num_components) {
      name <- paste0("partition_umap", num_components, "d_", i)
      colData(cds)[[name]] <- reducedDims(cds)[["UMAP"]][, i]
    }

    # if it doesn't have an error, then it hasn't gone through combine cds
    has_partitions <- tryCatch(monocle3::partitions(cds),
      error = function(e) {
        NULL
      }
    )

    if (!is.null(has_partitions)) {
      if (is.null(partition_name)) {
        colData(cds)$partition <- monocle3::partitions(cds)
        colData(cds)$cluster <- monocle3::clusters(cds)
        colData(cds)$cell_state <- monocle3::clusters(cds)
      } else {
        colData(cds)$partition <- as.character(partition_name) #
        # colData(cds)$partition = paste0(partition_name, "_",  monocle3::partitions(cds))
        colData(cds)$cluster <- monocle3::clusters(cds)
        colData(cds)$cell_state <- paste0(partition_name, "-", monocle3::clusters(cds))
      }
    }
  }

  return(cds)
}


# Need to convert the cluster IDs in each graph to cell_state IDs (which are what we'll use in the final model)
convert_graph_ids <- function(state_graph, partition, cluster_to_state_id_tbl, support_col = "support") {
  if (is.null(state_graph)) {
    return(NA)
  }
  # if (is.na(state_graph))
  #  return (NA)
  if (igraph::is.igraph(state_graph) == FALSE) {
    return(NA)
  }

  # cluster_to_state_id_tbl = cluster_to_state_id_tbl %>% filter(grepl(paste("^",partition,"-", sep=""), cell_state))

  # cluster_to_state_id_tbl = colData(cds)[,c("cluster", "cell_state")] %>% as_tibble() %>% distinct()
  state_graph <- igraph::as_data_frame(state_graph)
  state_graph$to <- as.character(state_graph$to)
  state_graph$from <- as.character(state_graph$from)
  state_graph <- left_join(state_graph, cluster_to_state_id_tbl, by = c("from" = "cluster")) %>%
    dplyr::select(-from) %>%
    dplyr::rename(from = cell_state)
  state_graph <- left_join(state_graph, cluster_to_state_id_tbl, by = c("to" = "cluster")) %>%
    dplyr::select(-to) %>%
    dplyr::rename(to = cell_state)
  state_graph <- state_graph[, c("from", "to", support_col)]
  state_graph <- igraph::graph_from_data_frame(state_graph, directed = TRUE)
}

# if you have saved partition coords, can reconstruct sub cdss
# still need to cluster them

get_partition_cds <- function(cds,
                              partition_id,
                              partition_col = "pcor_cluster",
                              umap_prefix = "partition_umap3d_") {
  coldata <- colData(cds) %>% as.data.frame()
  cell_ids <- rownames(coldata %>% filter(!!sym(partition_col) == partition_id))
  # subset to just that partition
  partition_cds <- cds[, cell_ids]
  # change to sub umap coords
  curr_umap_matrix <- reducedDim(partition_cds, type = "UMAP")
  umap_dims <- dim(curr_umap_matrix)[2]
  umap_names <- paste0(umap_prefix, 1:umap_dims)

  new_umap_matrix <- coldata %>%
    select(all_of(umap_names)) %>%
    as.matrix()
  reducedDims(partition_cds)[["UMAP"]] <- new_umap_matrix[rownames(curr_umap_matrix), ]


  return(partition_cds)
}


adjust_time_stage <- function(cds) {
  staging_df <- colData(cds) %>%
    as.data.frame() %>%
    select(cell, embryo, expt, timepoint, mean_nn_time) %>%
    as_tibble()
  staging_df <- staging_df %>% mutate(timepoint = as.numeric(timepoint))

  staging_model <- lm(mean_nn_time ~ as.numeric(timepoint) * expt, data = staging_df)
  staging_df$predicted_timepoint <- predict(staging_model, newdata = staging_df)

  colData(cds)$adjusted_timepoint <- staging_df$predicted_timepoint
  return(cds)
}

assign_cell_states <- function(comb_res) {
  cell_state_assignments <- comb_res %>%
    select(data) %>%
    tidyr::unnest(data)
  cell_state_assignments$cluster <- as.character(cell_state_assignments$cluster)
  cell_state_assignments <- cell_state_assignments %>%
    distinct() %>%
    as.data.frame(stringsAsFactors = FALSE)
  # row.names(cell_state_assignments) = cell_state_assignments$cell
  # comb_res = sub_partition_cds(comb_sample_cds)

  cell_state_assignments <- cell_state_assignments %>% mutate(partition = stringr::str_split_fixed(cell_state, "-", n = 2)[, 1])
  row.names(cell_state_assignments) <- cell_state_assignments$cds_row_id

  return(cell_state_assignments)
}

# this cds was clustered/has partitions, but i only want to run assembly on 1
# partition
# wrapper function for purrr use cases
run_assembly <- function(cds,
                         partition_group,
                         interval_col = "timepoint",
                         recluster = FALSE,
                         recursive_subcluster = FALSE,
                         ...) {
  part_cds <- get_partition_cds(cds, partition_group)
  partition_results <- assemble_partition(part_cds,
    recluster = recluster,
    recursive_subcluster = recursive_subcluster,
    interval_col = interval_col,
    partition_name = partition_group,
    ...
  )
  return(partition_results)
}



# wrapper function to split it up
# because i can't hold comb_cds in memory

run_partition_assembly <- function(wt_cds,
                                   mt_cds,
                                   partition,
                                   coemb = FALSE,
                                   resolution_fun = NULL,
                                   max_num_cells = NULL,
                                   min_res = 5e-6,
                                   max_res = 1e-5,
                                   ...) {
  wt_i_cds <- wt_cds[, colData(wt_cds)$partition == partition]
  mt_i_cds <- mt_cds[, replace_na(colData(mt_cds)$partition == partition, F)]

  # my mt cds is wrong

  comb_i_cds <- combine_cds(list(wt_i_cds, mt_i_cds), keep_reduced_dims = T)

  if (coemb) {
    comb_i_cds <- comb_i_cds %>%
      preprocess_cds(num_dim = 50) %>%
      align_cds(residual_model_formula_str = "~log.n.umi") %>%
      reduce_dimension(max_components = 3)
  } else {
    comb_i_cds <- get_partition_cds(comb_i_cds, partition)
  }

  # remove outliers
  # comb_i_cds = drop_outlier_cells(comb_i_cds)

  # recursively sub cluster
  comb_i_cds <- subcluster_cds(comb_i_cds,
    recursive_subcluster = T,
    partition_name = partition,
    num_dim = NULL,
    max_components = 3,
    resolution_fun = resolution_fun,
    max_num_cells = max_num_cells,
    min_res = min_res,
    max_res = max_res
  )

  partitions <- unique(colData(comb_i_cds)$partition)

  comb_res <- lapply(partitions, function(p) {
    comb_p_cds <- get_partition_cds(comb_i_cds, p)

    partition_results <- assemble_partition(comb_p_cds,
      recluster = T,
      recursive_subcluster = F,
      interval_col = "timepoint",
      partition_name = p,
      max_num_cells = max_num_cells,
      min_res = min_res,
      max_res = max_res,
      ...
    )
    partition_results
  })

  comb_res <- bind_rows(comb_res)

  return(comb_res)
}


run_cds_assembly <- function(cds,
                             resolution_fun = NULL,
                             max_num_cells = NULL,
                             min_res = 5e-6,
                             max_res = 1e-5,
                             recluster = TRUE,
                             ...) {
  partitions <- unique(colData(cds)$partition)

  comb_res <- lapply(partitions, function(p) {
    comb_p_cds <- get_partition_cds(cds, p)

    partition_results <- assemble_partition(comb_p_cds,
      recluster = recluster,
      recursive_subcluster = F,
      interval_col = "timepoint",
      cell_group = "cell_state",
      partition_name = p,
      max_num_cells = max_num_cells,
      min_res = min_res,
      max_res = max_res,
      ...
    )
    partition_results
  })

  comb_res <- bind_rows(comb_res)

  return(comb_res)
}

#'
#' @export
collect_genotype_effects <- function(ccm, newdata = tibble()) {
  if (nrow(newdata) > 0) {
    newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
    newdata_mt <- cross_join(tibble(knockout = FALSE), newdata)
  } else {
    newdata_mt <- newdata_wt <- tibble(knockout = FALSE)
  }

  control_abund <- estimate_abundances(ccm, newdata_wt)
  knockout_abund <- estimate_abundances(ccm, newdata_mt)
  genotype_comparison_tbl <- compare_abundances(ccm, control_abund, knockout_abund)
}


#' This function looks at the effects of each perturbation to assemble a list
#' of genes required by each cell type
#' @export
categorize_genetic_requirements <- function(perturb_ccm_tbl, state_graph) {
  lost_cell_groups <- perturb_ccm_tbl %>%
    tidyr::unnest(perturb_summary_tbl) %>%
    group_by(perturb_name) %>%
    arrange(loss_when_present_p_value) %>%
    filter(is_lost_when_present) %>%
    ungroup()

  lost_cell_groups <- lost_cell_groups %>%
    dplyr::select(perturb_name, cell_group, loss_when_present_p_value) %>%
    group_by(perturb_name) %>%
    tidyr::nest(lost_cell_groups = c(cell_group, loss_when_present_p_value))

  get_dir_losses <- function(lcgs, sg) {
    lost_cell_groups <- lcgs$cell_group
    lost_cell_groups_in_graph <- intersect(igraph::V(sg)$name, lost_cell_groups)
    lost_cell_groups_not_in_graph <- setdiff(lost_cell_groups, igraph::V(sg)$name)
    lost_subgraph <- igraph::subgraph(sg, lost_cell_groups_in_graph)
    if (length(igraph::V(lost_subgraph)) > 0) {
      directly_lost <- igraph::V(lost_subgraph)[igraph::degree(lost_subgraph,
        mode =
          "in"
      ) == 0]$name
      directly_lost
    } else {
      directly_lost <- c()
    }
    directly_lost <- union(directly_lost, lost_cell_groups_not_in_graph)
  }
  # debug(get_dir_losses)

  lost_cell_groups <- lost_cell_groups %>%
    mutate(
      directly_lost_cell_groups = purrr::map(
        .f = get_dir_losses,
        .x = lost_cell_groups,
        state_graph
      ),
      indirectly_lost_cell_groups = purrr::map2(
        .f = function(x, y) {
          setdiff(x$cell_group, y)
        },
        .x = lost_cell_groups,
        .y = directly_lost_cell_groups
      )
    )
  # return (lost_cell_groups)

  direct_requirements <- lost_cell_groups %>%
    select(directly_lost_cell_groups, perturb_name) %>%
    tidyr::unnest(directly_lost_cell_groups) %>%
    dplyr::rename(id = directly_lost_cell_groups) %>%
    mutate(perturb_effect = "direct")
  indirect_requirements <- lost_cell_groups %>%
    select(indirectly_lost_cell_groups, perturb_name) %>%
    tidyr::unnest(indirectly_lost_cell_groups) %>%
    dplyr::rename(id = indirectly_lost_cell_groups) %>%
    mutate(perturb_effect = "indirect")

  requirements <- bind_rows(direct_requirements, indirect_requirements) %>% arrange(id)
  return(requirements)
  # node_direct_perturbs = dplyr::setdiff(node_direct_perturbs, node_indirect_perturbs)
}
# debug(categorize_genetic_requirements)


#' extract a single ccm from perturb_models_tbl
#' @param perturb_models_tbl tibble of perturbation results
#' @param perturb_name which perturbation to select out
get_perturb_ccm <- function(perturb_models_tbl, perturb_name) {
  perturb_ccm <- perturb_models_tbl %>%
    filter(perturb_name == perturb_name) %>%
    pull(perturb_ccm)
  perturb_ccm <- perturb_ccm[[1]]
  return(perturb_ccm)
}


#' Fit a Subset of Genotype Cell Cycle Models (CCM)
#'
#' This function fits a subset of genotype cell cycle models (CCM) based on the
#' provided parameters. It allows for switching to a specified UMAP space and
#' subsetting the cell cycle states (CCS) before fitting the genotype CCM.
#'
#' @param ccm A CCM object containing the cell cycle model and associated metadata.
#' @param umap_space A character string specifying the UMAP space to switch to.
#'   If `NULL`, the function attempts to retrieve the UMAP space from the CCM metadata.
#' @param ... Additional arguments passed to the `subset_ccs` function for subsetting CCS.
#'
#' @return A CCM object fitted to the subset of genotype data.
#'
#' @details
#' - If `umap_space` is not provided, the function tries to retrieve it from
#'   `ccm@ccs@cds@metadata$umap_space`.
#' - If a valid UMAP space is found, the function switches the CCM to that space
#'   using `switch_ccm_space`.
#' - The function subsets the CCS using `subset_ccs` and fits the genotype CCM
#'   using `fit_genotype_ccm`.
#'
#' @seealso
#' - `switch_ccm_space` for switching the CCM to a specified UMAP space.
#' - `subset_ccs` for subsetting cell cycle states.
#' - `fit_genotype_ccm` for fitting genotype cell cycle models.
#'
#' @examples
#' # Example usage:
#' # Assuming `ccm` is a valid CCM object:
#' sub_ccm <- fit_subset_genotype_ccm(ccm, umap_space = "UMAP_1", some_filter = TRUE)
#'
#' @export
fit_subset_genotype_ccm <- function(ccm, umap_space = NULL, ...) {
  # if i didn't specify a umap space, try to find one
  if (is.null(umap_space)) {
    umap_space <- ccm@ccs@cds@metadata$umap_space
  }
  # if it exists, switch
  if (is.null(umap_space) == FALSE) {
    ccm <- switch_ccm_space(ccm, umap_space = umap_space)
  }

  sub_ccs <- subset_ccs(ccm@ccs, ...)
  sub_ccm <- fit_genotype_ccm(ccm@info$genotype,
    sub_ccs,
    perturbation_col = ccm@info$perturbation_col,
    ctrl_ids = ccm@info$ctrl_ids
  )

  return(sub_ccm)
}





construct_platt_graph <- function(cds,
                                  partition_name,
                                  sample_group,
                                  comb_cds,
                                  perturbation_col = "perturbation",
                                  interval_col = "timepoint",
                                  component_col = "assembly_group",
                                  ctrl_ids = c(),
                                  num_threads = 6,
                                  k = 15,
                                  resolution = NULL,
                                  edge_allowlist = NULL,
                                  batch_col = "expt",
                                  newdata = tibble(),
                                  batches_excluded_from_assembly = c()) {
  cds <- cluster_cells(cds, random_seed = 42, res = resolution, k = k)
  colData(cds)$cell_state <- monocle3:::clusters(cds)

  res_cluster <- assemble_partition(cds,
    partition_name = partition_name,
    sample_group = sample_group,
    cell_group = cell_group,
    interval_col = interval_col,
    component_col = component_col,
    perturbation_col = perturbation_col,
    ctrl_ids = ctrl_ids,
    num_threads = num_threads,
    batch_col = batch_col,
    newdata = newdata,
    batches_excluded_from_assembly = batches_excluded_from_assembly
  )

  contract_graph <- platt::contract_state_graph(ccs,
    mt_graph = res_cluster$mt_graph[[1]],
    group_nodes_by = "cell_type"
  )

  global_wt_graph_edge_allowlist <- igraph::as_data_frame(contract_graph)
  global_wt_graph_edge_allowlist <- global_wt_graph_edge_allowlist %>%
    select(from, to) %>%
    distinct()

  res_cell_type <- assemble_partition(cds,
    partition_name = partition_name,
    sample_group = sample_group,
    cell_group = cell_group,
    interval_col = interval_col,
    component_col = component_col,
    perturbation_col = perturbation_col,
    ctrl_ids = ctrl_ids,
    num_threads = num_threads,
    batch_col = batch_col,
    newdata = newdata,
    batches_excluded_from_assembly = batches_excluded_from_assembly,
    edge_allowlist = global_wt_graph_edge_allowlist
  )



  return(res_cell_type)
}
