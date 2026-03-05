#' Get perturbation time window from a CCS object
#'
#' Returns the start and stop times (based on `interval_col`) for samples
#' matching a perturbation or genotype within a cell count dataset.
#'
#' @param genotype Character or character vector. Perturbation/genotype IDs to match.
#' @param ccs A cell count set object with sample metadata in `colData(ccs)`.
#' @param interval_col Character. Column name in `colData(ccs)` with time values.
#' @param perturbation_col Character. Column name in `colData(ccs)` with
#'   perturbation/genotype labels. Default is "perturbation".
#'
#' @return A tibble with `start_time` and `stop_time` columns.
#'
#' @export
get_time_window <- function(genotype, ccs, interval_col, perturbation_col = "perturbation") {
  subset_ccs <- ccs[, replace_na(colData(ccs)[[perturbation_col]] %in% genotype, F)]
  colData(subset_ccs)$knockout <- colData(subset_ccs)[[perturbation_col]] %in% genotype
  knockout_time_start <- min(colData(subset_ccs)[[interval_col]][colData(subset_ccs)$knockout])
  knockout_time_stop <- max(colData(subset_ccs)[[interval_col]][colData(subset_ccs)$knockout])
  return(tibble(start_time = knockout_time_start, stop_time = knockout_time_stop))
}

#' Summarize perturbation effects across timepoints
#'
#' Builds a tidy table of perturbation effects by combining per-timepoint
#' contrasts from a fitted CCM.
#'
#' @param ccm A fitted CCM object with `ccs` and model information.
#' @param interval_col Character. Column name in `colData(ccm@ccs)` that defines
#'   timepoints. Default is "timepoint".
#' @param newdata Tibble. Optional metadata to cross-join with timepoints before
#'   computing contrasts.
#' @param adjust_q_values Logical. Whether to adjust q-values in the contrast
#'   calculation. Default is FALSE.
#'
#' @return A tibble with one row per timepoint (and any `newdata` combinations)
#'   plus contrast results.
#'
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

#' Build WT vs KO contrasts from a fitted CCM
#'
#' Convenience wrapper around `estimate_abundances()` and `compare_abundances()`
#' to compute WT vs KO contrasts for a fitted CCM, optionally across additional
#' covariates supplied in `newdata`.
#'
#' @param ccm A fitted cell count model object.
#' @param newdata Tibble. Optional covariates to cross-join with knockout status
#'   before estimating abundances.
#' @param adjust_q_values Logical. Whether to adjust q-values in
#'   `compare_abundances()`. Default is FALSE.
#'
#' @return A tibble of contrast results as returned by `compare_abundances()`.
#'
#' @details
#' When `newdata` is provided, this function evaluates contrasts for each
#' combination of covariates and knockout status by cross-joining with a
#' `knockout` column set to FALSE (WT) or TRUE (KO).
#'
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
                         min_penalty = 0.01,
                         max_penalty = 1e+06,
                         num_bootstraps = 10,
                         ...) {
  if (is.null(ctrl_ids)) {
    ctrl_ids <- unique(colData(cds)[[perturbation_col]])
    ctrl_ids <- ctrl_ids[grepl("wt|ctrl|reference", ctrl_ids)]
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
    min_penalty = min_penalty,
    max_penalty = max_penalty,
    num_bootstraps = num_bootstraps,
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
                              links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                              log_abund_detection_thresh = -5,
                              q_val = 0.1,
                              min_interval = 4,
                              max_interval = 24,
                              min_pathfinding_lfc = 0,
                              break_cycles = TRUE,
                              edge_allowlist = NULL,
                              edge_denylist = NULL,
                              component_col = "partition",
                              force_allowlist = FALSE,
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
    newdata = newdata,
    force_allowlist = force_allowlist
  )
  if (break_cycles) {
    print("breaking cycles in control timeseries graph...")
    wt_state_transition_graph <- platt:::break_cycles_in_state_transition_graph(wt_state_transition_graph, "support")
  }

  # ensure edges again
  if (!is.null(edge_allowlist)) {
    # Add missing edges from edge_allowlist
    current_edges <- igraph::as_data_frame(wt_state_transition_graph, what = "edges") %>% select(from, to)
    missing_edges <- anti_join(edge_allowlist, current_edges, by = c("from", "to"))
    if (nrow(missing_edges) > 0) {
      # Add missing edges with default support value
      for (i in seq_len(nrow(missing_edges))) {
        wt_state_transition_graph <- igraph::add_edges(
          wt_state_transition_graph,
          c(missing_edges$from[i], missing_edges$to[i]),
          attr = list(support = NA)
        )
      }
    }
  }

  if (!is.null(edge_denylist)) {
    # Remove edges in edge_denylist
    edges_to_remove <- igraph::as_data_frame(wt_state_transition_graph, what = "edges") %>%
      inner_join(edge_denylist, by = c("from", "to"))
    if (nrow(edges_to_remove) > 0) {
      vp <- as.vector(t(as.matrix(edges_to_remove[, c("from", "to")])))
      edge_ids <- igraph::get_edge_ids(wt_state_transition_graph, vp)
      wt_state_transition_graph <- igraph::delete_edges(wt_state_transition_graph, edge_ids)
    }
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

#' Assemble mutant transition graph from perturbation models
#'
#' Constructs a mutant supergraph by combining a WT graph with perturbation
#' models, with optional cycle breaking.
#'
#' @param ref_ccs Reference CCS object used for WT context.
#' @param wt_graph WT transition graph (igraph).
#' @param perturb_models_tbl Tibble of perturbation models, typically from
#'   `fit_mt_models()` / `assess_perturbation_effects()`.
#' @param interval_col Character. Column name in `colData(ref_ccs@cds)` with
#'   time values. Default is "timepoint".
#' @param start_time Numeric. Start time for graph assembly; defaults to min
#'   observed time.
#' @param stop_time Numeric. Stop time for graph assembly; defaults to max
#'   observed time.
#' @param interval_step Numeric. Step size for time discretization. Default is 2.
#' @param links_between_components Character. Strategy for linking components.
#' @param log_abund_detection_thresh Numeric. Log abundance threshold.
#' @param q_val Numeric. Q-value threshold for perturbation effects.
#' @param newdata Tibble. Optional covariates for effect estimation.
#' @param break_cycles Logical. Whether to remove cycles in the assembled graph.
#' @param component_col Character. Column name defining components.
#' @param edge_allowlist Optional allowlist for edges.
#' @param edge_denylist Optional denylist for edges.
#' @param verbose Logical. Emit progress messages.
#'
#' @return An igraph mutant supergraph.
#'
#' @export
assemble_mt_graph <- function(ref_ccs,
                              wt_graph,
                              perturb_models_tbl,
                              interval_col = "timepoint",
                              # perturbation_col = "knockout",
                              start_time = NULL,
                              stop_time = NULL,
                              interval_step = 2,
                              links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                              log_abund_detection_thresh = -5,
                              q_val = 0.1,
                              newdata = tibble(),
                              break_cycles = TRUE,
                              component_col = "partition",
                              edge_allowlist = NULL,
                              edge_denylist = NULL,
                              discordant_pruning_mode = c("none", "greedy"),
                              discordant_pruning_k = 1,
                              discordant_pruning_power_threshold = 0,
                              discordant_pruning_cost_attr = "total_path_score_supporting",
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
    discordant_pruning_mode = discordant_pruning_mode,
    discordant_pruning_k = discordant_pruning_k,
    discordant_pruning_power_threshold = discordant_pruning_power_threshold,
    discordant_pruning_cost_attr = discordant_pruning_cost_attr,
    components = component_col,
    verbose = verbose
  )
  if (break_cycles & inherits(mutant_supergraph, "igraph")) {
    print("breaking cycles in perturbation graph...")
    mutant_supergraph <- platt:::break_cycles_in_state_transition_graph(mutant_supergraph, "total_perturb_path_score_supporting")
  }

  return(mutant_supergraph)
}


#' Categorize genetic requirements by cell group
#'
#' Summarizes perturbation effects to identify cell groups that are directly or
#' indirectly lost for each perturbation, using a state transition graph.
#'
#' @param perturb_ccm_tbl Tibble of perturbation models with a
#'   `perturb_summary_tbl` column.
#' @param state_graph State transition graph (igraph) defining relationships
#'   among cell groups.
#'
#' @return A tibble with per-perturbation lists of directly and indirectly lost
#'   cell groups.
#'
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


#' Return a dataframe that describes which cell types are present in a time interval
#'
#' @export
get_extant_cell_types <- function(ccm,
                                  start,
                                  stop,
                                  interval_col = "timepoint",
                                  interval_step = 2,
                                  log_abund_detection_thresh = 0,
                                  pct_dynamic_range = 0.25,
                                  pct_range_detection_thresh = pct_dynamic_range,
                                  min_cell_range = 2,
                                  newdata = tibble()) {
  timepoint_pred_df <- estimate_abundances_over_interval(ccm, start, stop,
    interval_col = interval_col, interval_step = interval_step, newdata = newdata
  )

  norm_mat <- normalized_counts(ccm@ccs, "size_only")
  norm_mat[norm_mat == 0] <- NA
  count_quantiles <- sparseMatrixStats::rowQuantiles(norm_mat, probs = seq(from = 0, to = 1, by = pct_range_detection_thresh), na.rm = T)
  count_ranges <- sparseMatrixStats::rowRanges(norm_mat, na.rm = T)
  row.names(count_ranges) <- row.names(count_quantiles)

  cell_type_thresh_df <- tibble(
    cell_group = timepoint_pred_df %>% pull(cell_group) %>% unique(),
    cell_group_pct_range_detection_thresh = count_quantiles[cell_group, 2],
    min_count = count_ranges[cell_group, 1],
    max_count = count_ranges[cell_group, 2]
  )

  timepoint_pred_df <- timepoint_pred_df %>% left_join(cell_type_thresh_df)

  if (is.null(log_abund_detection_thresh)) {
    log_abund_detection_thresh <- abund_range[1] + pct_dynamic_range * dynamic_range
  }

  timepoint_pred_df <- timepoint_pred_df %>%
    group_by(cell_group) %>%
    mutate(
      max_abundance = max(exp(log_abund)),
      percent_max_abund = exp(log_abund) / max_abundance,
      cell_type_prediction_range = max(log_abund) - (min(log_abund)),
      percent_cell_type_range = (log_abund - min(log_abund)) / cell_type_prediction_range,
      # above_log_abund_thresh = (log_abund - 2*log_abund_se > log_abund_detection_thresh & log_abund - 2*log_abund_se > log(cell_group_pct_range_detection_thresh)) | cell_type_prediction_range < min_cell_range,
      above_log_abund_thresh = log_abund > log_abund_detection_thresh,
      present_flag = ifelse(above_log_abund_thresh, TRUE, NA)
    ) %>%
    ungroup()

  longest_present_interval <- function(tps_df) {
    tryCatch(
      {
        delta_t <- as.numeric(tps_df[2, 1] - tps_df[1, 1])
        ts_la <- ts(tps_df$present_flag,
          start = min(tps_df[, 1]),
          # end=max(tps_df[,1]),
          deltat = delta_t
        )
        longest_contig <- na.contiguous(ts_la)
        return(tibble(longest_contig_start = start(longest_contig)[1], longest_contig_end = end(longest_contig)[1]))
      },
      error = function(e) {
        return(tibble(longest_contig_start = NA, longest_contig_end = NA))
      }
    )
  }

  # undebug(longest_present_interval)
  nested_timepoints_df <- timepoint_pred_df %>%
    select(cell_group, !!sym(interval_col), present_flag) %>%
    group_by(cell_group)

  nested_timepoints_df <- nested_timepoints_df %>%
    group_modify(~ longest_present_interval(.x))
  # mutate(cg_ts = purrr:::map2(.f=purrr::possibly(longest_present_interval, NA_real_),
  #                            .x=!!sym(interval_col),
  #                            .y=present_flag))

  timepoint_pred_df <- left_join(timepoint_pred_df, nested_timepoints_df)

  timepoint_pred_df <- timepoint_pred_df %>%
    mutate(
      present_above_thresh = !!sym(interval_col) >= longest_contig_start & !!sym(interval_col) <= longest_contig_end,
      present_above_thresh = ifelse(is.na(present_above_thresh), FALSE, present_above_thresh)
    )

  extant_cell_type_df <- timepoint_pred_df %>%
    select(
      !!sym(interval_col),
      cell_group,
      log_abund,
      max_abundance,
      percent_max_abund,
      percent_cell_type_range,
      longest_contig_start,
      longest_contig_end,
      present_above_thresh
    )
  return(extant_cell_type_df)
}
# undebug(get_extant_cell_types)
# get_extant_cell_types(wt_ccm_wl, 72, 96) %>% filter(cell_group == "21")
