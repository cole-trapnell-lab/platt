# FIXME: we need to standardize notation (use either "partition" or "component") throughout the code
#' @noRd
add_cross_component_pathfinding_links <- function(ccm,
                                                  pathfinding_graph,
                                                  extant_cell_type_df,
                                                  type = c("strongest-pcor", "strong-pcor", "ctp"),
                                                  surprise_thresh = 2,
                                                  components = "partition") {
  assertthat::assert_that(
    components == "partition" ||
      tryCatch(
        expr = components %in% colnames(colData(ccm@ccs@cds)),
        error = function(e) FALSE
      ),
    msg = paste0(components, "not found in colData")
  )

  pcor_graph <- tibble(cell_group = row.names(counts(ccm@ccs)))
  pcor_graph <- pcor_graph %>%
    select(cell_group) %>%
    tidyr::expand(cell_group, cell_group)
  colnames(pcor_graph) <- c("from", "to")
  pcor_graph <- pcor_graph %>% filter(from != to)

  nz_pcors <- igraph::graph_from_adjacency_matrix(model(ccm, "reduced")$latent_network(type = "partial_cor"), mode = "undirected", weighted = TRUE, diag = FALSE) %>% igraph::as.directed()
  nz_pcor_graph_edges <- igraph::as_data_frame(nz_pcors, what = "edges") %>%
    dplyr::rename(pcor = weight) %>%
    select(from, to, pcor) # %>% dplyr::filter(pcor != 0.00)

  pcor_graph <- dplyr::left_join(pcor_graph, nz_pcor_graph_edges)

  if (components == "partition") {
    component_assignments <- partitions(ccm@ccs@cds)
  } else {
    component_assignments <- colData(ccm@ccs@cds)[, components]
    names(component_assignments) <- row.names(colData(ccm@ccs@cds))
  }

  clusters_by_partition <- tibble(
    cell_group = ccm@ccs@metadata$cell_group_assignments$cell_group,
    partition = component_assignments[row.names(ccm@ccs@metadata$cell_group_assignments)]
  )
  clusters_by_partition <- clusters_by_partition %>%
    group_by(cell_group, partition) %>%
    summarize(cells_from_group_in_partition = n())
  clusters_by_partition <- clusters_by_partition %>%
    group_by(cell_group) %>%
    slice_max(cells_from_group_in_partition, n = 1)
  clusters_by_partition <- clusters_by_partition %>% select(cell_group, partition)

  # FIXME: we should be indexing extant_cell_type_df by name, but this function doesn't know what the
  # timepoint column is called yet
  start_time <- min(extant_cell_type_df[, 1])
  clusters_by_partition <- clusters_by_partition %>% left_join(extant_cell_type_df, by = "cell_group")
  partitions_present_at_start <- clusters_by_partition %>%
    filter(timepoint == start_time & present_above_thresh) %>%
    pull(partition) %>%
    unique()
  partitions_absent_at_start <- setdiff(clusters_by_partition$partition, partitions_present_at_start)
  cell_groups_in_emergent_partitions <- clusters_by_partition %>%
    filter(partition %in% partitions_absent_at_start) %>%
    pull(cell_group) %>%
    unique()

  # only add cross-partition links between cell groups that aren't present
  # at the start of the assembly time window
  cross_partition_map <- clusters_by_partition %>%
    ungroup() %>%
    select(cell_group) %>%
    tidyr::expand(cell_group, cell_group)
  colnames(cross_partition_map) <- c("from", "to")
  cross_partition_map <- cross_partition_map %>% filter(from != to)
  # Only allow cross partition links when one cell group is in an emergent partition
  cross_partition_map <- cross_partition_map %>% filter(from %in% cell_groups_in_emergent_partitions |
    to %in% cell_groups_in_emergent_partitions)
  cross_partition_map <- dplyr::left_join(cross_partition_map, clusters_by_partition %>% setNames(paste0("to_", names(.))), by = c("to" = "to_cell_group"), relationship = "many-to-many") # %>%
  cross_partition_map <- dplyr::left_join(cross_partition_map, clusters_by_partition %>% setNames(paste0("from_", names(.))), by = c("from" = "from_cell_group"), relationship = "many-to-many") # %>%
  cross_partition_map <- cross_partition_map %>% filter(from_partition != to_partition)

  # only add cross-partition links between cell groups that are present
  # at the same time
  cross_partition_map <- cross_partition_map %>% filter(from_timepoint == to_timepoint &
    from_present_above_thresh &
    to_present_above_thresh)

  # pcor_graph = pcor_graph %>% tidyr::replace_na(list(pcor = 0))

  # pcor_graph = pcor_graph %>% left_join()
  pcor_graph <- hooke:::weigh_edges_by_umap_dist(ccm, pcor_graph) %>% dplyr::rename(umap_dist = weight)

  # pcor_graph = pcor_graph %>% filter(pcor != 0)
  pcor_graph <- pcor_graph %>% mutate(abs_pcor = abs(pcor))

  # pcor_vs_dist_model = VGAM::vglm(pcor ~ I(1/umap_dist), data=pcor_graph, family=VGAM::gamma2(zero=NULL, lmu="identitylink", lshape="loglink"))
  # mod_predict = predict(pcor_vs_dist_model, newdata=pcor_graph) %>% as.data.frame()
  # pcor_graph$model_fit = mod_predict$mu
  # pcor_graph$model_sd = sqrt((mod_predict$mu)^2 / exp(mod_predict$`loglink(shape)`))

  # pcor_vs_dist_model = VGAM::vglm(pcor ~ umap_dist, data=pcor_graph, family=VGAM::gaussianff(zero=NULL, lmean="identitylink", lsd="loglink"))
  pcor_vs_dist_model <- VGAM::vglm(pcor ~ umap_dist, data = pcor_graph, family = VGAM::uninormal(zero = NULL, lmean = "identitylink", lsd = "loglink"))
  mod_predict <- predict(pcor_vs_dist_model, newdata = pcor_graph) %>% as.data.frame()
  pcor_graph$model_fit <- mod_predict$mean
  # pcor_graph$model_sd = sqrt((mod_predict$mean)^2 / exp(mod_predict$`loglink(sd)`))
  # pcor_graph$model_sd = sqrt((mod_predict$mean)^2 / mod_predict$sd)
  pcor_graph$model_sd <- mod_predict$sd
  pcor_graph$model_sd <- exp(mod_predict$`loglink(sd)`)


  pcor_graph <- pcor_graph %>%
    mutate(
      model_fit_lower = model_fit - surprise_thresh * model_sd,
      # model_fit_lower = ifelse(model_fit_lower > 0, model_fit_lower, min(abs_pcor)),
      model_fit_upper = model_fit + surprise_thresh * model_sd
    )

  cross_partition_map <- cross_partition_map %>% left_join(pcor_graph, by = c("from" = "from", "to" = "to"))



  if (type == "strongest-pcor") {
    cross_partition_edges <- cross_partition_map %>%
      group_by(to_partition) %>%
      group_by(from_partition, to_partition) %>%
      slice_max(abs_pcor, n = 1) %>%
      ungroup() %>%
      # filter (pcor < model_fit_lower | pcor > model_fit_upper) %>%
      select(from, to, weight = umap_dist)
  } else if (type == "strong-pcor") {
    cross_partition_edges <- cross_partition_map %>%
      filter(pcor < model_fit_lower | pcor > model_fit_upper) %>%
      select(from, to, weight = umap_dist)
  } else if (type == "ctp") {
    cross_partition_edges <- cross_partition_map %>%
      group_by(to_partition) %>%
      slice_min(to_timepoint) %>%
      slice_min(umap_dist) %>%
      ungroup() %>%
      select(from, to, weight = umap_dist)
  }

  cross_partition_graph <- igraph::graph_from_data_frame(cross_partition_edges, directed = FALSE, vertices = data.frame(id = row.names(counts(ccm@ccs)))) %>% igraph::as.directed()


  # unexpectedly_strong_pcor_edges = pcor_graph %>%
  #  filter (pcor < model_fit_lower | pcor > model_fit_upper) %>%
  #  select(from, to, weight=umap_dist)
  # uspe_graph = igraph::graph_from_data_frame(unexpectedly_strong_pcor_edges, directed=TRUE, vertices=data.frame(id=row.names(counts(ccm@ccs))))
  updated_pathfinding_graph <- igraph::union(pathfinding_graph, cross_partition_graph)
  updated_pathfinding_graph <- igraph::simplify(updated_pathfinding_graph)
  return(updated_pathfinding_graph)
}



#' Initialize a graph over which cells can transition
#'
#' @param ccm A cell count model
#' @param extant_cell_type_df A data frame describing when cell types are present, generated by get_extant_cell_types()
#' @param allow_links_between_components Whether cells in separate partitions of the UMAP can be linked by paths
#' @export
init_pathfinding_graph <- function(ccm,
                                   extant_cell_type_df,
                                   links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                                   components = "partition",
                                   weigh_by_pcor = F,
                                   edge_allowlist = NULL,
                                   edge_denylist = NULL,
                                   force_allowlist = FALSE) {
  # There are a number of different ways we could set up this "pathfinding graph" but for now
  # Let's just use the PAGA (weighed by distance in UMAP space), subtracting edges between which
  # there is zero partial correlation in the nuisance cell count model.

  links_between_components <- match.arg(links_between_components)

  cell_groups <- ccm@ccs@metadata[["cell_group_assignments"]] %>%
    pull(cell_group) %>%
    unique()
  node_metadata <- data.frame(id = cell_groups)


  if (force_allowlist == FALSE) {
    message("Initializing pathfinding graph from partially correlated pairs linked in PAGA")
    paga_graph <- initial_pcor_graph(ccm@ccs) %>%
      igraph::graph_from_data_frame(directed = FALSE, vertices = node_metadata) %>%
      igraph::as.directed()

    cov_graph <- hooke:::return_igraph(model(ccm, "reduced"))
    cov_graph_edges <- igraph::as_data_frame(cov_graph, what = "edges")

    if (nrow(cov_graph_edges) > 0 & "weight" %in% colnames(cov_graph_edges)) {
      cov_graph_edges <- cov_graph_edges %>%
        dplyr::rename(pcor = weight) %>%
        dplyr::filter(pcor != 0.00)
    }

    cov_graph_edges$to <- as.character(cov_graph_edges$to)
    cov_graph_edges$from <- as.character(cov_graph_edges$from)

    weighted_edges <- hooke:::weigh_edges_by_umap_dist(ccm, cov_graph_edges)

    if (is.null(edge_allowlist) == FALSE) {
      weighted_edges_allow <- hooke:::weigh_edges_by_umap_dist(ccm, edge_allowlist) %>% mutate(pcor = 1)
      weighted_edges <- rbind(weighted_edges, weighted_edges_allow) %>%
        select(from, to, weight) %>%
        distinct()
    }

    paga_components <- igraph::components(paga_graph)
    same_partition_mat <- outer(paga_components$membership, paga_components$membership, FUN = "==")
    weighted_edges <- weighted_edges %>%
      group_by(from, to) %>%
      mutate(
        adjacent_in_paga = igraph::are_adjacent(paga_graph, from, to),
        same_partition = same_partition_mat[from, to]
      ) %>%
      ungroup()

    if (is.null(edge_allowlist) == FALSE) {
      # if in the allowlist, we want to keep the edge even if it's not adjacent in paga
      edge_allowlist$in_allowlist <- TRUE
      weighted_edges <- left_join(weighted_edges, edge_allowlist, by = c("from", "to"))
      weighted_edges$in_allowlist <- ifelse(is.na(weighted_edges$in_allowlist), F, weighted_edges$in_allowlist)

      weighted_edges <- weighted_edges %>% dplyr::filter(adjacent_in_paga | in_allowlist)
    } else {
      weighted_edges <- weighted_edges %>% dplyr::filter(adjacent_in_paga)
    }

    pathfinding_graph <- weighted_edges %>%
      transmute(
        from = pmin(from, to),
        to   = pmax(from, to),
        weight
      ) %>%
      group_by(from, to) %>%
      summarise(weight = mean(weight), .groups = "drop") %>% # or first/max/min
      filter(from != to) %>%
      igraph::graph_from_data_frame(directed = FALSE, vertices = node_metadata) %>%
      igraph::as.directed()
  } else {
    weighted_edges <- hooke:::weigh_edges_by_umap_dist(ccm, edge_allowlist)
    pathfinding_graph <- weighted_edges %>%
      select(from, to, weight) %>%
      igraph::graph_from_data_frame(directed = TRUE, vertices = node_metadata)
  }

  if (links_between_components != "none") {
    pathfinding_graph <- add_cross_component_pathfinding_links(ccm,
      pathfinding_graph,
      extant_cell_type_df,
      type = links_between_components,
      components = components
    )
    pathfinding_graph <- igraph::as_data_frame(pathfinding_graph) %>% select(from, to)
    pathfinding_graph <- hooke:::weigh_edges_by_umap_dist(ccm, pathfinding_graph) %>%
      igraph::graph_from_data_frame(directed = FALSE, vertices = node_metadata) %>%
      igraph::as.directed() %>%
      igraph::simplify()
  }

  if (is.null(edge_denylist) == FALSE) {
    denylist_graph <- igraph::graph_from_data_frame(edge_denylist, directed = TRUE, vertices = node_metadata)
    pathfinding_graph <- pathfinding_graph - denylist_graph
  }

  if (is.null(edge_allowlist) == FALSE) {
    # subtract the opposite approvelist
    edge_allowlist_opp <- edge_allowlist %>% select(from = to, to = from)
    edge_allowlist_opp_graph <- igraph::graph_from_data_frame(edge_allowlist_opp, directed = TRUE, vertices = node_metadata)
    pathfinding_graph <- pathfinding_graph - edge_allowlist_opp_graph
  }


  # if (is.null(edge_allowlist) == FALSE){
  #   # Ensuring allowlisted edges remain in pathfinding graph
  #   edge_allowlist = edge_allowlist[,c(1,2)] %>% as_tibble()
  #   colnames(edge_allowlist) = c("from", "to")
  #   pathfinding_graph = igraph::as_data_frame(pathfinding_graph) %>%  select(from, to)
  #   pathfinding_graph = pathfinding_graph %>% bind_rows(edge_allowlist) %>% distinct()
  #   pathfinding_graph = hooke:::weigh_edges_by_umap_dist(ccm, pathfinding_graph)
  # }


  return(pathfinding_graph)
}



#' Generate a denylist of state transition relationships based on a perturbation
#' experiment.
#'
#' This function generates a denylist of edges between cell states based on
#' the idea that if cell state X is lost, and Y is not ever lost in the experiment
#' Y can't come from X. This function is very simplistic right now. We could get
#' more sophisticated by looking at the relative timing of losses, etc. We are
#' also not accounting for power at all. We should be asking whether we have
#' power to detect a change in state Y, and if not, exclude (X,Y) from the
#' denylist
#'
#' @noRd
get_discordant_loss_pairs <- function(perturbation_ccm,
                                      perturb_time_window,
                                      control_timeseries_ccm,
                                      control_time_window,
                                      interval_step,
                                      interval_col,
                                      log_abund_detection_thresh,
                                      q_val,
                                      model_for_pcors = "reduced",
                                      min_pathfinding_lfc = 0,
                                      power_threshold = 0,
                                      newdata = tibble()) {
  # print ("getting perturbation paths")
  # print (time_window)
  if (is.null(perturbation_ccm) || is.na(perturbation_ccm)) {
    message("Error: no perturbation model")
    return(NA)
  }

  perturb_start_time <- min(as.numeric(perturb_time_window$start_time))
  perturb_stop_time <- min(as.numeric(perturb_time_window$stop_time))

  control_start_time <- min(as.numeric(control_time_window$start_time))
  control_stop_time <- min(as.numeric(control_time_window$stop_time))

  # If the perturbation experiment model covers a winder time interval than the
  # control timeseries, just limit testing to the perturbation experiment's
  # time interval
  if (perturb_start_time < control_start_time) {
    perturb_start_time <- control_start_time
  }

  if (perturb_stop_time < control_stop_time) {
    perturb_stop_time <- control_stop_time
  }

  if (nrow(newdata) > 0) {
    newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
  } else {
    newdata_wt <- tibble(knockout = FALSE)
  }

  wt_timepoint_pred_df <- estimate_abundances_over_interval(control_timeseries_ccm, control_start_time, control_stop_time,
    interval_col = interval_col, interval_step = interval_step, newdata = newdata_wt
  )
  peak_wt_abundance <- wt_timepoint_pred_df %>%
    group_by(cell_group) %>%
    slice_max(log_abund, n = 1)
  peak_outside_perturbation_window <- peak_wt_abundance %>%
    filter(!!sym(interval_col) > perturb_stop_time | !!sym(interval_col) < perturb_start_time) %>%
    pull(cell_group)

  message("\tEstimating loss timing")
  earliest_loss_tbl <- estimate_loss_timing(perturbation_ccm,
    start_time = perturb_start_time,
    stop_time = perturb_stop_time,
    interval_step = interval_step,
    interval_col = interval_col,
    control_ccm = control_timeseries_ccm,
    control_start_time = control_start_time,
    control_stop_time = control_stop_time,
    log_abund_detection_thresh = log_abund_detection_thresh,
    q_val = q_val,
    delta_log_abund_loss_thresh = min_pathfinding_lfc,
    newdata = newdata
  )

  assert_power_columns_present(earliest_loss_tbl, context = "discordant-loss pair construction")

  discordant_loss_pairs <- compute_discordant_loss_pairs(
    earliest_loss_tbl = earliest_loss_tbl,
    peak_outside_perturbation_window = peak_outside_perturbation_window,
    power_threshold = power_threshold
  )

  return(discordant_loss_pairs)
}

#' @noRd
assert_power_columns_present <- function(tbl, context = "perturbation summary") {
  power_cols <- colnames(tbl)[stringr::str_detect(colnames(tbl), "power")]
  if (length(power_cols) == 0) {
    stop(sprintf("Expected at least one power column in %s, but none were found.", context))
  }
  invisible(power_cols)
}

#' @noRd
compute_discordant_loss_pairs <- function(earliest_loss_tbl,
                                          peak_outside_perturbation_window = character(),
                                          power_threshold = 0) {
  lost_cell_groups <- earliest_loss_tbl %>%
    dplyr::filter(is_lost_at_peak) %>%
    dplyr::pull(cell_group) %>%
    unique()

  unaffected_cell_groups <- earliest_loss_tbl %>%
    dplyr::filter(peak_time_in_ctrl_within_perturb_time_range) %>%
    dplyr::pull(cell_group) %>%
    unique()
  unaffected_cell_groups <- setdiff(unaffected_cell_groups, lost_cell_groups)

  # Exclude states that peak outside the window of this perturbation experiment,
  # as we may simply have not yet seen their loss.
  unaffected_cell_groups <- setdiff(unaffected_cell_groups, peak_outside_perturbation_window)

  if (power_threshold > 0) {
    power_cols <- assert_power_columns_present(earliest_loss_tbl, context = "discordant-loss pair construction")

    cell_group_power_tbl <- earliest_loss_tbl %>%
      dplyr::select(cell_group, dplyr::all_of(power_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(power_cols),
        names_to = "power_metric",
        values_to = "power_value"
      ) %>%
      dplyr::group_by(cell_group) %>%
      dplyr::summarize(
        max_power = ifelse(all(is.na(power_value)), NA_real_, max(power_value, na.rm = TRUE)),
        .groups = "drop"
      )

    powered_unaffected_cell_groups <- cell_group_power_tbl %>%
      dplyr::filter(!is.na(max_power), max_power >= power_threshold) %>%
      dplyr::pull(cell_group)

    unaffected_cell_groups <- intersect(unaffected_cell_groups, powered_unaffected_cell_groups)
  }

  tidyr::expand_grid(lost_cell_groups = lost_cell_groups, unaffected_cell_groups = unaffected_cell_groups)
}

#' Compute edge discordance penalties from perturbation summaries
#'
#' Aggregates discordant evidence for each directed edge in a state graph across
#' perturbations, weighted by downstream-state power.
#'
#' @param perturbation_ccm_tbl Tibble with `perturb_name` and nested
#'   `perturb_summary_tbl` columns.
#' @param state_transition_graph Directed igraph whose edges are scored.
#' @param power_threshold Numeric. Minimum downstream power required for
#'   unaffected states to contribute discordant evidence.
#'
#' @return Tibble with one row per edge, including `discordance_penalty` and
#'   summary columns; existing edge-support columns can be preserved by joining
#'   this table on `from`/`to`.
#' @noRd
compute_edge_discordance_penalty <- function(perturbation_ccm_tbl,
                                             state_transition_graph,
                                             power_threshold = 0) {
  edge_tbl <- igraph::as_data_frame(state_transition_graph, what = "edges") %>%
    dplyr::select(from, to) %>%
    dplyr::distinct()

  if (nrow(edge_tbl) == 0) {
    return(tibble::tibble(
      from = character(),
      to = character(),
      discordance_penalty = numeric(),
      num_discordant_perturbs = integer(),
      mean_discordant_power = numeric()
    ))
  }

  perturb_summaries <- perturbation_ccm_tbl %>%
    dplyr::filter(!is.na(perturb_summary_tbl)) %>%
    dplyr::select(perturb_name, perturb_summary_tbl) %>%
    tidyr::unnest(cols = c(perturb_summary_tbl))

  if (nrow(perturb_summaries) == 0) {
    return(edge_tbl %>%
      dplyr::mutate(
        discordance_penalty = 0,
        num_discordant_perturbs = 0L,
        mean_discordant_power = 0
      ))
  }

  power_cols <- colnames(perturb_summaries)[stringr::str_detect(colnames(perturb_summaries), "power")]
  if (length(power_cols) > 0) {
    perturb_cell_power <- perturb_summaries %>%
      dplyr::select(perturb_name, cell_group, dplyr::all_of(power_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(power_cols),
        names_to = "power_metric",
        values_to = "power_value"
      ) %>%
      dplyr::group_by(perturb_name, cell_group) %>%
      dplyr::summarize(
        max_power = ifelse(all(is.na(power_value)), NA_real_, max(power_value, na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    perturb_cell_power <- perturb_summaries %>%
      dplyr::select(perturb_name, cell_group) %>%
      dplyr::distinct() %>%
      dplyr::mutate(max_power = 1)
  }

  lost_states <- perturb_summaries %>%
    dplyr::filter(is_lost_when_present) %>%
    dplyr::transmute(
      perturb_name,
      from = cell_group,
      loss_effect = pmax(-loss_when_present, 0)
    )

  unaffected_powered_states <- perturb_summaries %>%
    dplyr::filter(!is_lost_when_present) %>%
    dplyr::select(perturb_name, to = cell_group) %>%
    dplyr::left_join(perturb_cell_power, by = c("perturb_name", "to" = "cell_group")) %>%
    dplyr::filter(!is.na(max_power), max_power >= power_threshold)

  edge_discordant_evidence <- edge_tbl %>%
    dplyr::inner_join(lost_states, by = "from", relationship = "many-to-many") %>%
    dplyr::inner_join(unaffected_powered_states, by = c("perturb_name", "to"), relationship = "many-to-many")

  if (nrow(edge_discordant_evidence) == 0) {
    return(edge_tbl %>%
      dplyr::mutate(
        discordance_penalty = 0,
        num_discordant_perturbs = 0L,
        mean_discordant_power = 0
      ))
  }

  # Formula: discordance_penalty(edge) = sum_over_perturbations(max_power_to * max(0, -loss_when_present_from))
  discordance_tbl <- edge_discordant_evidence %>%
    dplyr::mutate(discordant_weight = max_power * loss_effect) %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarize(
      discordance_penalty = sum(discordant_weight, na.rm = TRUE),
      num_discordant_perturbs = dplyr::n_distinct(perturb_name),
      mean_discordant_power = mean(max_power, na.rm = TRUE),
      .groups = "drop"
    )

  edge_tbl %>%
    dplyr::left_join(discordance_tbl, by = c("from", "to")) %>%
    dplyr::mutate(
      discordance_penalty = ifelse(is.na(discordance_penalty), 0, discordance_penalty),
      num_discordant_perturbs = ifelse(is.na(num_discordant_perturbs), 0L, as.integer(num_discordant_perturbs)),
      mean_discordant_power = ifelse(is.na(mean_discordant_power), 0, mean_discordant_power)
    )
}

#' Build discordant source-target pairs from perturbation summaries
#' @noRd
collect_discordant_pairs_from_perturbation_summaries <- function(perturbation_ccm_tbl,
                                                                 power_threshold = 0) {
  perturb_summaries <- perturbation_ccm_tbl %>%
    dplyr::filter(!is.na(perturb_summary_tbl)) %>%
    dplyr::select(perturb_name, perturb_summary_tbl) %>%
    tidyr::unnest(cols = c(perturb_summary_tbl))

  if (nrow(perturb_summaries) == 0) {
    return(tibble::tibble(from = character(), to = character(), pair_weight = numeric()))
  }

  power_cols <- colnames(perturb_summaries)[stringr::str_detect(colnames(perturb_summaries), "power")]
  if (length(power_cols) > 0) {
    perturb_cell_power <- perturb_summaries %>%
      dplyr::select(perturb_name, cell_group, dplyr::all_of(power_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(power_cols),
        names_to = "power_metric",
        values_to = "power_value"
      ) %>%
      dplyr::group_by(perturb_name, cell_group) %>%
      dplyr::summarize(
        max_power = ifelse(all(is.na(power_value)), NA_real_, max(power_value, na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    perturb_cell_power <- perturb_summaries %>%
      dplyr::select(perturb_name, cell_group) %>%
      dplyr::distinct() %>%
      dplyr::mutate(max_power = 1)
  }

  lost_states <- perturb_summaries %>%
    dplyr::filter(is_lost_when_present) %>%
    dplyr::select(perturb_name, from = cell_group) %>%
    dplyr::distinct()

  unaffected_states <- perturb_summaries %>%
    dplyr::filter(!is_lost_when_present) %>%
    dplyr::select(perturb_name, to = cell_group) %>%
    dplyr::distinct() %>%
    dplyr::left_join(perturb_cell_power, by = c("perturb_name", "to" = "cell_group")) %>%
    dplyr::filter(!is.na(max_power), max_power >= power_threshold)

  discordant_pairs <- lost_states %>%
    dplyr::inner_join(unaffected_states, by = "perturb_name", relationship = "many-to-many") %>%
    dplyr::filter(from != to) %>%
    dplyr::transmute(from, to, pair_weight = max_power) %>%
    dplyr::distinct()

  return(discordant_pairs)
}

#' Break forbidden source-target paths by greedy edge removals
#'
#' Greedy objective: pick the edge with highest
#' (total broken forbidden-path weight) / (edge deletion cost),
#' where forbidden paths are up to `k_paths` shortest paths for each discordant pair.
#'
#' @param state_graph Directed igraph to prune.
#' @param discordant_pairs Tibble with at least `from` and `to`; optional `pair_weight`.
#' @param k_paths Integer. Number of candidate shortest paths to consider per pair.
#' @param deletion_cost_attr Edge attribute used as deletion cost.
#' @param traversal_weight_attr Edge attribute used to rank shortest paths.
#'
#' @return A list with `graph` (pruned igraph) and `removed_edges` tibble.
#' @noRd
prune_discordant_paths_greedily <- function(state_graph,
                                            discordant_pairs,
                                            k_paths = 1,
                                            deletion_cost_attr = "total_path_score_supporting",
                                            traversal_weight_attr = "weight",
                                            lambda_edge = 0) {
  if (is.null(state_graph) || igraph::gsize(state_graph) == 0 || nrow(discordant_pairs) == 0) {
    return(list(
      graph = state_graph,
      removed_edges = tibble::tibble(iteration = integer(), from = character(), to = character(), score_ratio = numeric())
    ))
  }

  k_paths <- max(1, as.integer(k_paths))

  if (!"pair_weight" %in% colnames(discordant_pairs)) {
    discordant_pairs <- discordant_pairs %>% dplyr::mutate(pair_weight = 1)
  }

  graph <- state_graph
  removed_edges <- tibble::tibble(iteration = integer(), from = character(), to = character(), score_ratio = numeric())
  max_iters <- igraph::gsize(graph)

  collect_forbidden_paths <- function(g, pairs_tbl, k, weight_attr) {
    edge_df <- igraph::as_data_frame(g, what = "edges") %>% dplyr::select(from, to)
    edge_df$edge_key <- paste(edge_df$from, edge_df$to, sep = "||")
    edge_weights <- igraph::edge_attr(g, weight_attr)
    if (is.null(edge_weights)) {
      edge_df$edge_weight <- 1
    } else {
      edge_df$edge_weight <- edge_weights
      edge_df$edge_weight[is.na(edge_df$edge_weight)] <- 1
    }

    pair_path_rows <- lapply(seq_len(nrow(pairs_tbl)), function(i) {
      p <- pairs_tbl[i, ]
      from_node <- as.character(p$from[[1]])
      to_node <- as.character(p$to[[1]])
      pair_weight <- as.numeric(p$pair_weight[[1]])
      if (!(from_node %in% igraph::V(g)$name && to_node %in% igraph::V(g)$name)) {
        return(NULL)
      }
      all_paths <- igraph::all_simple_paths(g, from = from_node, to = to_node, mode = "out")
      if (length(all_paths) == 0) {
        return(NULL)
      }
      path_tbl <- tibble::tibble(path_vertices = all_paths) %>%
        dplyr::mutate(
          path_vertices = lapply(path_vertices, names),
          path_edges = lapply(path_vertices, function(vs) {
            if (length(vs) < 2) {
              return(tibble::tibble(from = character(), to = character()))
            }
            tibble::tibble(from = head(vs, -1), to = tail(vs, -1))
          }),
          path_hops = vapply(path_edges, nrow, integer(1)),
          path_weight_total = vapply(path_edges, function(pe) {
            if (nrow(pe) == 0) {
              return(0)
            }
            edge_df %>%
              dplyr::inner_join(pe, by = c("from", "to")) %>%
              dplyr::summarize(w = sum(edge_weight, na.rm = TRUE)) %>%
              dplyr::pull(w) %>%
              as.numeric()
          }, numeric(1))
        ) %>%
        dplyr::arrange(path_weight_total, path_hops) %>%
        dplyr::slice_head(n = k) %>%
        dplyr::mutate(
          pair_from = from_node,
          pair_to = to_node,
          pair_weight = pair_weight,
          path_score = pair_weight,
          edge_keys = lapply(path_edges, function(pe) paste(pe$from, pe$to, sep = "||"))
        ) %>%
        dplyr::select(pair_from, pair_to, path_score, edge_keys)
      return(path_tbl)
    })

    path_tbl <- dplyr::bind_rows(pair_path_rows)
    if (nrow(path_tbl) == 0) {
      return(path_tbl)
    }
    return(path_tbl %>% dplyr::mutate(path_id = dplyr::row_number()))
  }

  for (iter in seq_len(max_iters)) {
    forbidden_paths <- collect_forbidden_paths(graph, discordant_pairs, k_paths, traversal_weight_attr)
    if (nrow(forbidden_paths) == 0) {
      break
    }

    edge_scores <- forbidden_paths %>%
      dplyr::select(path_id, path_score, edge_keys) %>%
      tidyr::unnest(edge_keys) %>%
      dplyr::group_by(edge_keys) %>%
      dplyr::summarize(broken_path_weight = sum(path_score, na.rm = TRUE), .groups = "drop")

    curr_edges <- igraph::as_data_frame(graph, what = "edges") %>% dplyr::select(from, to)
    curr_edges$edge_key <- paste(curr_edges$from, curr_edges$to, sep = "||")
    deletion_cost <- igraph::edge_attr(graph, deletion_cost_attr)
    if (is.null(deletion_cost)) {
      curr_edges$deletion_cost <- 1
    } else {
      curr_edges$deletion_cost <- deletion_cost
      curr_edges$deletion_cost[is.na(curr_edges$deletion_cost)] <- 1
      curr_edges$deletion_cost <- pmax(abs(curr_edges$deletion_cost), 1e-8)
    }

    candidate_edges <- curr_edges %>%
      dplyr::inner_join(edge_scores, by = c("edge_key" = "edge_keys")) %>%
      dplyr::mutate(score_ratio = broken_path_weight / (deletion_cost + max(0, lambda_edge))) %>%
      dplyr::arrange(dplyr::desc(score_ratio), deletion_cost)

    if (nrow(candidate_edges) == 0) {
      break
    }

    edge_to_remove <- candidate_edges[1, ]
    graph <- igraph::delete_edges(
      graph,
      igraph::E(graph)[.from(edge_to_remove$from[[1]]) & .to(edge_to_remove$to[[1]])]
    )

    removed_edges <- dplyr::bind_rows(
      removed_edges,
      tibble::tibble(
        iteration = iter,
        from = edge_to_remove$from[[1]],
        to = edge_to_remove$to[[1]],
        score_ratio = edge_to_remove$score_ratio[[1]]
      )
    )
  }

  return(list(graph = graph, removed_edges = removed_edges))
}

#' Normalize discordant pruning configuration
#'
#' Supported config keys:
#' - `power_threshold`
#' - `prune_mode` ("existing" or "greedy")
#' - `K_paths`
#' - `lambda_edge`
#'
#' @noRd
normalize_discordant_pruning_config <- function(discordant_config = NULL,
                                                discordant_pruning_mode = "none",
                                                discordant_pruning_k = 1,
                                                discordant_pruning_power_threshold = 0,
                                                discordant_pruning_lambda_edge = 0) {
  cfg <- list(
    power_threshold = discordant_pruning_power_threshold,
    prune_mode = ifelse(identical(discordant_pruning_mode, "greedy"), "greedy", "existing"),
    K_paths = discordant_pruning_k,
    lambda_edge = discordant_pruning_lambda_edge
  )

  if (!is.null(discordant_config)) {
    if (!is.list(discordant_config)) {
      stop("discordant_config must be a list when provided.")
    }
    cfg[names(discordant_config)] <- discordant_config
  }

  cfg$prune_mode <- match.arg(as.character(cfg$prune_mode), c("existing", "greedy"))
  cfg$K_paths <- max(1, as.integer(cfg$K_paths))
  cfg$power_threshold <- as.numeric(cfg$power_threshold)
  cfg$lambda_edge <- as.numeric(cfg$lambda_edge)

  return(cfg)
}

#' Count remaining discordant reachability violations in a graph
#' @noRd
count_discordant_reachability_violations <- function(state_graph, discordant_pairs) {
  if (is.null(state_graph) || igraph::gsize(state_graph) == 0 || nrow(discordant_pairs) == 0) {
    return(0L)
  }

  sum(vapply(seq_len(nrow(discordant_pairs)), function(i) {
    from_node <- as.character(discordant_pairs$from[[i]])
    to_node <- as.character(discordant_pairs$to[[i]])
    if (!(from_node %in% igraph::V(state_graph)$name && to_node %in% igraph::V(state_graph)$name)) {
      return(FALSE)
    }
    d <- igraph::distances(state_graph, v = from_node, to = to_node, mode = "out")
    is.finite(d[1, 1])
  }, logical(1)))
}

#' Summarize edge support-related distributions for a graph
#' @noRd
summarize_edge_support_distributions <- function(state_graph, mode_label = "graph") {
  if (is.null(state_graph) || igraph::gsize(state_graph) == 0) {
    return(tibble::tibble(
      mode = character(),
      metric = character(),
      n_edges = integer(),
      mean = numeric(),
      median = numeric(),
      p90 = numeric(),
      max = numeric()
    ))
  }

  edge_tbl <- igraph::as_data_frame(state_graph, what = "edges")
  numeric_support_cols <- names(edge_tbl)[vapply(edge_tbl, is.numeric, logical(1))]
  numeric_support_cols <- numeric_support_cols[grepl("support|penalty", numeric_support_cols)]

  if (length(numeric_support_cols) == 0) {
    return(tibble::tibble(
      mode = mode_label,
      metric = "none",
      n_edges = nrow(edge_tbl),
      mean = NA_real_,
      median = NA_real_,
      p90 = NA_real_,
      max = NA_real_
    ))
  }

  purrr::map_dfr(numeric_support_cols, function(col_nm) {
    vals <- edge_tbl[[col_nm]]
    tibble::tibble(
      mode = mode_label,
      metric = col_nm,
      n_edges = length(vals),
      mean = mean(vals, na.rm = TRUE),
      median = stats::median(vals, na.rm = TRUE),
      p90 = as.numeric(stats::quantile(vals, probs = 0.9, na.rm = TRUE, names = FALSE)),
      max = max(vals, na.rm = TRUE)
    )
  })
}

#' Compare existing vs greedy discordant pruning modes
#'
#' Runs perturbation graph assembly twice with the same inputs:
#' - existing mode (`prune_mode = "existing"`)
#' - greedy mode (`prune_mode = "greedy"`)
#'
#' Returns edge-change counts/tables, remaining discordant reachability
#' violations, and support distribution summaries.
#'
#' @noRd
compare_discordant_pruning_modes <- function(ref_ccs,
                                             timeseries_graph,
                                             perturbation_ccm_tbl,
                                             discordant_config = list(),
                                             q_val = 0.01,
                                             start_time = NULL,
                                             stop_time = NULL,
                                             interval_col = "timepoint",
                                             interval_step = 2,
                                             min_interval = 4,
                                             max_interval = 24,
                                             log_abund_detection_thresh = -5,
                                             min_pathfinding_lfc = 0,
                                             links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                                             components = "partition",
                                             verbose = FALSE,
                                             edge_allowlist = NULL,
                                             edge_denylist = NULL,
                                             discordant_pruning_cost_attr = "total_path_score_supporting",
                                             newdata = tibble()) {
  base_cfg <- normalize_discordant_pruning_config(discordant_config = discordant_config)
  existing_cfg <- base_cfg
  existing_cfg$prune_mode <- "existing"
  greedy_cfg <- base_cfg
  greedy_cfg$prune_mode <- "greedy"

  g_existing <- assemble_transition_graph_from_perturbations(
    ref_ccs = ref_ccs,
    timeseries_graph = timeseries_graph,
    perturbation_ccm_tbl = perturbation_ccm_tbl,
    q_val = q_val,
    start_time = start_time,
    stop_time = stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    min_interval = min_interval,
    max_interval = max_interval,
    log_abund_detection_thresh = log_abund_detection_thresh,
    min_pathfinding_lfc = min_pathfinding_lfc,
    links_between_components = links_between_components,
    components = components,
    verbose = verbose,
    edge_allowlist = edge_allowlist,
    edge_denylist = edge_denylist,
    discordant_config = existing_cfg,
    discordant_pruning_cost_attr = discordant_pruning_cost_attr,
    newdata = newdata
  )

  g_greedy <- assemble_transition_graph_from_perturbations(
    ref_ccs = ref_ccs,
    timeseries_graph = timeseries_graph,
    perturbation_ccm_tbl = perturbation_ccm_tbl,
    q_val = q_val,
    start_time = start_time,
    stop_time = stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    min_interval = min_interval,
    max_interval = max_interval,
    log_abund_detection_thresh = log_abund_detection_thresh,
    min_pathfinding_lfc = min_pathfinding_lfc,
    links_between_components = links_between_components,
    components = components,
    verbose = verbose,
    edge_allowlist = edge_allowlist,
    edge_denylist = edge_denylist,
    discordant_config = greedy_cfg,
    discordant_pruning_cost_attr = discordant_pruning_cost_attr,
    newdata = newdata
  )

  edges_existing <- igraph::as_data_frame(g_existing, what = "edges") %>% dplyr::select(from, to) %>% dplyr::distinct()
  edges_greedy <- igraph::as_data_frame(g_greedy, what = "edges") %>% dplyr::select(from, to) %>% dplyr::distinct()

  edges_removed <- dplyr::anti_join(edges_existing, edges_greedy, by = c("from", "to"))
  edges_added <- dplyr::anti_join(edges_greedy, edges_existing, by = c("from", "to"))

  discordant_pairs <- collect_discordant_pairs_from_perturbation_summaries(
    perturbation_ccm_tbl = perturbation_ccm_tbl,
    power_threshold = base_cfg$power_threshold
  ) %>% dplyr::select(from, to, pair_weight) %>% dplyr::distinct()

  violations_existing <- count_discordant_reachability_violations(g_existing, discordant_pairs)
  violations_greedy <- count_discordant_reachability_violations(g_greedy, discordant_pairs)

  support_summary <- dplyr::bind_rows(
    summarize_edge_support_distributions(g_existing, "existing"),
    summarize_edge_support_distributions(g_greedy, "greedy")
  )

  list(
    config = base_cfg,
    edge_changes = list(
      n_removed = nrow(edges_removed),
      n_added = nrow(edges_added),
      removed = edges_removed,
      added = edges_added
    ),
    discordant_violations = tibble::tibble(
      mode = c("existing", "greedy"),
      n_remaining = c(violations_existing, violations_greedy)
    ),
    support_summary = support_summary,
    graphs = list(existing = g_existing, greedy = g_greedy)
  )
}

#' @export
get_perturbation_paths <- function(perturbation_ccm,
                                   perturb_summary_tbl,
                                   pathfinding_graph,
                                   delta_log_abund_loss_thresh = 0) {
  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  # old_omp_num_threads = single_thread_omp()
  # old_blas_num_threads = single_thread_blas()
  Sys.setenv("OMP_NUM_THREADS" = 1)
  Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

  tryCatch({
    lost_cell_groups <- perturb_summary_tbl %>%
      filter(loss_when_present < -delta_log_abund_loss_thresh) %>%
      filter(cell_group %in% igraph::V(pathfinding_graph)$name) %>%
      pull(cell_group) %>%
      unique()
    gained_cell_groups <- perturb_summary_tbl %>%
      filter(gain_when_present > delta_log_abund_loss_thresh) %>%
      filter(cell_group %in% igraph::V(pathfinding_graph)$name) %>%
      pull(cell_group) %>%
      unique()

    # lost_cell_groups = perturb_summary_tbl %>% filter(is_lost_when_present) %>% pull(cell_group) %>% unique
    # gained_cell_groups = perturb_summary_tbl %>% filter(is_gained_when_present) %>% pull(cell_group) %>% unique

    perturbed_cell_groups <- union(lost_cell_groups, gained_cell_groups)
    lost_subgraph <- igraph::subgraph(pathfinding_graph, lost_cell_groups)
    nodes_without_lost_parent <- igraph::V(lost_subgraph)[igraph::degree(lost_subgraph, mode = "in") == 0]$name

    # We will do pathfinding to nodes that are lost, starting from nodes that
    # are either lost and have no parents that are also lost, or nodes that are
    # gained.
    if (length(lost_cell_groups) == 0) {
      return(NA)
    }
    loss_neighborhood_graph <- igraph::subgraph(pathfinding_graph, lost_cell_groups)
    perturb_pathfinding_graph <- loss_neighborhood_graph

    if (length(gained_cell_groups) > 0) {
      gain_neighborhood_graph <- do.call(igraph::union, igraph::make_ego_graph(pathfinding_graph,
        order = 1,
        nodes = gained_cell_groups,
        mode = "out"
      ))
      perturb_pathfinding_graph <- igraph::union(perturb_pathfinding_graph, gain_neighborhood_graph)
    }

    perturb_pathfinding_graph <- igraph::subgraph(pathfinding_graph, igraph::V(perturb_pathfinding_graph)$name)

    perturb_pathfinding_graph <- igraph::as_data_frame(perturb_pathfinding_graph) %>%
      select(from, to) %>%
      distinct()

    perturb_pathfinding_graph <- left_join(perturb_pathfinding_graph,
      igraph::as_data_frame(pathfinding_graph),
      by = c("from", "to")
    )
    perturb_pathfinding_graph <- igraph::graph_from_data_frame(perturb_pathfinding_graph)

    perturb_pathfinding_graph_nodes <- igraph::V(perturb_pathfinding_graph)$name

    possible_origins <- union(nodes_without_lost_parent, gained_cell_groups)

    # origin_dest_pairs_to_test = tibble(cell_group=loss_neighborhood_graph_nodes) %>%
    #  tidyr::expand(cell_group, cell_group)
    origin_dest_pairs_to_test <- expand.grid(possible_origins, lost_cell_groups) %>% as_tibble()
    colnames(origin_dest_pairs_to_test) <- c("from", "to")
    origin_dest_pairs_to_test$from <- as.character(origin_dest_pairs_to_test$from)
    origin_dest_pairs_to_test$to <- as.character(origin_dest_pairs_to_test$to)
    origin_dest_pairs_to_test <- origin_dest_pairs_to_test %>% filter(from != to)

    origin_dest_pairs_to_test <- dplyr::left_join(origin_dest_pairs_to_test, perturb_summary_tbl %>% setNames(paste0("to_", names(.))), by = c("to" = "to_cell_group")) # %>%
    origin_dest_pairs_to_test <- dplyr::left_join(origin_dest_pairs_to_test, perturb_summary_tbl %>% setNames(paste0("from_", names(.))), by = c("from" = "from_cell_group")) # %>%

    concordant_fwd_loss_pairs <- origin_dest_pairs_to_test %>%
      # dplyr::filter(pcor < 0) %>% # do we just want negative again?
      dplyr::filter(is.na(to_is_lost_when_present) == FALSE & to_is_lost_when_present & from_peak_wt_time <= to_peak_wt_time)

    # Compute the shortest paths between each of those node pairs in the pathfinding graph
    concordant_fwd_loss_pairs <- concordant_fwd_loss_pairs %>%
      select(from, to) %>%
      distinct()
    message(paste("\tfinding shortest paths between", nrow(concordant_fwd_loss_pairs), "possible loss pairs"))

    paths_between_concordant_fwd_loss_pairs <- concordant_fwd_loss_pairs %>%
      mutate(path = furrr::future_map2(
        .f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
        .x = from, .y = to,
        perturb_pathfinding_graph,
        .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
        .progress = TRUE
      ))

    paths_between_concordant_fwd_loss_pairs <- paths_between_concordant_fwd_loss_pairs %>%
      filter(is.na(path) == FALSE)

    if (nrow(paths_between_concordant_fwd_loss_pairs) == 0) {
      return(NA)
    }
    message(paste("\tscoring", nrow(paths_between_concordant_fwd_loss_pairs), "paths between loss pairs"))
    paths_between_concordant_fwd_loss_pairs <- score_paths_for_perturbations(perturbation_ccm,
      paths_between_concordant_fwd_loss_pairs,
      loss_tbl = perturb_summary_tbl
    )
    return(paths_between_concordant_fwd_loss_pairs)
  }, error = function(e) {
    print(e)
    return(NA)
  }, finally = {
    # RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    # RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  })
}


# 2) Make sure this graph is acyclic by deleting problematic edges
# from https://github.com/sachsmc/causaloptim/blob/master/R/graph-utilities.R

#' Find cycles in a graph
#'
#' @param g an igraph object
#' @return A list of vectors of integers, indicating the vertex sequences for the cycles found in the graph
#' @export
find_cycles <- function(g) {
  Cycles <- NULL
  for (v1 in igraph::V(g)) {
    if (igraph::degree(g, v1, mode = "in") == 0) {
      next
    }
    GoodNeighbors <- igraph::neighbors(g, v1, mode = "out")
    GoodNeighbors <- GoodNeighbors[GoodNeighbors > v1]
    for (v2 in GoodNeighbors) {
      TempCyc <- lapply(igraph::all_simple_paths(g, v2, v1, mode = "out"), function(p) c(v1, p))
      TempCyc <- TempCyc[which(sapply(TempCyc, length) > 2)]
      TempCyc <- TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles <- c(Cycles, TempCyc)
    }
  }
  Cycles
}


#' score a path based on fitting a linear model of time ~ geodeisic distance
#' @param ccs cell count set
#' @param path_df path df
#' @noRd
measure_time_delta_along_path <- function(path_df, ccs, cells_along_path_df, interval_col = "timepoint") {
  # cds = ccs@cds
  vertices <- union(path_df$to, path_df$from) %>% unique()
  path_df <- path_df %>%
    arrange(distance_from_root) %>%
    mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df <- cells_along_path_df %>% filter(cell_group %in% vertices)

  # cells_along_path_df = normalized_counts(ccs, "size_only", pseudocount = 0)[vertices,] %>%
  #   as.matrix() %>%
  #   Matrix::t() %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column() %>%
  #   tidyr::pivot_longer(!matches("rowname")) %>%
  #   rename(sample=rowname, cell_group=name, num_cells=value)

  cells_along_path_df <- cells_along_path_df %>%
    left_join(colData(ccs) %>% as.data.frame() %>%
      select(sample, !!sym(interval_col)) %>%
      as_tibble(), by = c("sample" = "sample"))

  cells_along_path_df <- cells_along_path_df %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  cells_along_path_df <- cells_along_path_df %>% filter(num_cells > 0)
  cells_along_path_df$distance_from_root <- tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  cells_along_path_df$geodesic_dist <- tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$y <- cells_along_path_df[[interval_col]]

  path_model <- lm(y ~ geodesic_dist, data = cells_along_path_df, weights = cells_along_path_df$num_cells)

  path_model_tidied <- broom::tidy(path_model)
  path_model_glanced <- broom::glance(path_model)
  pm_stats <- tibble(
    time_dist_effect = unlist(path_model_tidied[2, "estimate"]),
    time_dist_effect_pval = unlist(path_model_tidied[2, "p.value"]),
    time_dist_model_adj_rsq = unlist(path_model_glanced[1, "adj.r.squared"]),
    time_dist_model_ncells = unlist(path_model_glanced[1, "nobs"]),
    path_length = max(path_df$geodesic_dist)
  )
  return(pm_stats)
  # return(coef(path_model)[["distance_from_root"]])
}
# debug(cells_along_path)

# # ' score a path based on fitting a linear model of perturbation ~ geodesic distance
# # ' @param ccs
# # ' @param path_df
# # ' @noRd
# measure_perturbation_freq_along_path <- function(path_df, ccs, cells_along_path_df, perturbation_col="knockout", interval_col="timepoint", batch_col=NULL) {
#
#
#   #print ("measuring perturbation score")
#   #cds = ccs@cds
#   vertices = union(path_df$to, path_df$from) %>% unique()
#   path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))
#
#   cells_along_path_df = cells_along_path_df %>% filter(cell_group %in% vertices)
#
#   num_batches = 1
#   if (is.null(batch_col)){
#     cells_along_path_df = cells_along_path_df %>%
#       left_join(colData(ccs) %>% as.data.frame %>%
#                   select(sample, !!sym(perturbation_col), !!(interval_col)) %>%
#                   as_tibble(), by = c("sample" = "sample"))
#   }else{
#     cells_along_path_df = cells_along_path_df %>%
#       left_join(colData(ccs) %>% as.data.frame %>%
#                   select(sample, !!sym(perturbation_col), !!(interval_col), !!(batch_col)) %>%
#                   as_tibble(), by = c("sample" = "sample"))
#     #print(head(cells_along_path_df))
#     num_batches = cells_along_path_df %>% pull(!!sym(batch_col)) %>% unique %>% length
#   }
#
#   cells_along_path_df = cells_along_path_df %>%
#     left_join(path_df %>% select(-from), by = c("cell_group" = "to"))
#
#   #print (cells_along_path_df)
#   cells_along_path_df = cells_along_path_df %>% mutate_if(is.numeric, tidyr::replace_na, replace = 0)
#   cells_along_path_df = cells_along_path_df %>% dplyr::filter (num_cells > 0)
#   #cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
#   #cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)
#
#   cells_along_path_df$perturb = cells_along_path_df[[perturbation_col]]
#   cells_along_path_df$timepoint = cells_along_path_df[[interval_col]]
#
#   if (is.null(batch_col) || num_batches == 1){
#     path_model = glm(#perturb ~ timepoint + geodesic_dist,
#       perturb ~ geodesic_dist,
#       data=cells_along_path_df,
#       weights=cells_along_path_df$num_cells,
#       family=binomial(),
#       singular.ok = FALSE)
#   }else{
#     cells_along_path_df$batch = cells_along_path_df[[batch_col]]
#     #print (head(cells_along_path_df))
#     path_model = glm(#perturb ~ timepoint + geodesic_dist,
#       perturb ~ batch + geodesic_dist,
#       #perturb ~ geodesic_dist,
#       data=cells_along_path_df,
#       weights=cells_along_path_df$num_cells,
#       family=binomial(),
#       singular.ok = FALSE)
#     #print (summary(path_model))
#   }
#
#   path_model_tidied = broom::tidy(path_model)
#   path_model_glanced = broom::glance(path_model)
#   path_model_glanced$pseudo.r.squared = 1 - (path_model_glanced$deviance / path_model_glanced$null.deviance)
#
#   dist_param_row = nrow(path_model_tidied)
#   pm_stats = tibble(perturb_dist_effect = unlist(path_model_tidied[dist_param_row, "estimate"]),
#                     perturb_dist_effect_pval = unlist(path_model_tidied[dist_param_row, "p.value"]),
#                     perturb_dist_model_adj_rsq = unlist(path_model_glanced[1, "pseudo.r.squared"]),
#                     perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"]))
#   # if (unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0) { # should not happen
#   #   print (cells_along_path_df %>% as.data.frame)
#   #   print (summary(path_model))
#   # }
#
#   if (path_model$converged == FALSE | unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0)
#   {
#     return(tibble(perturb_dist_effect = 0,
#                   perturb_dist_effect_pval = 1,
#                   perturb_dist_model_adj_rsq = 0,
#                   perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"])))
#   }
#   return(pm_stats)
#   #return(coef(path_model)[["distance_from_root"]])
# }
# #debug(cells_along_path)

#' @noRd
measure_perturbation_freq_along_path <- function(path_df, ccs, cells_along_path_df, perturbation_col = "knockout", interval_col = "timepoint", batch_col = NULL) {
  # print ("measuring perturbation score")
  # cds = ccs@cds
  vertices <- union(path_df$to, path_df$from) %>% unique()
  path_df <- path_df %>%
    arrange(distance_from_root) %>%
    mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df <- cells_along_path_df %>% filter(cell_group %in% vertices)

  num_batches <- 1
  if (is.null(batch_col)) {
    cells_along_path_df <- cells_along_path_df %>%
      left_join(colData(ccs) %>% as.data.frame() %>%
        select(sample, !!sym(perturbation_col), !!(interval_col)) %>%
        as_tibble(), by = c("sample" = "sample"))
  } else {
    cells_along_path_df <- cells_along_path_df %>%
      left_join(colData(ccs) %>% as.data.frame() %>%
        select(sample, !!sym(perturbation_col), !!(interval_col), !!(batch_col)) %>%
        as_tibble(), by = c("sample" = "sample"))
    # print(head(cells_along_path_df))
    num_batches <- cells_along_path_df %>%
      pull(!!sym(batch_col)) %>%
      unique() %>%
      length()
  }

  cells_along_path_df <- cells_along_path_df %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  # print (cells_along_path_df)
  cells_along_path_df <- cells_along_path_df %>% mutate_if(is.numeric, tidyr::replace_na, replace = 0)
  cells_along_path_df <- cells_along_path_df %>% dplyr::filter(num_cells > 0)
  # cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  # cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$perturb <- cells_along_path_df[[perturbation_col]]
  cells_along_path_df$timepoint <- cells_along_path_df[[interval_col]]

  if (is.null(batch_col) || num_batches == 1) {
    path_model <- glm( # perturb ~ timepoint + geodesic_dist,
      num_cells ~ cell_group + perturb,
      data = cells_along_path_df,
      # weights=cells_along_path_df$num_cells,
      family = quasipoisson(),
      singular.ok = FALSE
    )
  } else {
    cells_along_path_df$batch <- cells_along_path_df[[batch_col]]
    # print (head(cells_along_path_df))
    path_model <- glm( # perturb ~ timepoint + geodesic_dist,
      num_cells ~ batch + cell_group + perturb,
      # perturb ~ geodesic_dist,
      data = cells_along_path_df,
      # weights=cells_along_path_df$num_cells,
      family = quasipoisson(),
      singular.ok = FALSE
    )
    # print (summary(path_model))
  }

  path_model_tidied <- broom::tidy(path_model)
  path_model_glanced <- broom::glance(path_model)
  path_model_glanced$pseudo.r.squared <- 1 - (path_model_glanced$deviance / path_model_glanced$null.deviance)

  dist_param_row <- nrow(path_model_tidied)
  pm_stats <- tibble(
    perturb_dist_effect = unlist(path_model_tidied[dist_param_row, "estimate"]),
    perturb_dist_effect_pval = unlist(path_model_tidied[dist_param_row, "p.value"]),
    perturb_dist_model_adj_rsq = unlist(path_model_glanced[1, "pseudo.r.squared"]),
    perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"])
  )
  # if (unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0) { # should not happen
  #   print (cells_along_path_df %>% as.data.frame)
  #   print (summary(path_model))
  # }

  if (path_model$converged == FALSE | unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0) {
    return(tibble(
      perturb_dist_effect = 0,
      perturb_dist_effect_pval = 1,
      perturb_dist_model_adj_rsq = 0,
      perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"])
    ))
  }
  return(pm_stats)
  # return(coef(path_model)[["distance_from_root"]])
}
# debug(cells_along_path)


#' Identify the possible origins for each destination
#'
#' @noRd

build_timeseries_transition_graph <- function(ccm,
                                              extant_cell_type_df,
                                              pathfinding_graph,
                                              q_val = 0.01,
                                              start_time = NULL,
                                              stop_time = NULL,
                                              interval_col = "timepoint",
                                              interval_step = 2,
                                              min_interval = 4,
                                              max_interval = 24,
                                              min_pathfinding_lfc = 0,
                                              make_dag = FALSE,
                                              newdata = tibble(),
                                              log_abund_detection_thresh = -5,
                                              edge_allowlist = NULL,
                                              edge_denylist = NULL) {
  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  # old_omp_num_threads = single_thread_omp()
  # old_blas_num_threads = single_thread_blas()
  Sys.setenv("OMP_NUM_THREADS" = 1)
  Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

  # First, let's figure out when each cell type is present and
  # which ones emerge over the course of the caller's time interval
  if (is.null(start_time)) {
    start_time <- min(colData(ccm@ccs)[, interval_col])
  }
  if (is.null(stop_time)) {
    stop_time <- max(colData(ccm@ccs)[, interval_col])
  }

  timepoints <- seq(start_time, stop_time, interval_step)

  message("Estimating abundances over time interval")
  timepoint_pred_df <- estimate_abundances_over_interval(ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    newdata = newdata,
    min_log_abund = log_abund_detection_thresh
  )

  #' @noRd
  select_timepoints <- function(timepoint_pred_df, t1, t2, interval_col) {
    cond_x <- timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
    cond_y <- timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
    return(compare_abundances(ccm, cond_x, cond_y))
  }

  time_contrasts <- expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2 - t1) >= min_interval & (t2 - t1) <= max_interval) %>%
    tibble::as_tibble()

  message(paste("Comparing abundances over", nrow(time_contrasts), "timepoint contrasts"))
  relevant_comparisons <- time_contrasts %>%
    as_tibble() %>%
    mutate(comp_abund = furrr::future_map2(
      .f = select_timepoints,
      .x = t1,
      .y = t2,
      interval_col = interval_col,
      timepoint_pred_df = timepoint_pred_df,
      .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
      .progress = TRUE
    ))

  message(paste("Collecting relevant PLN network graph edges"))
  relevant_comparisons <- relevant_comparisons %>%
    mutate(rec_edges = purrr::map(
      .f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
      .x = comp_abund,
      ccm = ccm,
      log_abundance_thresh = log_abund_detection_thresh,
      edge_allowlist = edge_allowlist,
      edge_denylist = edge_denylist
    ))

  # mutate(rec_edges = furrr::future_map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
  #                                      .x = comp_abund,
  #                                      ccm = ccm,
  #                                      .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
  #                                      .progress=TRUE))

  message(paste("Finding reciprocal pairs"))
  relevant_comparisons <- relevant_comparisons %>%
    tidyr::unnest(rec_edges) %>%
    # dplyr::filter(pcor < 0) %>% # do we just want negative again?
    dplyr::filter((from_delta_log_abund > abs(min_pathfinding_lfc) & to_delta_log_abund < -abs(min_pathfinding_lfc)) |
      (to_delta_log_abund > abs(min_pathfinding_lfc) & from_delta_log_abund < -abs(min_pathfinding_lfc))) %>%
    dplyr::filter(from_delta_q_value < q_val & to_delta_q_value < q_val)

  if (nrow(relevant_comparisons) == 0) {
    stop("No reciprocal pairs of nodes found")
  }

  edge_union <- relevant_comparisons %>%
    select(from, to) %>%
    distinct()

  if (!is.null(edge_allowlist)) {
    edge_union <- edge_union %>%
      rbind(edge_allowlist) %>%
      select(from, to) %>%
      distinct()
  }
  if (!is.null(edge_denylist)) {
    edge_union <- edge_union %>%
      dplyr::anti_join(edge_denylist, by = c("from", "to"))
  }

  print(paste("finding shortest paths between ", nrow(edge_union), "pairs of nodes"))
  # print(head(relevant_comparisons))
  paths_for_relevant_edges <- edge_union %>%
    mutate(path = furrr::future_map2(
      .f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
      .x = from, .y = to,
      traversal_graph = pathfinding_graph,
      .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
      .progress = TRUE
    ))

  cells_along_path_df <- normalized_counts(ccm@ccs, "size_only", pseudocount = 0) %>%
    as.matrix() %>%
    Matrix::t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(!matches("rowname")) %>%
    dplyr::rename(sample = rowname, cell_group = name, num_cells = value)

  paths_for_relevant_edges <- paths_for_relevant_edges %>%
    filter(!is.na(path))

  if (nrow(paths_for_relevant_edges) == 0) {
    stop("No time-forward paths found")
  }

  print(paste("scoring", nrow(paths_for_relevant_edges), "paths for time flow..."))
  # print(head(paths_for_relevant_edges))

  paths_for_relevant_edges <- paths_for_relevant_edges %>%
    mutate(time_vs_distance_model_stats = furrr::future_map(
      .f = purrr::possibly(measure_time_delta_along_path, NA_character_),
      .x = path,
      ccs = ccm@ccs,
      cells_along_path_df = cells_along_path_df,
      interval_col = interval_col,
      .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
      .progress = TRUE
    )) %>%
    tidyr::unnest(time_vs_distance_model_stats)
  # cells_along_path(ccs, paths_for_relevant_edges$path[[1]] %>% dplyr::select(from, to))


  print("tabulating model stats...")
  # print(head(paths_for_relevant_edges))
  paths_for_relevant_edges <- paths_for_relevant_edges %>% mutate( # time_dist_model_score = time_dist_effect * time_dist_model_adj_rsq,
    time_dist_model_score = time_dist_model_ncells * time_dist_model_adj_rsq,
    path_score = time_dist_model_score
  )
  paths_for_relevant_edges <- paths_for_relevant_edges %>% mutate(time_dist_effect_qval = p.adjust(time_dist_effect_pval, method = "BH"))

  # print (paths_for_relevant_edges)
  selected_paths <- paths_for_relevant_edges %>%
    filter(time_dist_effect > 0 &
      time_dist_effect_qval < q_val &
      time_dist_model_adj_rsq > 0) %>%
    group_by(to) %>%
    # slice_max(dist_model_score, n=3) %>%
    ungroup() %>%
    arrange(desc(time_dist_model_score)) %>%
    dplyr::mutate(path_contrast = "time")

  # selected_paths %>% select(origin=from, destination=to, path, path_score) %>% tidyr::unnest(path) %>% arrange(origin, destination) %>% filter(from %in% c("20", "38") & to %in% c("20", "38"))  %>% print(n=1000)

  # G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = FALSE)

  # Add allowlisted edges if they aren't already present in selected_paths
  if (!is.null(edge_allowlist)) {
    # Get all edges from selected_paths
    selected_edges <- selected_paths %>%
      dplyr::select(from, to) %>%
      distinct()
    # Find allowlist edges not in selected_paths
    missing_allowlist_edges <- edge_allowlist %>%
      anti_join(selected_edges, by = c("from", "to"))
    if (nrow(missing_allowlist_edges) > 0) {
      # Add missing allowlist edges as new paths (single-edge paths)
      allowlist_paths <- missing_allowlist_edges %>%
        mutate(
          path = purrr::map2(from, to, ~ tibble(from = .x, to = .y, weight = NA, distance_from_root = 0)),
          path_contrast = "allowlist",
          path_score = 0,
          time_dist_model_score = 0,
          time_dist_effect = NA,
          time_dist_effect_pval = NA,
          time_dist_model_adj_rsq = NA,
          time_dist_model_ncells = NA,
          time_dist_effect_qval = NA
        )
      selected_paths <- bind_rows(selected_paths, allowlist_paths)
    }
  }

  if (!is.null(edge_denylist)) {
    # Remove any paths that contain an edge in the denylist
    deny_edges <- edge_denylist %>%
      select(from, to) %>%
      distinct()
    selected_paths <- selected_paths %>%
      mutate(
        has_deny_edge = purrr::map_lgl(
          path,
          function(p) {
            any(
              dplyr::semi_join(
                p %>% select(from, to),
                deny_edges,
                by = c("from", "to")
              ) %>% nrow() > 0
            )
          }
        )
      ) %>%
      filter(!has_deny_edge) %>%
      select(-has_deny_edge)
  }

  print("combining paths...")
  G <- select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = TRUE)

  # if (!is.null(G) && !is.na(G) && make_dag){
  if (!is.null(G) && make_dag) {
    # igraph::edge_attr(G, "support") = igraph::edge_attr(G, "total_perturb_path_score_supporting")
    cycle_breaking_scores <- igraph::edge_attr(G, "total_path_score_supporting")
    cycle_breaking_scores[is.na(cycle_breaking_scores)] <- 0
    igraph::edge_attr(G, "support") <- cycle_breaking_scores

    print(igraph::edge_attr(G, "support"))

    print("breaking cycles...")
    G <- platt:::break_cycles_in_state_transition_graph(G)

    G <- compute_min_path_cover(ccm, G, weight_attribute = "support")

    edge_support <- igraph::as_data_frame(G) %>% select(from, to)

    edge_support <- left_join(
      edge_support,
      selected_paths %>% select(-from, -to) %>% tidyr::unnest(path)
    )

    # print (edge_support)
    edge_support <- edge_support %>%
      group_by(from, to) %>%
      slice_max(time_dist_model_adj_rsq, n = 1)

    # print (edge_support)
    G <- igraph::graph_from_data_frame(edge_support, directed = TRUE, vertices = data.frame(id = igraph::V(G)$name))
  }

  print("Finished building timeseries graph")
  return(G)
}


#' @noRd
get_paths_between_recip_time_nodes <- function(ccm,
                                               pathfinding_graph,
                                               timepoint_pred_df,
                                               timepoints,
                                               min_interval,
                                               max_interval,
                                               q_val,
                                               min_pathfinding_lfc,
                                               interval_col,
                                               verbose = FALSE) {
  # timepoints = sort(unique(timepoint_pred_df[[interval_col]]))

  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  # old_omp_num_threads = single_thread_omp()
  # old_blas_num_threads = single_thread_blas()

  tryCatch({
    compare_timepoints <- function(timepoint_pred_df, t1, t2, interval_col) {
      cond_x <- timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
      cond_y <- timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
      return(compare_abundances(ccm, cond_x, cond_y))
    }

    time_contrasts <- expand.grid("t1" = timepoints, "t2" = timepoints) %>%
      filter(t1 < t2 & (t2 - t1) >= min_interval & (t2 - t1) <= max_interval)

    if (verbose) {
      message(paste("\tLooking at", nrow(time_contrasts), "time contrasts"))
    }

    # Find the nodes that undergo reciprocal fold changes between successive time points
    recip_time_node_pairs <- time_contrasts %>%
      mutate(comp_abund = furrr::future_map2(
        .f = compare_timepoints,
        .x = t1,
        .y = t2,
        interval_col = interval_col,
        timepoint_pred_df = timepoint_pred_df,
        .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
        .progress = TRUE
      )) %>%
      mutate(rec_edges = furrr::future_map(
        .f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
        .x = comp_abund,
        ccm = ccm,
        .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
        .progress = TRUE
      ))
    recip_time_node_pairs <- recip_time_node_pairs %>%
      tidyr::unnest(rec_edges)

    # Let's adjust all the p values to account for the fact that we've done many contrasts
    recip_time_node_pairs <- recip_time_node_pairs %>%
      dplyr::mutate(
        from_delta_q_value = p.adjust(from_delta_p_value),
        to_delta_p_value = p.adjust(to_delta_p_value)
      )

    if (verbose) {
      num_distinct_node_pairs <- recip_time_node_pairs %>%
        select(from, to) %>%
        distinct()
      message(paste("\t Found", nrow(num_distinct_node_pairs), "reciprocal node pairs"))
    }

    recip_time_node_pairs <- recip_time_node_pairs %>%
      # dplyr::filter(pcor < 0) %>% # do we just want negative again?
      dplyr::filter((from_delta_log_abund > abs(min_pathfinding_lfc) & to_delta_log_abund < -abs(min_pathfinding_lfc)) |
        (to_delta_log_abund > abs(min_pathfinding_lfc) & from_delta_log_abund < -abs(min_pathfinding_lfc))) %>%
      dplyr::filter(from_delta_q_value < q_val & to_delta_q_value < q_val)

    if (verbose) {
      num_distinct_node_pairs <- recip_time_node_pairs %>%
        select(from, to) %>%
        distinct()
      message(paste("\tOf these", nrow(num_distinct_node_pairs), "survived thresholding"))
    }


    if (verbose) {
      message("Computing shortest paths between reciprocal time nodes")
    }

    # Compute the shortest paths between each of those node pairs in the pathfinding graph
    recip_time_node_pairs <- recip_time_node_pairs %>%
      select(from, to) %>%
      distinct()
    paths_between_recip_time_nodes <- recip_time_node_pairs %>%
      mutate(path = furrr::future_map2(
        .f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
        .x = from, .y = to,
        pathfinding_graph,
        .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
        .progress = TRUE
      ))

    # paths_between_recip_time_nodes %>% filter (from == "15") %>% print

    cells_along_path_df <- normalized_counts(ccm@ccs, "size_only", pseudocount = 0) %>%
      as.matrix() %>%
      Matrix::t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(!matches("rowname")) %>%
      rename(sample = rowname, cell_group = name, num_cells = value)

    # Score each shortest path for correlation between time and distance along it
    paths_between_recip_time_nodes <- paths_between_recip_time_nodes %>%
      filter(!is.na(path))

    if (verbose) {
      message(paste("Found", nrow(paths_between_recip_time_nodes), " paths through pathfinding graph. Scoring for time-flow..."))
    }

    paths_between_recip_time_nodes <- paths_between_recip_time_nodes %>%
      mutate(time_vs_distance_model_stats = furrr::future_map(
        .f = purrr::possibly(measure_time_delta_along_path, NA_character_),
        .x = path,
        ccs = ccm@ccs,
        cells_along_path_df = cells_along_path_df,
        interval_col = interval_col,
        .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
        .progress = TRUE
      )) %>%
      tidyr::unnest(time_vs_distance_model_stats)

    paths_between_recip_time_nodes <- paths_between_recip_time_nodes %>% mutate(time_dist_model_score = time_dist_effect * time_dist_model_adj_rsq)
    paths_between_recip_time_nodes <- paths_between_recip_time_nodes %>% mutate(time_dist_effect_qval = p.adjust(time_dist_effect_pval, method = "BH"))
    # paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% filter (dist_effect > 0 & dist_effect_q_val < q_val & dist_model_adj_rsq > min_dist_vs_time_r_sq) %>%
    #  group_by(to) %>%
    #  #slice_max(dist_model_score, n=3) %>%
    #  ungroup() %>% arrange(desc(dist_model_score))
    return(paths_between_recip_time_nodes)
  }, error = function(e) {
    print(e)
    return(NA)
  }, finally = {
    # RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    # RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  })
  return(NA)
}



#' @noRd
get_timeseries_paths <- function(ccm,
                                 extant_cell_type_df,
                                 pathfinding_graph,
                                 start_time,
                                 stop_time,
                                 # perturbation_col,
                                 interval_col,
                                 interval_step,
                                 min_interval,
                                 max_interval,
                                 q_val = 0.01,
                                 min_pathfinding_lfc = 0,
                                 verbose = FALSE,
                                 newdata = tibble()) {
  # First, let's figure out when each cell type is present and
  # which ones emerge over the course of the caller's time interval
  if (is.null(start_time)) {
    start_time <- min(colData(ccm@ccs)[, interval_col])
  }
  if (is.null(stop_time)) {
    stop_time <- max(colData(ccm@ccs)[, interval_col])
  }

  # FIXME: shouldn't use hardcoded knockout here
  wt_timepoint_pred_df <- estimate_abundances_over_interval(ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    newdata = newdata
  )

  timepoints <- seq(start_time, stop_time, interval_step)


  paths_between_recip_time_nodes <- get_paths_between_recip_time_nodes(ccm,
    pathfinding_graph,
    wt_timepoint_pred_df,
    timepoints,
    min_interval,
    max_interval,
    q_val,
    interval_col,
    min_pathfinding_lfc = min_pathfinding_lfc,
    verbose = verbose
  )
  return(paths_between_recip_time_nodes)
}

#' @noRd
select_paths_from_pathfinding_graph <- function(pathfinding_graph, selected_paths, allow_cycles = FALSE) {
  support_tbl <- selected_paths %>%
    dplyr::select(path_contrast, path_score, path) %>%
    tidyr::unnest(path)
  G <- igraph::graph_from_data_frame(support_tbl %>% select(from, to) %>% distinct(), vertices = data.frame(id = igraph::V(pathfinding_graph)$name))
  edge_support_summary <- support_tbl %>%
    group_by(from, to) %>%
    dplyr::select(from, to, path_score, path_contrast) %>%
    distinct() %>%
    summarize(support = sum(path_score))
  edge_support_summary <- edge_support_summary %>% mutate(support = ifelse(is.na(support), 0, support))
  annotated_G <- igraph::as_data_frame(G) %>%
    as_tibble() %>%
    dplyr::select(from, to) %>%
    distinct() %>%
    left_join(edge_support_summary)
  annotated_G <- igraph::graph_from_data_frame(annotated_G, directed = TRUE, vertices = data.frame(id = igraph::V(G)$name))
  return(annotated_G)

  # G = pathfinding_graph
  # G = igraph::delete_edges(G, igraph::E(G))
  #
  # for (i in (1:nrow(selected_paths))){
  #   next_path = selected_paths$path[[i]] %>% select(from, to)
  #   next_path_graph = next_path %>% igraph::graph_from_data_frame(directed=TRUE)
  #   G_prime = igraph::union(G, next_path_graph)
  #   if (allow_cycles | length(find_cycles(G_prime)) == 0){
  #     G = G_prime
  #   }else{
  #     # Debug:
  #     #print ("skipping due to cycle:")
  #     #print (selected_paths[i,] %>% select(-path))
  #   }
  # }
  # igraph::E(G)$transcriptome_dist = igraph::E(pathfinding_graph)[igraph::E(G)]$weight
  #
  # #print (selected_paths)
  # support_tbl= selected_paths %>% dplyr::select(path_contrast, path_score, path) %>% tidyr::unnest(path)
  # edge_support_summary = support_tbl %>% group_by(from, to) %>% dplyr::select(from, to, path_score, path_contrast)  %>%
  #   distinct() %>% summarize(support=sum(path_score))
  #
  # edge_support_summary = edge_support_summary %>% mutate(support = ifelse(is.na(support), 0, support))
  # #print (edge_support_summary %>% as.data.frame)
  #
  # #annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% left_join(edge_support)
  # annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% distinct() %>% left_join(edge_support_summary)
  # annotated_G = igraph::graph_from_data_frame(annotated_G, directed=TRUE, vertices=data.frame(id=igraph::V(G)$name))
  #
  # return (annotated_G)
}

#' Break cycles in a state graph
#'
#' @export
#'
break_cycles_in_state_transition_graph <- function(state_graph, support_attribute) {
  # print ("FAS weights")
  # print (igraph::E(state_graph)$support)

  cycle_breaking_scores <- igraph::edge_attr(state_graph, support_attribute)
  cycle_breaking_scores[is.na(cycle_breaking_scores)] <- 0

  # NOTE: squaring the scores should ensure that the cycle breaking procedure prefers to keep a few well supported edges even if that means removing lots of edges with lower support
  cycle_breaking_scores <- cycle_breaking_scores^2
  # igraph::edge_attr(state_graph, "support") = cycle_breaking_scores

  # print (igraph::edge_attr(state_graph, "support"))

  fas <- igraph::feedback_arc_set(state_graph, weights = cycle_breaking_scores)
  # print ("feedback arc set")
  # print (fas)
  state_graph <- igraph::delete_edges(state_graph, fas)
  # print ("checking for cycles")
  # curr_cycles = find_cycles(state_graph)
  # if (length(curr_cycles) > 0){
  # print ("warning: cycles remain")
  # print (curr_cycles)
  # }
  return(state_graph)
}

#' @noRd
get_pcor_between_pair <- function(perturb_model_pcor_matrix, from, to) {
  # edges_to_test = mapply(function(x, y) { return(c(x,y))}, from, to)
  pcor_vals <- mapply(function(x, y) {
    return(perturb_model_pcor_matrix[x, y])
  }, from, to)
  # pcor_val = igraph::E(pcor_graph)[igraph::get.edge.ids(pcor_graph, edges_to_test)]$weight
  # pcor_val = as.vector(perturb_model_pcor_matrix[edges_to_test])
  return(pcor_vals)
}


# score_path_for_loss = function(path, loss_tbl, effect_threshold=1){
#   cell_groups_on_path = tibble(cell_group=unique(c(path$from, path$to)))
#   cell_groups_on_path = inner_join(cell_groups_on_path, loss_tbl, by="cell_group")
#
#   # FIXME: we should find a way to incorporate effects from all the nodes along the path.
#   # cell_groups_on_path = cell_groups_on_path %>% filter(is_lost_when_present)
#
#   # except for the first node, find any cells that are above threshold of one
#   # even if they have a loss value, consider these positive
#   pos_cells = cell_groups_on_path[-1,] %>% filter(gain_when_present > effect_threshold)
#
#
#
#   cell_groups_on_path = cell_groups_on_path[-1,] #%>%
#   # group_by(cell_group) %>%
#   # mutate(when_present = sum(loss_when_present,  gain_when_present,na.rm=T),
#   #        when_present_tvalue = sum(loss_when_present_tvalue, gain_when_present_tvalue,na.rm=T))
#
#   # if there are lost cells, calculate the significance of the loss path
#   # should this be above 1 ?
#   if (nrow(cell_groups_on_path) > 1 & nrow(pos_cells) == 0){
#     cell_loss_stats_on_path = cell_groups_on_path %>%
#       summarize(perturb_dist_effect = sum(loss_when_present, na.rm=TRUE),
#                 perturb_dist_effect_pval = pnorm(sum(loss_when_present_tvalue, na.rm=TRUE), mean = 0,
#                                                  sd = nrow(cell_groups_on_path)),
#                 perturb_dist_model_adj_rsq = NA,
#                 perturb_dist_model_ncells = NA)
#   }else{
#     cell_loss_stats_on_path = tibble(perturb_dist_effect = NA,
#                                      perturb_dist_effect_pval = 1,
#                                      perturb_dist_model_adj_rsq = NA,
#                                      perturb_dist_model_ncells = NA)
#   }
#
#
#   return(cell_loss_stats_on_path)
# }


score_path_for_loss <- function(path, loss_tbl) {
  cell_groups_on_path <- tibble(cell_group = unique(c(path$from, path$to)))
  cell_groups_on_path <- inner_join(cell_groups_on_path, loss_tbl, by = "cell_group")

  # FIXME: we should find a way to incorporate effects from all the nodes along the path.
  # cell_groups_on_path = cell_groups_on_path %>% filter(is_lost_when_present)

  # except for the first node, make sure it isnt significantly positive
  # pos_cells = cell_groups_on_path[-1,] %>% filter(is_gain_when_present)
  pos_cells <- cell_groups_on_path[-1, ] %>% filter(!is.na(gain_when_present))
  # print('hi')
  cell_groups_on_path <- cell_groups_on_path[-1, ] %>%
    group_by(cell_group) %>%
    mutate(
      when_present = sum(loss_when_present, gain_when_present, na.rm = T),
      when_present_tvalue = sum(loss_when_present_tvalue, gain_when_present_tvalue, na.rm = T)
    )

  if (nrow(cell_groups_on_path) > 0 & nrow(pos_cells) == 0) {
    cell_loss_stats_on_path <- cell_groups_on_path %>%
      summarize(
        perturb_dist_effect = sum(when_present, na.rm = TRUE),
        perturb_dist_effect_pval = pnorm(sum(when_present_tvalue, na.rm = TRUE),
          mean = 0,
          sd = nrow(cell_groups_on_path)
        ),
        perturb_dist_model_adj_rsq = NA,
        perturb_dist_model_ncells = NA
      )
  } else {
    cell_loss_stats_on_path <- tibble(
      perturb_dist_effect = NA,
      perturb_dist_effect_pval = 1,
      perturb_dist_model_adj_rsq = NA,
      perturb_dist_model_ncells = NA
    )
  }


  return(cell_loss_stats_on_path)
}


#' @noRd
score_paths_for_perturbations <- function(perturbation_ccm,
                                          paths_between_recip_time_nodes,
                                          loss_tbl) {
  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  # old_omp_num_threads = single_thread_omp()
  # old_blas_num_threads = single_thread_blas()

  tryCatch({
    paths_between_perturb_vs_wt_node_pairs <- paths_between_recip_time_nodes %>%
      filter(!is.na(path))

    paths_between_perturb_vs_wt_node_pairs <- paths_between_recip_time_nodes %>%
      # mutate(pcor_between_node_pair =  get_pcor_between_pair(perturb_model_pcor_matrix, from, to)) %>%
      mutate(perturb_vs_distance_model_stats = purrr::map(
        .f = purrr::possibly(score_path_for_loss, NA_character_),
        .x = path,
        loss_tbl # ,
        # .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
        # .progress=TRUE
      )) %>%
      filter(!is.na(perturb_vs_distance_model_stats)) %>%
      tidyr::unnest(perturb_vs_distance_model_stats)
    paths_between_perturb_vs_wt_node_pairs <- paths_between_perturb_vs_wt_node_pairs %>% mutate(perturb_dist_effect_qval = p.adjust(perturb_dist_effect_pval, method = "BH"))
    return(paths_between_perturb_vs_wt_node_pairs)
  }, error = function(e) {
    print(e)
    return(NA)
  }, finally = {
    # RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    # RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  })
  return(NA)
}

#' @noRd
compare_ko_to_wt_at_timepoint <- function(tp, perturbation_ccm, wt_pred_df, ko_pred_df, interval_col) {
  cond_wt <- wt_pred_df %>% filter(!!sym(interval_col) == tp)
  cond_ko <- ko_pred_df %>% filter(!!sym(interval_col) == tp)
  return(compare_abundances(perturbation_ccm, cond_wt, cond_ko))
}

#'
#' @noRd
estimate_loss_timing <- function(perturbation_ccm,
                                 start_time,
                                 stop_time,
                                 interval_step,
                                 control_ccm = perturbation_ccm,
                                 control_start_time = start_time,
                                 control_stop_time = stop_time,
                                 log_abund_detection_thresh = -5,
                                 delta_log_abund_loss_thresh = 0,
                                 interval_col = "timepoint",
                                 q_val = 0.01,
                                 with_ties = FALSE,
                                 newdata = tibble()) {
  fraction_of_presence_window_lost_thresh <- 0.5

  if (nrow(newdata) > 0) {
    newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
    newdata_mt <- cross_join(tibble(knockout = TRUE), newdata)
  } else {
    newdata_wt <- tibble(knockout = FALSE)
    newdata_mt <- tibble(knockout = TRUE)
  }

  wt_timepoint_pred_df <- hooke:::estimate_abundances_over_interval(perturbation_ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    newdata = newdata_wt
  )
  ko_timepoint_pred_df <- hooke:::estimate_abundances_over_interval(perturbation_ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    newdata = newdata_mt
  )

  timepoints <- seq(start_time, stop_time, interval_step)

  perturb_vs_wt_nodes <- tibble(t1 = timepoints) %>%
    mutate(comp_abund = purrr::map(
      .f = compare_ko_to_wt_at_timepoint,
      .x = t1,
      perturbation_ccm = perturbation_ccm,
      interval_col = interval_col,
      wt_pred_df = wt_timepoint_pred_df,
      ko_pred_df = ko_timepoint_pred_df
    )) %>%
    tidyr::unnest(comp_abund)

  extant_wt_tbl <- get_extant_cell_types(control_ccm,
    control_start_time,
    control_stop_time,
    log_abund_detection_thresh = log_abund_detection_thresh,
    newdata = newdata_wt
  )

  changes_when_present_in_wt <- left_join(
    perturb_vs_wt_nodes %>% select(cell_group,
      wt_time_present = t1,
      delta_log_abund_when_present = delta_log_abund,
      delta_log_abund_when_present_se = delta_log_abund_se,
      power_when_present = power,
      log_abund_wt = log_abund_x,
      delta_q_value
    ),
    extant_wt_tbl,
    by = c("cell_group" = "cell_group", "wt_time_present" = "timepoint")
  ) %>%
    # mutate(peak_time_in_ctrl_within_perturb_time_range = tidyr::replace_na(peak_time_in_ctrl_within_perturb_time_range, FALSE)) %>%
    mutate(
      delta_q_value = ifelse(is.na(delta_q_value), 1, delta_q_value),
      delta_log_abund_when_present = ifelse(is.na(delta_log_abund_when_present), 0, delta_log_abund_when_present)
    ) %>%
    mutate(is_lost_when_present = present_above_thresh & delta_log_abund_when_present < -abs(delta_log_abund_loss_thresh)) %>%
    mutate(is_gained_when_present = present_above_thresh & delta_log_abund_when_present > -abs(delta_log_abund_loss_thresh))
  # loss_when_present_in_wt = loss_when_present_in_wt %>% group_by(cell_group) %>% slice_min(peak_wt_time, n=1, with_ties=with_ties)


  # FIXME: maybe should refactor the code below into another function that summarizes contrast over intervals:
  # num samples
  n <- nrow(model(perturbation_ccm)$fitted)
  # num parameters
  k <- length(rownames(coef(perturbation_ccm@best_full_model)))
  df.r <- n - k - 1
  # df_correction = sqrt(n / (n - k - 1))
  loss_summary_tbl <- changes_when_present_in_wt %>%
    filter(is_lost_when_present) %>%
    # filter(wt_time_present %in% timepoints) %>%
    group_by(cell_group) %>%
    summarize(
      loss_when_present = weighted.mean(delta_log_abund_when_present, percent_max_abund, na.rm = T),
      loss_when_present_se = weighted.mean(delta_log_abund_when_present_se, percent_max_abund, na.rm = T),
      loss_when_present_power = weighted.mean(power_when_present, percent_max_abund, na.rm = T),
      loss_when_present_tvalue = weighted.mean(delta_log_abund_when_present / delta_log_abund_when_present_se, percent_max_abund, na.rm = T),
      loss_when_present_tvalue_df = df.r,
      loss_when_present_p_value = 2 * pt(-abs(loss_when_present_tvalue), loss_when_present_tvalue_df)
    ) %>%
    mutate(loss_when_present_q_val = p.adjust(loss_when_present_p_value, method = "bonferroni"))

  gain_summary_tbl <- changes_when_present_in_wt %>%
    filter(is_gained_when_present) %>%
    group_by(cell_group) %>%
    summarize(
      gain_when_present = weighted.mean(delta_log_abund_when_present, percent_max_abund, na.rm = T),
      gain_when_present_se = weighted.mean(delta_log_abund_when_present_se, percent_max_abund, na.rm = T),
      gain_when_present_power = weighted.mean(power_when_present, percent_max_abund, na.rm = T),
      gain_when_present_tvalue = weighted.mean(delta_log_abund_when_present / delta_log_abund_when_present_se, percent_max_abund, na.rm = T),
      gain_when_present_tvalue_df = df.r,
      gain_when_present_p_value = 2 * pt(-abs(gain_when_present_tvalue), gain_when_present_tvalue_df)
    ) %>%
    mutate(gain_when_present_q_val = p.adjust(gain_when_present_p_value, method = "bonferroni"))


  change_summary_tbl <- changes_when_present_in_wt %>%
    select(cell_group) %>%
    distinct()
  change_summary_tbl <- left_join(change_summary_tbl, loss_summary_tbl, by = "cell_group")
  change_summary_tbl <- left_join(change_summary_tbl, gain_summary_tbl, by = "cell_group")
  change_summary_tbl <- change_summary_tbl %>%
    mutate(
      loss_when_present_q_val = ifelse(is.na(loss_when_present_q_val), 1, loss_when_present_q_val),
      loss_when_present = ifelse(is.na(loss_when_present), NA, loss_when_present),
      is_lost_when_present = loss_when_present_q_val < q_val,
      gain_when_present_q_val = ifelse(is.na(gain_when_present_q_val), 1, gain_when_present_q_val),
      gain_when_present = ifelse(is.na(gain_when_present), NA, gain_when_present),
      is_gained_when_present = gain_when_present_q_val < q_val
    )

  peak_wt_abundance <- estimate_abundances_over_interval(control_ccm,
    control_start_time,
    control_stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    newdata = newdata_wt
  ) %>%
    group_by(cell_group) %>%
    slice_max(log_abund, n = 1) %>%
    select(cell_group, peak_wt_time = !!sym(interval_col))

  change_summary_tbl <- left_join(change_summary_tbl,
    peak_wt_abundance,
    by = c("cell_group")
  )

  return(change_summary_tbl)
}


#' @noRd
compute_min_path_cover <- function(ccm, G, weight_attribute = "weight") {
  igraph::edge_attr(G, "weight") <- igraph::edge_attr(G, weight_attribute)
  transitive.closure <- function(g, mat = FALSE, loops = TRUE) {
    g <- igraph::as_adjacency_matrix(g, attr = "weight")

    n <- ncol(g)

    matExpIterativ <- function(x, pow, y = x, z = x, i = 1) {
      while (i < pow) {
        z <- z %*% x
        y <- y + z
        i <- i + 1
      }
      return(y)
    }

    h <- matExpIterativ(g, n)
    h <- (h > 0) * 1
    dimnames(h) <- dimnames(g)
    if (!loops) diag(h) <- rep(0, n) else diag(h) <- rep(1, n)
    if (!mat) h <- igraph::graph_from_adjacency_matrix(h, weighted = TRUE) # h <- as(h,"graphNEL")
    return(h)
  }
  G_tr <- transitive.closure(G, loops = F)

  G_split <- G_tr %>% igraph::as_data_frame(what = "edges")
  G_split$weight <- NULL

  cov_graph <- hooke:::return_igraph(model(ccm, "reduced"))
  pcor_mat <- cov_graph %>% igraph::as_adjacency_matrix(attr = "weight")
  pcor_mat_summ <- Matrix::summary(pcor_mat)
  pcor_mat <- data.frame(
    from = rownames(pcor_mat)[pcor_mat_summ$i],
    to = colnames(pcor_mat)[pcor_mat_summ$j],
    weight = pcor_mat_summ$x
  )
  G_split <- left_join(G_split, pcor_mat) %>% tidyr::replace_na(list(weight = 0))
  G_split$weight <- abs(G_split$weight)

  G_nodes <- union(G_split$from, G_split$to)
  split_node_metadata <- data.frame(id = c(
    stringr::str_c("left_", union(G_split$from, G_split$to)),
    stringr::str_c("right_", union(G_split$from, G_split$to))
  ))
  G_split$from <- stringr::str_c("left_", G_split$from)
  G_split$to <- stringr::str_c("right_", G_split$to)
  G_split <- igraph::graph_from_data_frame(G_split %>% dplyr::select(from, to, weight), directed = FALSE, vertices = split_node_metadata)

  igraph::V(G_split)$type <- grepl("left", igraph::V(G_split)$name)

  mbm <- maxmatching::maxmatching(G_split, weighted = TRUE)
  mbm$matching
  matching <- data.frame(dest_node = mbm$matching, orig_node = names(mbm$matching))
  matching <- subset(matching, grepl("left", orig_node) & grepl("right", dest_node))
  matching <- matching %>%
    dplyr::select(orig_node, dest_node) %>%
    mutate(
      orig_node = stringr::str_replace_all(orig_node, "left_", ""),
      dest_node = stringr::str_replace_all(dest_node, "right_", "")
    ) %>%
    dplyr::rename(from = orig_node, to = dest_node)

  # node_dag = G_tr
  # FIXME: use real weights here:
  # igraph::E(node_dag)$weight = 1

  node_dag <- G

  possible_origins <- names(which(igraph::degree(node_dag, mode = "in") == 0))
  possible_termini <- names(which(igraph::degree(node_dag, mode = "out") == 0))
  node_dag <- igraph::add_vertices(node_dag, 2, attr = list("name" = c("source", "sink")))
  source_edge_df <- data.frame(from = "source", to = possible_origins)
  node_dag <- igraph::union(node_dag, source_edge_df %>% igraph::graph_from_data_frame())
  sink_edge_df <- data.frame(from = possible_termini, to = "sink")
  node_dag <- igraph::union(node_dag, sink_edge_df %>% igraph::graph_from_data_frame())

  # replace NA values with weight 0
  igraph::edge_attr(node_dag, "weight") <- tidyr::replace_na(igraph::edge_attr(node_dag, "weight"), 0)
  # print(node_dag)
  paths_from_chains <- matching %>%
    as_tibble() %>%
    mutate(chain_leg = furrr::future_map2(
      .f = purrr::possibly(hooke:::get_shortest_path, NULL),
      .x = from, .y = to,
      node_dag,
      .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
      .progress = TRUE
    ))

  covered_graph <- paths_from_chains %>%
    select(chain_leg) %>%
    tidyr::unnest(c(chain_leg)) %>%
    igraph::graph_from_data_frame(vertices = data.frame(id = igraph::V(G_tr)$name))

  chain_heads <- names(which(igraph::degree(covered_graph, mode = "in") == 0))
  chain_tails <- names(which(igraph::degree(covered_graph, mode = "out") == 0))

  source_edge_df <- data.frame(from = "source", to = chain_heads)
  sink_edge_df <- data.frame(from = chain_tails, to = "sink")


  paths_to_chain_heads <- source_edge_df %>%
    as_tibble() %>%
    mutate(chain_leg = furrr::future_map2(
      .f = purrr::possibly(hooke:::get_shortest_path, NULL),
      .x = from, .y = to,
      node_dag,
      .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
      .progress = TRUE
    ))

  paths_to_chain_tails <- sink_edge_df %>%
    as_tibble() %>%
    mutate(chain_leg = furrr::future_map2(
      .f = purrr::possibly(hooke:::get_shortest_path, NULL),
      .x = from, .y = to,
      node_dag,
      .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
      .progress = TRUE
    ))

  covered_graph <- paths_to_chain_heads %>%
    select(chain_leg) %>%
    tidyr::unnest(c(chain_leg)) %>%
    igraph::graph_from_data_frame(vertices = data.frame(id = igraph::V(node_dag)$name)) %>%
    igraph::union(covered_graph)

  # covered_graph = paths_to_chain_tails %>% select(chain_leg) %>%
  #   tidyr::unnest() %>%
  #   igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(node_dag)$name)) %>%
  #   igraph::union(covered_graph)

  covered_graph <- igraph::simplify(covered_graph)

  covered_graph <- igraph::delete_vertices(covered_graph, c("source", "sink"))

  return(covered_graph)
}


#' assemble a state transition graph from a timeseries
#' @export
assemble_timeseries_transitions <- function(ccm,
                                            q_val = 0.01,
                                            start_time = NULL,
                                            stop_time = NULL,
                                            interval_col = "timepoint",
                                            interval_step = 2,
                                            min_interval = 4,
                                            max_interval = 24,
                                            log_abund_detection_thresh = -5,
                                            min_pathfinding_lfc = 0,
                                            make_dag = FALSE,
                                            links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                                            components = "partition",
                                            edge_allowlist = NULL,
                                            edge_denylist = NULL,
                                            force_allowlist = FALSE,
                                            newdata = tibble()) {
  message("Determining extant cell types")
  extant_cell_type_df <- get_extant_cell_types(ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    log_abund_detection_thresh = log_abund_detection_thresh,
    newdata = newdata
  )

  # Now let's set up a directed graph that links the states between which cells *could* directly
  # transition. If we don't know the direction of flow, add edges in both directions. The idea is
  # that we will find shortest paths over this graph between destination states and their plausible
  # origin states, then choose the best origins for each destination.

  message("Initializing pathfinding graph")
  pathfinding_graph <- init_pathfinding_graph(ccm,
    extant_cell_type_df,
    links_between_components = links_between_components,
    components = components,
    edge_allowlist = edge_allowlist,
    edge_denylist = edge_denylist,
    force_allowlist = force_allowlist
  )


  G <- build_timeseries_transition_graph(ccm,
    extant_cell_type_df,
    pathfinding_graph,
    q_val,
    start_time,
    stop_time,
    interval_col,
    interval_step,
    min_interval,
    max_interval,
    min_pathfinding_lfc = min_pathfinding_lfc,
    make_dag = make_dag,
    newdata = newdata,
    log_abund_detection_thresh = log_abund_detection_thresh,
    edge_allowlist = edge_allowlist,
    edge_denylist = edge_denylist
  )


  # FIXME: Consider moving this step into build_timeseries_transition_graph()?
  if (!is.null(G)) {
    igraph::V(G)$cell_group <- ccm@ccs@info$cell_group
  }

  return(G)
}

#' helper function for assessing the gains and losses of each cell type in a
#' perturbation
#' @noRd
collect_perturb_effects <- function(perturbation_ccm,
                                    time_window,
                                    interval_col,
                                    q_val = 0.01,
                                    interval_step = 2,
                                    min_interval = 4,
                                    max_interval = 24,
                                    log_abund_detection_thresh = -5,
                                    min_lfc = 0,
                                    newdata = tibble()) {
  start_time <- min(as.numeric(time_window$start_time))
  stop_time <- min(as.numeric(time_window$stop_time))

  perturb_effect_summary <- estimate_loss_timing(perturbation_ccm,
    start_time = start_time,
    stop_time = stop_time,
    interval_step = interval_step,
    interval_col = interval_col,
    log_abund_detection_thresh = log_abund_detection_thresh,
    q_val = q_val,
    delta_log_abund_loss_thresh = min_lfc,
    newdata = newdata
  )
  return(perturb_effect_summary)
}

#' Assess the cell type gains and losses following a perturbation.
#'
#' Returns a table of summary stats for gains and losses that's used to assemble
#' the cell types into a depencency graph
#'
#' @export
assess_perturbation_effects <- function(perturbation_ccm_tbl,
                                        q_val = 0.01,
                                        start_time = NULL,
                                        stop_time = NULL,
                                        interval_col = "timepoint",
                                        interval_step = 2,
                                        min_interval = 4,
                                        max_interval = 24,
                                        log_abund_detection_thresh = -5,
                                        min_lfc = 0,
                                        verbose = FALSE,
                                        newdata = tibble()) {
  perturbation_ccm_tbl <- perturbation_ccm_tbl %>%
    dplyr::mutate(perturb_summary_tbl = purrr::map2(
      .f = purrr::possibly(
        collect_perturb_effects, NA_real_
      ),
      .x = perturb_ccm,
      .y = perturb_time_window,
      interval_col = interval_col,
      interval_step = interval_step,
      q_val = q_val,
      min_lfc = min_lfc,
      log_abund_detection_thresh = log_abund_detection_thresh,
      newdata = newdata
    ))


  # Perform a global correction for multiple testing
  perturbs <- perturbation_ccm_tbl %>%
    filter(!is.na(perturb_summary_tbl)) %>%
    dplyr::select(perturb_name, perturb_summary_tbl)

  # Start from the q values, as these are already corrected for the number of cell types in the model
  perturbs <- perturbs %>%
    filter(!is.na(perturb_summary_tbl)) %>%
    tidyr::unnest(cols = c(perturb_summary_tbl)) %>%
    ungroup() %>%
    mutate(
      loss_when_present_q_val = p.adjust(loss_when_present_q_val, method = "bonferroni"),
      loss_when_present_q_val = ifelse(is.na(loss_when_present_q_val), 1, loss_when_present_q_val),
      is_lost_when_present = loss_when_present_p_value < q_val,
      gain_when_present_q_val = p.adjust(gain_when_present_q_val, method = "bonferroni"),
      gain_when_present_q_val = ifelse(is.na(gain_when_present_q_val), 1, gain_when_present_q_val),
      is_gained_when_present = gain_when_present_p_value < q_val
    )
  perturbs <- perturbs %>% tidyr::nest(perturb_summary_tbl = !perturb_name)

  perturbation_ccm_tbl$perturb_summary_tbl <- NULL
  perturbation_ccm_tbl <- left_join(perturbation_ccm_tbl, perturbs, by = "perturb_name")
  # perturbation_ccm_tbl$perturb_summary_tbl = perturbs$perturb_summary_tbl

  return(perturbation_ccm_tbl)
}

#' Assemble a state transition graph from a set of perturbations
#'
#' This function takes as input a control timeseries Hooke model, and a set of timeseries perturbation models.
#' Each perturbation model describes a separate experimental perturbation that may eliminate one or more
#' cell states in the experiment. The tibble must have columns "perturb_name" and "perturb_ccm", and each row
#' must have a perturbation model with a unique name.
#'
#' The function returns a state graph with edges annotated by the level of support from the perturbations.
#' @param discordant_config Optional list with keys `power_threshold`,
#'   `prune_mode` (`"existing"` or `"greedy"`), `K_paths`, and `lambda_edge`.
#'   Defaults preserve existing behavior (`prune_mode = "existing"`).
#' @export
assemble_transition_graph_from_perturbations <- function(ref_ccs,
                                                         timeseries_graph,
                                                         perturbation_ccm_tbl,
                                                         q_val = 0.01,
                                                         start_time = NULL,
                                                         stop_time = NULL,
                                                         # perturbation_col="knockout",
                                                         interval_col = "timepoint",
                                                         interval_step = 2,
                                                         min_interval = 4,
                                                         max_interval = 24,
                                                         log_abund_detection_thresh = -5,
                                                         min_pathfinding_lfc = 0,
                                                         links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
                                                         components = "partition",
                                                         verbose = FALSE,
                                                         edge_allowlist = NULL,
                                                         edge_denylist = NULL,
                                                         discordant_config = NULL,
                                                         discordant_pruning_mode = c("none", "greedy"),
                                                         discordant_pruning_k = 1,
                                                         discordant_pruning_power_threshold = 0,
                                                         discordant_pruning_cost_attr = "total_path_score_supporting",
                                                         discordant_pruning_lambda_edge = 0,
                                                         newdata = tibble()) {
  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  # old_omp_num_threads = single_thread_omp()
  # old_blas_num_threads = single_thread_blas()

  tryCatch({
    discordant_pruning_mode <- match.arg(discordant_pruning_mode)
    discordant_cfg <- normalize_discordant_pruning_config(
      discordant_config = discordant_config,
      discordant_pruning_mode = discordant_pruning_mode,
      discordant_pruning_k = discordant_pruning_k,
      discordant_pruning_power_threshold = discordant_pruning_power_threshold,
      discordant_pruning_lambda_edge = discordant_pruning_lambda_edge
    )

    # Get a table of the cell types that are in the control
    # FIXME: "knockout" is hard coded and should be a user-defined term in the model

    if (nrow(newdata) > 0) {
      newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
    } else {
      newdata_wt <- tibble(knockout = FALSE)
    }

    if (is.null(perturbation_ccm_tbl$perturb_summary_tbl)) {
      perturbation_ccm_tbl <- assess_perturbation_effects(perturbation_ccm_tbl,
        q_val = q_val,
        start_time = start_time,
        stop_time = stop_time,
        interval_col = interval_col,
        interval_step = interval_step,
        min_interval = min_interval,
        max_interval = max_interval,
        log_abund_detection_thresh = log_abund_detection_thresh,
        min_lfc = min_pathfinding_lfc,
        verbose = verbose,
        newdata = newdata
      )
    }

    # pathfinding_graph = igraph::intersection(pathfinding_graph, timeseries_graph)
    pathfinding_graph <- timeseries_graph

    if (verbose) {
      message("Computing paths between cell state losses")
    }

    # Now let's find paths between nodes that are lost following each perturbation. This step
    # returns a tibble of paths, each of which links lost nodes in the pathfinding graph.
    # Paths are also consistent with the flow of time in the control, at least over the window
    # of time measured in the corresponding perturbation experiment. Note that the same path may
    # recur in multiple perturbation experiments. Such paths are returned separately in this table.
    path_tbl <- perturbation_ccm_tbl %>%
      dplyr::mutate(paths_between_concordant_loss_nodes = purrr::map2(
        .f = purrr::possibly(
          get_perturbation_paths, NA_real_
        ),
        .x = perturb_ccm,
        .y = perturb_summary_tbl,
        pathfinding_graph
      ))

    path_tbl <- path_tbl %>%
      filter(!is.na(paths_between_concordant_loss_nodes)) %>%
      select(-perturb_summary_tbl)


    if (nrow(path_tbl) == 0) {
      stop("No loss paths")
    } else {
      path_tbl <- path_tbl %>%
        tidyr::unnest(paths_between_concordant_loss_nodes) %>%
        filter(!is.na(path))
      message(paste("Found ", nrow(path_tbl), " loss paths"))
    }


    cells_along_path_df <- normalized_counts(ref_ccs, "size_only", pseudocount = 0) %>%
      as.matrix() %>%
      Matrix::t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(!matches("rowname")) %>%
      dplyr::rename(sample = rowname, cell_group = name, num_cells = value)

    if (verbose) {
      message(paste("Assessing time flows across", nrow(path_tbl), " loss paths"))
    }

    # Measure the flow of time along each perturbation loss path
    path_tbl <- path_tbl %>%
      mutate(time_vs_distance_model_stats = furrr::future_map(
        .f = purrr::possibly(measure_time_delta_along_path, NA_character_),
        .x = path,
        ccs = ref_ccs,
        cells_along_path_df = cells_along_path_df,
        interval_col = interval_col,
        .options = furrr::furrr_options(stdout = FALSE, conditions = character()),
        .progress = TRUE
      )) %>%
      tidyr::unnest(time_vs_distance_model_stats)

    # Compute some scores that summarize the level of support for each path by the various perturbations.
    path_tbl <- path_tbl %>% mutate(time_dist_effect_qval = p.adjust(time_dist_effect_pval, method = "BH"))
    path_tbl <- path_tbl %>% mutate(timeseries_model_path_score = ifelse(time_dist_effect > 0 & time_dist_effect_qval < q_val, time_dist_model_adj_rsq, 0))
    path_tbl <- path_tbl %>% mutate(perturb_model_path_score = ifelse(perturb_dist_effect < 0 & perturb_dist_effect_qval < q_val, abs(perturb_dist_effect), 0))
    path_tbl <- path_tbl %>% mutate(path_score = ifelse(timeseries_model_path_score > 0 & perturb_model_path_score > 0, perturb_model_path_score * timeseries_model_path_score, 0))

    # path_tbl = path_tbl %>% mutate(path_score = ifelse(timeseries_model_path_score > 0, timeseries_model_path_score, 0))

    # Exclude paths that have no support from perturbations or go against the flow of time
    path_tbl <- path_tbl %>% filter(path_score > 0)

    message(paste("Found ", nrow(path_tbl), "significant loss paths with forward time-flow"))

    if (nrow(path_tbl) == 0) {
      stop("No significant loss paths")
    }


    if (verbose) {
      message("Constructing transition graph")
    }

    # Now let's start to build up a state transition graph from the paths by taking their union across all
    # perturbations, provided they survived the above filtering steps.
    selected_paths <- path_tbl %>%
      # group_by(to) %>%
      ungroup() %>%
      arrange(desc(path_score)) %>%
      dplyr::rename(path_contrast = perturb_name)
    G <- select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = TRUE)

    # Now go back and score each edge in the state graph for support from perturbations, as well as support
    # the full control timeseries
    if (verbose) {
      message("Assessing perturbation support for transition graph")
    }
    G <- assess_support_for_transition_graph(
      perturbation_ccm_tbl,
      path_tbl,
      G,
      q_val,
      start_time,
      stop_time,
      # perturbation_col,
      interval_col,
      interval_step,
      min_interval,
      max_interval,
      log_abund_detection_thresh
    )

    if (discordant_cfg$prune_mode == "greedy") {
      discordant_pairs <- collect_discordant_pairs_from_perturbation_summaries(
        perturbation_ccm_tbl = perturbation_ccm_tbl,
        power_threshold = discordant_cfg$power_threshold
      )
      pruning_res <- prune_discordant_paths_greedily(
        state_graph = G,
        discordant_pairs = discordant_pairs,
        k_paths = discordant_cfg$K_paths,
        deletion_cost_attr = discordant_pruning_cost_attr,
        traversal_weight_attr = "weight",
        lambda_edge = discordant_cfg$lambda_edge
      )
      G <- pruning_res$graph
      igraph::graph_attr(G, "discordant_pruning_removed_edges") <- list(pruning_res$removed_edges)
      igraph::graph_attr(G, "discordant_pruning_mode") <- discordant_cfg$prune_mode
    }


    # FIXME: this is gross and there is probably a cleaner way, but what we're doing here
    # is annotating that these edges from the timeseries are not supported by the perturbations
    G_num_perturbs_supporting <- igraph::edge_attr(G, "num_perturbs_supporting")
    G_num_perturbs_supporting[is.na(G_num_perturbs_supporting)] <- 0
    igraph::edge_attr(G, "num_perturbs_supporting") <- G_num_perturbs_supporting

    G_max_timeseries_path_score_supporting <- igraph::edge_attr(G, "max_timeseries_path_score_supporting")
    G_max_timeseries_path_score_supporting[is.na(G_max_timeseries_path_score_supporting)] <- 0
    igraph::edge_attr(G, "max_timeseries_path_score_supporting") <- G_max_timeseries_path_score_supporting

    G_total_timeseries_path_score_supporting <- igraph::edge_attr(G, "total_timeseries_path_score_supporting")
    G_total_timeseries_path_score_supporting[is.na(G_total_timeseries_path_score_supporting)] <- 0
    igraph::edge_attr(G, "total_timeseries_path_score_supporting") <- G_total_timeseries_path_score_supporting

    G_max_perturb_path_score_supporting <- igraph::edge_attr(G, "max_perturb_path_score_supporting")
    G_max_perturb_path_score_supporting[is.na(G_max_perturb_path_score_supporting)] <- 0
    igraph::edge_attr(G, "max_perturb_path_score_supporting") <- G_max_perturb_path_score_supporting

    G_total_perturb_path_score_supporting <- igraph::edge_attr(G, "total_perturb_path_score_supporting")
    G_total_perturb_path_score_supporting[is.na(G_total_perturb_path_score_supporting)] <- 0
    igraph::edge_attr(G, "total_perturb_path_score_supporting") <- G_total_perturb_path_score_supporting

    G_max_path_score_supporting <- igraph::edge_attr(G, "max_path_score_supporting")
    G_max_path_score_supporting[is.na(G_max_path_score_supporting)] <- 0
    igraph::edge_attr(G, "max_path_score_supporting") <- G_max_path_score_supporting

    G_total_path_score_supporting <- igraph::edge_attr(G, "total_path_score_supporting")
    G_total_path_score_supporting[is.na(G_total_path_score_supporting)] <- 0
    igraph::edge_attr(G, "total_path_score_supporting") <- G_total_path_score_supporting

    G_label <- igraph::edge_attr(G, "support_label")
    G_label[is.na(G_label)] <- ""
    igraph::edge_attr(G, "support_label") <- G_label

    return(G)
  }, error = function(e) {
    print(e)
    return(NA)
  }, finally = {
    # RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    # RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  })

  return(NA)
}

#' Assess support for a graph built via perturbations
#' @export
assess_support_for_transition_graph <- function(perturbation_ccm_tbl,
                                                perturbation_path_tbl,
                                                state_transition_graph,
                                                q_val = 0.01,
                                                start_time = NULL,
                                                stop_time = NULL,
                                                # perturbation_col="knockout",
                                                interval_col = "timepoint",
                                                interval_step = 2,
                                                min_interval = 4,
                                                max_interval = 24,
                                                log_abund_detection_thresh = -5) {
  # Flatten all the paths that the perturbation assembler used to link up the
  # cell states
  path_score_tbl <- perturbation_path_tbl %>%
    dplyr::rename(origin = from, destination = to) %>%
    tidyr::unnest(path)

  # print (path_score_tbl)

  # edge_support = left_join(edge_support, path_score_tbl)

  # Uplift some statistics collected during assembly to score edges based on
  # strength of the supporting evidence
  edge_support_summary <- path_score_tbl %>%
    ungroup() %>%
    dplyr::select(
      from, to, perturb_name,
      timeseries_model_path_score,
      perturb_model_path_score,
      path_score
    ) %>%
    group_by(from, to, perturb_name) %>%
    summarize(
      timeseries_model_path_score = max(timeseries_model_path_score),
      perturb_model_path_score = max(perturb_model_path_score),
      path_score = max(path_score)
    ) %>%
    distinct() %>%
    ungroup() %>%
    group_by(from, to) %>%
    summarize(
      num_perturbs_supporting = length(unique(perturb_name)),
      # num_perturb_intervals_supporting=sum(num_intervals_supported),
      max_timeseries_path_score_supporting = max(timeseries_model_path_score),
      total_timeseries_path_score_supporting = sum(timeseries_model_path_score),
      max_perturb_path_score_supporting = max(perturb_model_path_score),
      total_perturb_path_score_supporting = sum(perturb_model_path_score),
      max_path_score_supporting = max(path_score),
      total_path_score_supporting = sum(path_score)
    )

  # Which perturbations support each edge?
  edge_support_perturbs <- path_score_tbl %>%
    ungroup() %>%
    arrange(desc(perturb_model_path_score)) %>%
    dplyr::select(from, to, perturb_name) %>%
    group_by(from, to) %>%
    distinct() %>%
    tidyr::nest(supporting_perturbs = c(perturb_name))

  # Which perturbations lead to a direct loss of each node?
  # Direct loss is defined as either:
  # - A node that is lost at its peak in the WT (and where the perturbation time series covers that peak) and NOT part of a loss path
  # - The first node in a loss path (both the first and last nodes in a loss path are lost at peak WT abundance)
  # FIXME: right now the code below does not detect the first category of direct losses. Fixing this
  # will require that we pass in the arguments needed to call estimate_loss_timing()
  node_direct_perturbs <- path_score_tbl %>%
    ungroup() %>%
    arrange(desc(perturb_model_path_score)) %>%
    dplyr::select(id = origin, perturb_name) %>%
    distinct()

  # Which perturbations lead to an indirect loss of each node?
  # Indirect loss is all nodes on loss paths except the first.

  node_indirect_perturbs <- path_score_tbl %>%
    group_by(origin, destination) %>%
    ungroup() %>%
    dplyr::select(origin, id = from, perturb_name)
  node_indirect_perturbs <- node_indirect_perturbs %>% rbind(
    path_score_tbl %>%
      ungroup() %>%
      dplyr::select(origin, id = to, perturb_name)
  )
  node_indirect_perturbs <- node_indirect_perturbs %>%
    filter(origin != id) %>%
    select(id, perturb_name)
  node_indirect_perturbs <- node_indirect_perturbs %>% distinct()

  # If a node is the first in a path, but internal to another, count it as an indirect loss:
  node_direct_perturbs <- dplyr::setdiff(node_direct_perturbs, node_indirect_perturbs)

  node_indirect_perturbs <- node_indirect_perturbs %>%
    group_by(id) %>%
    tidyr::nest(indirect_perturb = perturb_name)

  node_direct_perturbs <- node_direct_perturbs %>%
    group_by(id) %>%
    tidyr::nest(direct_perturb = perturb_name)

  node_metadata <- tibble(id = igraph::V(state_transition_graph)$name)
  node_metadata <- node_metadata %>% left_join(node_direct_perturbs)
  node_metadata <- node_metadata %>% left_join(node_indirect_perturbs)

  # TODO: we should also annotate nodes as having "inferred" dependency on a given
  # perturbation if its all its ancestors depend on that perturbation, but the node's
  # peak abundance is outside the window of our collected measurements for that
  # perturbation

  # Construct nice labels for labeling edges according to support
  # FIXME: we should probably just move this to its own function to be
  # used in plotting functions
  edge_support_labels <- path_score_tbl %>%
    ungroup() %>%
    arrange(desc(perturb_model_path_score)) %>%
    dplyr::select(from, to, perturb_name) %>%
    group_by(from, to) %>%
    distinct() %>%
    summarize(
      edge_name = stringr::str_c(from, to, sep = "~"),
      # supporting_perturbs = perturb_name,
      support_label = ifelse(n() > 3, paste0(c(perturb_name[1:3], paste("+", n() - 3, " more", sep = "")), collapse = "\n"),
        paste0(perturb_name, collapse = "\n")
      )
    )

  # print (edge_support_labels)
  edge_support_summary <- edge_support_summary %>% left_join(edge_support_labels)
  edge_support_summary <- edge_support_summary %>% left_join(edge_support_perturbs)
  edge_discordance_penalty <- compute_edge_discordance_penalty(
    perturbation_ccm_tbl = perturbation_ccm_tbl,
    state_transition_graph = state_transition_graph
  )
  edge_support_summary <- edge_support_summary %>% left_join(edge_discordance_penalty, by = c("from", "to"))

  # edge_support_summary = edge_support_summary %>% mutate(support_weight = ifelse(is.na(support_weight), 0, support_weight))
  edge_support_summary <- edge_support_summary %>% dplyr::distinct()
  # print ("re-nesting")
  # edge_support = edge_support %>%
  # dplyr::select(-distance_from_root, -weight) %>%
  # tidyr::nest(support=c(perturb_name, perturb_dist_effect, perturb_dist_effect_q_val))
  #  tidyr::nest(support=c(perturb_name))

  # print (edge_support_summary)
  # print ("annotating the graph")
  annotated_state_transition_graph <- igraph::as_data_frame(state_transition_graph) %>%
    as_tibble() %>%
    dplyr::select(from, to) %>%
    distinct() # %>% left_join(edge_support)
  annotated_state_transition_graph <- annotated_state_transition_graph %>% left_join(edge_support_summary)
  # annotated_state_transition_graph = annotated_state_transition_graph %>% mutate(support_weight = ifelse(is.na(support_weight), 0, support_weight))
  annotated_state_transition_graph <- igraph::graph_from_data_frame(annotated_state_transition_graph, directed = TRUE, vertices = node_metadata)
  return(annotated_state_transition_graph)
}

#' Simplify a directed state transition graph by grouping nodes according to a specified label
#' assumes the graphs over nodes corresponding to groups in the ccm
#' @export
contract_state_graph <- function(ccs,
                                 state_graph,
                                 group_nodes_by,
                                 edge_attr_policy = list(
                                   "weight" = "sum",
                                   "name" = "concat",
                                   "num_perturbs_supporting" = "sum",
                                   "max_timeseries_path_score_supporting" = "sum",
                                   "total_timeseries_path_score_supporting" = "sum",
                                   "max_perturb_path_score_supporting" = "sum",
                                   "total_perturb_path_score_supporting" = "sum",
                                   "max_path_score_supporting" = "sum",
                                   "total_path_score_supporting" = "sum",
                                   "edge_name" = "ignore",
                                   "support_label" = "concat",
                                   "supporting_perturbs" = "concat",
                                   "ignore"
                                 )) {
  # Create simplified cell state graph just on cell type (not cluster):
  cell_groups <- ccs@metadata[["cell_group_assignments"]] %>%
    pull(cell_group) %>%
    unique()
  node_metadata <- tibble(id = cell_groups)

  # G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata <- colData(ccs@cds) %>%
    as.data.frame() %>%
    select(!!sym(group_nodes_by))

  cell_group_metadata$cell_group <- ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  group_by_metadata <- cell_group_metadata[, c("cell_group", group_nodes_by)] %>%
    as.data.frame() %>%
    dplyr::count(cell_group, !!sym(group_nodes_by)) %>%
    dplyr::group_by(cell_group) %>%
    slice_max(n, with_ties = FALSE) %>%
    dplyr::select(-n)
  colnames(group_by_metadata) <- c("cell_group", "group_nodes_by")

  node_metadata <- igraph::as_data_frame(state_graph, what = "vertices") %>%
    mutate(order = row_number()) %>%
    left_join(group_by_metadata, by = c("name" = "cell_group"))

  # node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  # node_metadata = left_join(node_metadata, df, by = "id")%>% arrange(sort_id)
  # node_metadata = node_metadata %>% mutate(sort_id = as.numeric(gsub("\\D", "", id))) %>% arrange(sort_id)
  contraction_mapping <- as.factor(node_metadata$group_nodes_by)
  contraction_mapping_names <- as.character(levels(contraction_mapping))
  contraction_mapping <- as.numeric(contraction_mapping)
  names(contraction_mapping) <- node_metadata$id
  contracted_state_graph <- igraph::contract(state_graph, mapping = contraction_mapping, vertex.attr.comb = "ignore")
  igraph::V(contracted_state_graph)$name <- unlist(contraction_mapping_names[as.numeric(igraph::V(contracted_state_graph))])
  contracted_state_graph <- igraph::simplify(contracted_state_graph,
    edge.attr.comb = edge_attr_policy
  )
  igraph::V(contracted_state_graph)$cell_group <- group_nodes_by
  return(contracted_state_graph)
}
