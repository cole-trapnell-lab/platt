#' Plot Cell Type Control Kinetics
#'
#' This function generates a kinetic plot of cell type control data over time.
#'
#' @param control_ccm A control cell count matrix object.
#' @param cell_groups A vector of cell groups to include in the plot. Default is NULL.
#' @param start_time The start time for the plot. Default is NULL, which uses the minimum timepoint in the data.
#' @param stop_time The stop time for the plot. Default is NULL, which uses the maximum timepoint in the data.
#' @param interval_step The step size for the time intervals. Default is 1.
#' @param log_abund_detection_thresh The log abundance detection threshold. Default is -3.
#' @param delta_log_abund_loss_thresh The delta log abundance loss threshold. Default is 0.
#' @param interval_col The column name for the time intervals. Default is "timepoint".
#' @param q_val The q-value threshold for significance. Default is 0.01.
#' @param min_log_abund The minimum log abundance. Default is -5.
#' @param log_scale A boolean indicating whether to use a log scale for the y-axis. Default is TRUE.
#' @param group_nodes_by The column name to group nodes by. Default is "cell_type".
#' @param nrow The number of rows in the facet wrap. Default is 1.
#' @param newdata A tibble containing new data to be used in the plot. Default is an empty tibble.
#' @param color_points_by The column name to color points by. Default is NULL.
#' @param size The size of the points in the plot. Default is 0.5.
#' @param alpha The alpha transparency of the points in the plot. Default is 0.5.
#' @param raw_counts A boolean indicating whether to use raw counts. Default is FALSE.
#'
#' @return A ggplot object representing the kinetic plot.
#'
#' @examples
#' # Example usage:
#' plot <- plot_cell_type_control_kinetics(control_ccm)
#' print(plot)
#'
#' @import assertthat
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @importFrom Matrix t
#' @importFrom purrr map
#' @importFrom tibble tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join select mutate filter group_by slice_max slice_min pull
#' @importFrom ggplot2 aes geom_point geom_line facet_wrap scale_y_log10 ylab scale_color_manual
#'
#' @export
plot_cell_type_control_kinetics <- function(control_ccm,
                                            cell_groups = NULL,
                                            start_time = NULL,
                                            stop_time = NULL,
                                            interval_step = 1,
                                            log_abund_detection_thresh = -3,
                                            delta_log_abund_loss_thresh = 0,
                                            interval_col = "timepoint",
                                            q_val = 0.01,
                                            min_log_abund = -5,
                                            log_scale = TRUE,
                                            group_nodes_by = "cell_type",
                                            nrow = 1,
                                            newdata = tidyr::tibble(),
                                            color_points_by = NULL,
                                            size = 0.5,
                                            alpha = 0.5,
                                            linewidth= 1,
                                            raw_counts = FALSE) {
  # assertthat::assert_that(nrow(newdata) == 1)

  colData(control_ccm@ccs)[, interval_col] <- as.numeric(colData(control_ccm@ccs)[, interval_col])

  if (is.null(start_time)) {
    start_time <- min(colData(control_ccm@ccs)[, interval_col])
  }
  if (is.null(stop_time)) {
    stop_time <- max(colData(control_ccm@ccs)[, interval_col])
  }

  if (nrow(newdata) > 0) {
    newdata <- dplyr::cross_join(tibble(knockout = FALSE), newdata)
  } else {
    newdata <- tidyr::tibble(knockout = FALSE)
  }


  wt_timepoint_pred_df <- estimate_abundances_over_interval(control_ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    interval_step = interval_step,
    newdata = newdata
  )


  # if (is.null(batch_col)){
  #   wt_timepoint_pred_df = estimate_abundances_over_interval(control_ccm, start_time, stop_time,
  #                                                            interval_col=interval_col, interval_step=interval_step,
  #                                                            newdata = newdata)
  # } else if (is.null(reference_batch)) {
  #
  #   batches = tibble(batch = unique(colData(control_ccm@ccs)[,batch_col]))
  #
  #   batches = batches %>% mutate(tp_preds = purrr::map(.f = function(expt) {
  #     newdata[[batch_col]] = expt
  #     estimate_abundances_over_interval(control_ccm,
  #                                               start_time,
  #                                               stop_time,
  #                                               interval_col=interval_col,
  #                                               interval_step=interval_step,
  #                                               min_log_abund = min_log_abund,
  #                                               newdata = newdata)
  #   }, .x=batch))
  #   wt_timepoint_pred_df = batches %>% select(tp_preds) %>% tidyr::unnest(tp_preds)
  #
  #   vibrant.colors =
  #     c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
  #       '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')
  #   my_colors = colorRampPalette(c(vibrant.colors))(nrow(batches))
  # } else {
  #
  #   newdata[[batch_col]] = reference_batch
  #   wt_timepoint_pred_df = estimate_abundances_over_interval(control_ccm,
  #                                                            start_time,
  #                                                            stop_time,
  #                                                            interval_col = interval_col,
  #                                                            interval_step=interval_step,
  #                                                            newdata = newdata)
  #   batches = tibble(batch = unique(colData(control_ccm@ccs)[,batch_col]))
  #   vibrant.colors =
  #     c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
  #       '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')
  #   my_colors = colorRampPalette(c(vibrant.colors))(nrow(batches))
  # }

  if (is.null(color_points_by) == FALSE) {
    colors <- tidyr::tibble(batch = unique(colData(control_ccm@ccs)[, color_points_by]))
    vibrant.colors <-
      c(
        "#EE7733", "#0077BB", "#228833", "#33BBEE", "#EE3377", "#CC3311",
        "#AA3377", "#009988", "#004488", "#DDAA33", "#99CC66", "#D590DD"
      )
    my_colors <- colorRampPalette(c(vibrant.colors))(nrow(colors))
  }


  timepoints <- seq(start_time, stop_time, interval_step)

  # print (earliest_loss_tbl)
  extant_wt_tbl <- get_extant_cell_types(control_ccm,
    start_time,
    stop_time,
    log_abund_detection_thresh = log_abund_detection_thresh,
    interval_col = interval_col, 
    newdata = newdata
  )

  cell_group_order <- extant_wt_tbl %>%
    dplyr::group_by(cell_group) %>%
    dplyr::slice_max(percent_max_abund, n = 1) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::slice_min(!!sym(interval_col), n = 1) %>%
    dplyr::pull(cell_group)

  sel_ccs_counts <- monocle3::normalized_counts(control_ccm@ccs, norm_method = "size_only", pseudocount = 0)
  use_latent <- T
  if (raw_counts == FALSE) {
    sample_metadata <- colData(control_ccm@ccs) %>% tidyr::as_tibble()

    # if (is.null(reference_batch) == FALSE)
    # sample_metadata$expt = reference_batch

    # override the columns with the newdata columns
    for (c in colnames(newdata)) {
      sample_metadata[[c]] <- newdata[[c]]
    }

    conditional_counts <- hooke::estimate_abundances_cond(control_ccm,
      newdata = sample_metadata,
      cond_responses = sel_ccs_counts,
      # cond_responses=Matrix::t(sel_ccs_counts),
      pln_model = "reduced"
    )
    sel_ccs_counts_long <- conditional_counts %>%
      dplyr::select(embryo = sample, cell_group, log_abund) %>%
      dplyr::mutate(num_cells = exp(log_abund)) %>%
      dplyr::select(-log_abund)
  } else if (use_latent) {
    latent_counts <- Matrix::t(control_ccm@full_model_family$latent)
    rownames(latent_counts) <- rownames(control_ccm@ccs)
    colnames(latent_counts) <- colnames(control_ccm@ccs)

    sel_ccs_counts_long <- tibble::rownames_to_column(as.matrix(latent_counts) %>%
      as.data.frame(), var = "cell_group") %>%
      tidyr::pivot_longer(!cell_group, names_to = "embryo", values_to = "num_cells")
  } else {
    sel_ccs_counts <- Matrix::t(Matrix::t(counts(control_ccm@ccs)) / exp(control_ccm@model_aux[["full_model_offsets"]]))
    sel_ccs_counts_long <- tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>%
      as.data.frame(), var = "cell_group") %>%
      tidyr::pivot_longer(!cell_group, names_to = "embryo", values_to = "num_cells")
  }

  cell_group_metadata <- collect_psg_node_metadata(control_ccm@ccs,
    group_nodes_by = group_nodes_by,
    color_nodes_by = group_nodes_by,
    label_nodes_by = group_nodes_by
  ) %>%
    dplyr::select(id, cell_group = group_nodes_by)

  sel_ccs_counts_long <- dplyr::left_join(sel_ccs_counts_long,
    cell_group_metadata,
    by = c("cell_group" = "id")
  )

  if (is.null(color_points_by)) {
    sel_ccs_counts_long <- dplyr::left_join(sel_ccs_counts_long,
      colData(control_ccm@ccs) %>% as.data.frame() %>%
        dplyr::select(sample, !!sym(interval_col)),
      by = c("embryo" = "sample")
    )
  } else {
    sel_ccs_counts_long <- dplyr::left_join(sel_ccs_counts_long,
      colData(control_ccm@ccs) %>% as.data.frame() %>%
        dplyr::select(sample, !!sym(interval_col), !!sym(color_points_by)),
      by = c("embryo" = "sample")
    )
  }


  wt_timepoint_pred_df$cell_group <- factor(as.character(wt_timepoint_pred_df$cell_group), levels = cell_group_order)
  sel_ccs_counts_long$cell_group <- factor(as.character(sel_ccs_counts_long$cell_group), levels = cell_group_order)

  if (is.null(cell_groups) == FALSE) {
    wt_timepoint_pred_df <- wt_timepoint_pred_df %>% dplyr::filter(cell_group %in% cell_groups)
    sel_ccs_counts_long <- sel_ccs_counts_long %>% dplyr::filter(cell_group %in% cell_groups)
  }

  wt_timepoint_pred_df[[interval_col]] <- as.numeric(wt_timepoint_pred_df[[interval_col]])
  sel_ccs_counts_long[[interval_col]] <- as.numeric(sel_ccs_counts_long[[interval_col]])

  sel_ccs_counts_long <- sel_ccs_counts_long %>% dplyr::filter(!!sym(interval_col) >= start_time, !!sym(interval_col) <= stop_time)

  kinetic_plot <-
    ggplot(wt_timepoint_pred_df, aes(x = !!sym(interval_col))) +
    geom_point(
      data = sel_ccs_counts_long,
      aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh), color = expt, alpha = alpha),
      position = "jitter", 
      size = size
    ) +
    facet_wrap(~cell_group, scales = "free_y", nrow = nrow) +
    monocle3:::monocle_theme_opts()

  if (is.null(color_points_by)) {
    kinetic_plot <- ggplot(wt_timepoint_pred_df, aes(x = !!sym(interval_col))) +
      geom_point(
        data = sel_ccs_counts_long,
        aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh)), alpha = alpha,
        position = "jitter", size = size
      ) +
      facet_wrap(~cell_group, scales = "free_y", nrow = nrow) +
      monocle3:::monocle_theme_opts()
    kinetic_plot <- kinetic_plot +
      # scale_color_manual(values = my_colors) +
      geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh)), linewidth = linewidth)
  } else {
    
    sel_ccs_counts_long[[color_points_by]] <- as.character(sel_ccs_counts_long[[color_points_by]])
    kinetic_plot <- ggplot(wt_timepoint_pred_df, aes(x = !!sym(interval_col))) +
      geom_point(
        data = sel_ccs_counts_long,
        aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh), color = !!sym(color_points_by)),
        alpha = alpha,
        position = "jitter", size = size
      ) +
      facet_wrap(~cell_group, scales = "free_y", nrow = nrow) +
      monocle3:::monocle_theme_opts()
    
    if (is.null(wt_timepoint_pred_df[[color_points_by]])) {
      kinetic_plot <- kinetic_plot +
        scale_color_manual(values = my_colors) +
        geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh)), linewidth = linewidth)
    } else {
      kinetic_plot <- kinetic_plot +
        scale_color_manual(values = my_colors) +
        geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh), color = !!sym(color_points_by)), linewidth = linewidth)
    }
   
  }

  # kinetic_plot = ggplot(wt_timepoint_pred_df, aes(x = !!sym(interval_col))) +
  #   geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh), color=!!sym(batch_col))) +
  #   facet_wrap(~cell_group, scales="free_y", ncol=1)
  #
  # kinetic_plot = kinetic_plot +
  #   geom_point(data=sel_ccs_counts_long,
  #              aes(x = !!sym(interval_col), y = num_cells +  exp(log_abund_detection_thresh), color=!!sym(batch_col)),
  #              position="jitter") +
  #   # facet_wrap(~cell_group, scales="free_y", ncol=1)
  #   facet_wrap(~cell_group, scales="free_y", nrow=1)

  if (log_scale) {
    kinetic_plot <- kinetic_plot + scale_y_log10()
  }

  kinetic_plot <- kinetic_plot + ylab("cells per sample")

  return(kinetic_plot)
}


#' Plot Cell Type Perturbation Kinetics
#'
#' This function generates a plot to visualize the kinetics of cell type perturbations over time.
#'
#' @param perturbation_ccm A perturbation cell count matrix (CCM) object.
#' @param cell_groups A vector of cell groups to include in the plot. If NULL, all cell groups are included.
#' @param start_time The start time for the interval. If NULL, the minimum timepoint in the data is used.
#' @param stop_time The stop time for the interval. If NULL, the maximum timepoint in the data is used.
#' @param interval_step The step size for the time intervals.
#' @param log_abund_detection_thresh The log abundance detection threshold.
#' @param delta_log_abund_loss_thresh The threshold for the change in log abundance loss.
#' @param interval_col The column name for the time intervals.
#' @param q_val The q-value threshold for significance.
#' @param log_scale A logical value indicating whether to use a log scale for the y-axis.
#' @param control_ccm A control cell count matrix (CCM) object. Defaults to the perturbation CCM.
#' @param control_start_time The start time for the control interval. Defaults to the start time.
#' @param control_stop_time The stop time for the control interval. Defaults to the stop time.
#' @param group_nodes_by The column name to group nodes by.
#' @param newdata A tibble containing new data for predictions.
#' @param nrow The number of rows for the facet wrap.
#' @param size The size of the points in the plot.
#' @param alpha The alpha transparency for the difference shading.
#' @param raw_counts A logical value indicating whether to use raw counts.
#'
#' @return A ggplot object representing the kinetics of cell type perturbations.
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr mutate group_by slice_max slice_min pull select left_join filter
#' @importFrom tidyr unnest pivot_longer
#' @importFrom purrr map
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap scale_y_log10 xlab ylab
#' @importFrom tibble tibble rownames_to_column
#' @importFrom Matrix t
#' @importFrom rlang sym
#'
#' @examples
#' # Example usage:
#' plot <- plot_cell_type_perturb_kinetics(perturbation_ccm)
#' print(plot)
#'
#' @export
plot_cell_type_perturb_kinetics <- function(perturbation_ccm,
                                            cell_groups = NULL,
                                            start_time = NULL,
                                            stop_time = NULL,
                                            interval_step = 1,
                                            log_abund_detection_thresh = -3,
                                            delta_log_abund_loss_thresh = 0,
                                            interval_col = "timepoint",
                                            q_val = 0.01,
                                            log_scale = TRUE,
                                            control_ccm = perturbation_ccm,
                                            control_start_time = start_time,
                                            control_stop_time = control_stop_time,
                                            group_nodes_by = "cell_type",
                                            newdata = tibble(),
                                            nrow = 1,
                                            size = 0.75,
                                            alpha = 0.5,
                                            raw_counts = FALSE) {
  #
  # assertthat::assert_that(nrow(newdata)==1)

  colData(perturbation_ccm@ccs)[, interval_col] <- as.numeric(colData(perturbation_ccm@ccs)[, interval_col])

  if (is.null(start_time)) {
    start_time <- min(colData(perturbation_ccm@ccs)[, interval_col])
  }
  if (is.null(stop_time)) {
    stop_time <- max(colData(perturbation_ccm@ccs)[, interval_col])
  }

  if (nrow(newdata) > 0) {
    newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
    newdata_mt <- cross_join(tibble(knockout = TRUE), newdata)
  } else {
    newdata_wt <- tibble(knockout = FALSE)
    newdata_mt <- tibble(knockout = TRUE)
  }


  # maybe this is actually handled by new data always

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

  # Find the pairs of nodes that are both lost in the perturbation at the same time
  perturb_vs_wt_nodes <- tidyr::tibble(t1 = timepoints) %>%
    mutate(comp_abund = purrr::map(
      .f = platt:::compare_ko_to_wt_at_timepoint,
      .x = t1,
      perturbation_ccm = perturbation_ccm,
      interval_col = interval_col,
      wt_pred_df = wt_timepoint_pred_df,
      ko_pred_df = ko_timepoint_pred_df
    )) %>%
    tidyr::unnest(comp_abund)

  extant_wt_tbl <- get_extant_cell_types(perturbation_ccm,
    start_time,
    stop_time,
    interval_col = interval_col,
    log_abund_detection_thresh = log_abund_detection_thresh,
    newdata = newdata_wt
  )

  cell_group_order <- extant_wt_tbl %>%
    dplyr::group_by(cell_group) %>%
    dplyr::slice_max(percent_max_abund, n = 1) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::slice_min(!!sym(interval_col), n = 1) %>%
    dplyr::pull(cell_group)

  sel_ccs_counts <- monocle3::normalized_counts(perturbation_ccm@ccs, norm_method = "size_only", pseudocount = 0)

  # this is new
  if (raw_counts == FALSE) {
    sample_metadata <- colData(perturbation_ccm@ccs) %>% tidyr::as_tibble()

    # if (is.null(reference_batch) == FALSE)
    #   sample_metadata$expt = reference_batch

    # override the columns with the newdata columns
    for (c in colnames(newdata)) {
      sample_metadata[[c]] <- newdata[[c]]
    }

    conditional_counts <- hooke::estimate_abundances_cond(perturbation_ccm,
      newdata = sample_metadata,
      cond_responses = sel_ccs_counts,
      # cond_responses=Matrix::t(sel_ccs_counts),
      pln_model = "reduced"
    )
    sel_ccs_counts_long <- conditional_counts %>%
      dplyr::select(embryo = sample, cell_group, log_abund) %>%
      mutate(num_cells = exp(log_abund)) %>%
      select(-log_abund)
  } else {
    sel_ccs_counts <- Matrix::t(Matrix::t(counts(perturbation_ccm@ccs)) / exp(perturbation_ccm@model_aux[["full_model_offsets"]]))
    sel_ccs_counts_long <- tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>%
      as.data.frame(), var = "cell_group") %>%
      tidyr::pivot_longer(!cell_group, names_to = "embryo", values_to = "num_cells")
  }

  # this is new

  sel_ccs_counts_long <- tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame(), var = "cell_group") %>%
    tidyr::pivot_longer(!cell_group, names_to = "embryo", values_to = "num_cells")

  cell_group_metadata <- collect_psg_node_metadata(perturbation_ccm@ccs,
    group_nodes_by = group_nodes_by,
    color_nodes_by = group_nodes_by,
    label_nodes_by = group_nodes_by
  ) %>%
    dplyr::select(id, cell_group = group_nodes_by)

  sel_ccs_counts_long <- dplyr::left_join(sel_ccs_counts_long,
    cell_group_metadata,
    by = c("cell_group" = "id")
  )

  sel_ccs_counts_long <- dplyr::left_join(sel_ccs_counts_long,
    colData(perturbation_ccm@ccs) %>% as.data.frame() %>% dplyr::select(sample, !!sym(interval_col), knockout, expt),
    by = c("embryo" = "sample")
  )

  # INSERT TIME by peak abundances
  perturb_vs_wt_nodes <- dplyr::left_join(perturb_vs_wt_nodes,
    extant_wt_tbl %>% dplyr::select(cell_group, !!sym(interval_col), present_above_thresh),
    by = c("cell_group" = "cell_group", "t1" = interval_col)
  )

  # if (is.null(cell_groups)){
  #   cell_groups = unique(as.character(perturb_vs_wt_nodes$cell_group)) %>% sort()
  # }
  perturb_vs_wt_nodes$cell_group <- factor(as.character(perturb_vs_wt_nodes$cell_group), levels = cell_group_order)
  sel_ccs_counts_long$cell_group <- factor(as.character(sel_ccs_counts_long$cell_group), levels = cell_group_order)

  if (is.null(cell_groups) == FALSE) {
    sel_ccs_counts_long <- sel_ccs_counts_long %>% filter(cell_group %in% cell_groups)
    perturb_vs_wt_nodes <- perturb_vs_wt_nodes %>% filter(cell_group %in% cell_groups)
  }


  perturb_vs_wt_nodes$t1 <- as.numeric(perturb_vs_wt_nodes$t1)
  sel_ccs_counts_long[[interval_col]] <- as.numeric(sel_ccs_counts_long[[interval_col]])

  kinetic_plot <- ggplot(perturb_vs_wt_nodes, aes(x = !!sym(paste(interval_col, "_x", sep = "")))) +
    geom_point(
      data = perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund > 0 & delta_q_value < q_val),
      aes(
        x = !!sym(paste(interval_col, "_x", sep = "")),
        y = 1.0 * (exp(log_abund_y) + exp(log_abund_detection_thresh))
      ), shape = 8
    ) +
    geom_point(
      data = perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund < 0 & delta_q_value < q_val),
      aes(
        x = !!sym(paste(interval_col, "_x", sep = "")),
        y = 1.0 * (exp(log_abund_y) + exp(log_abund_detection_thresh))
      ), shape = 8
    ) +
    geom_point(
      data = sel_ccs_counts_long %>% filter(knockout == T),
      aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh), shape = knockout),
      color = "black",
      position = "jitter",
      size = size
    ) +
    geom_point(
      data = sel_ccs_counts_long %>% filter(knockout == F),
      aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh), shape = knockout),
      color = "gray",
      position = "jitter",
      size = size
    ) +
    geom_line(aes(y = exp(log_abund_x) + exp(log_abund_detection_thresh), linetype = "Wild-type")) +
    geom_line(aes(y = exp(log_abund_y) + exp(log_abund_detection_thresh), linetype = "Knockout")) +
    ggh4x::stat_difference(aes(
      ymin = exp(log_abund_x) + exp(log_abund_detection_thresh),
      ymax = exp(log_abund_y) + exp(log_abund_detection_thresh)
    ), alpha = alpha) +
    scale_fill_manual(values = c("orangered3", "royalblue3", "lightgray")) +
    facet_wrap(~cell_group, scales = "free_y", nrow = nrow) +
    monocle3:::monocle_theme_opts() +
    geom_hline(yintercept = exp(log_abund_detection_thresh), color = "lightgrey") +
    xlab(interval_col)

  # kinetic_plot = ggplot(perturb_vs_wt_nodes, aes(x = !!sym(paste(interval_col, "_x", sep="")))) +
  #   geom_line(aes(y = exp(log_abund_x) +exp(log_abund_detection_thresh), linetype = "Wild-type")) +
  #   geom_line(aes(y = exp(log_abund_y) +exp(log_abund_detection_thresh), linetype = "Knockout")) +
  # ggh4x::stat_difference(aes(ymin = exp(log_abund_x)+exp(log_abund_detection_thresh), ymax = exp(log_abund_y) +exp(log_abund_detection_thresh)), alpha=0.3) +
  # geom_point(aes(x = !!sym(paste(interval_col, "_x", sep="")), y =1.0*(exp(log_abund_y) + exp(log_abund_detection_thresh))), shape=8, data=perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund > 0 & delta_q_value < q_val)) +
  # geom_point(aes(x = !!sym(paste(interval_col, "_x", sep="")), y =1.0*(exp(log_abund_y) + exp(log_abund_detection_thresh))), shape=8, data=perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund < 0 & delta_q_value < q_val)) +
  #   # facet_wrap(~cell_group, scales="free_y", ncol=1)
  #   facet_wrap(~cell_group, scales="free_y", nrow=1)
  #
  # kinetic_plot = kinetic_plot +
  #   geom_jitter(data=sel_ccs_counts_long,
  #               aes(x = !!sym(interval_col), y = num_cells+exp(log_abund_detection_thresh), color=gene_target, shape=expt),
  #               height=0,
  #               width=2)
  # kinetic_plot = kinetic_plot + geom_hline(yintercept=exp(log_abund_detection_thresh), color="lightgrey")

  if (log_scale) {
    kinetic_plot <- kinetic_plot + scale_y_log10()
  }

  kinetic_plot <- kinetic_plot + ylab("cells per sample")


  return(kinetic_plot)
}




# get information for kinetic plots

get_cell_type_perturb_kinetics_data <- function(perturbation_ccm,
                                                perturb_time_window,
                                                cell_groups = NULL,
                                                # start_time=NULL,
                                                # stop_time=NULL,
                                                interval_step = 1,
                                                log_abund_detection_thresh = -3,
                                                delta_log_abund_loss_thresh = 0,
                                                interval_col = "timepoint",
                                                q_val = 0.01,
                                                log_scale = TRUE,
                                                control_ccm = perturbation_ccm,
                                                control_start_time = start_time,
                                                control_stop_time = control_stop_time,
                                                group_nodes_by = "cell_type",
                                                newdata = tibble(),
                                                nrow = 1,
                                                size = 0.75,
                                                alpha = 0.5,
                                                raw_counts = FALSE) {
  start_time <- as.numeric(perturb_time_window$start_time)
  stop_time <- as.numeric(perturb_time_window$stop_time)
  # assertthat::assert_that(nrow(newdata)==1)

  colData(perturbation_ccm@ccs)[, interval_col] <- as.numeric(colData(perturbation_ccm@ccs)[, interval_col])

  if (is.null(start_time)) {
    start_time <- min(colData(perturbation_ccm@ccs)[, interval_col])
  }
  if (is.null(stop_time)) {
    stop_time <- max(colData(perturbation_ccm@ccs)[, interval_col])
  }

  if (nrow(newdata) > 0) {
    newdata_wt <- cross_join(tibble(knockout = FALSE), newdata)
    newdata_mt <- cross_join(tibble(knockout = TRUE), newdata)
  } else {
    newdata_wt <- tibble(knockout = FALSE)
    newdata_mt <- tibble(knockout = TRUE)
  }


  # maybe this is actually handled by new data always

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

  # Find the pairs of nodes that are both lost in the perturbation at the same time
  perturb_vs_wt_nodes <- tidyr::tibble(t1 = timepoints) %>%
    mutate(comp_abund = purrr::map(
      .f = platt:::compare_ko_to_wt_at_timepoint,
      .x = t1,
      perturbation_ccm = perturbation_ccm,
      interval_col = interval_col,
      wt_pred_df = wt_timepoint_pred_df,
      ko_pred_df = ko_timepoint_pred_df
    )) %>%
    tidyr::unnest(comp_abund)

  extant_wt_tbl <- get_extant_cell_types(perturbation_ccm,
    start_time,
    stop_time,
    log_abund_detection_thresh = log_abund_detection_thresh,
    newdata = newdata_wt
  )

  cell_group_order <- extant_wt_tbl %>%
    group_by(cell_group) %>%
    slice_max(percent_max_abund, n = 1) %>%
    group_by(cell_group) %>%
    slice_min(timepoint, n = 1) %>%
    pull(cell_group)

  sel_ccs_counts <- normalized_counts(perturbation_ccm@ccs, norm_method = "size_only", pseudocount = 0)

  # this is new
  if (raw_counts == FALSE) {
    sample_metadata <- colData(perturbation_ccm@ccs) %>% as_tibble()

    # if (is.null(reference_batch) == FALSE)
    #   sample_metadata$expt = reference_batch

    # override the columns with the newdata columns
    for (c in colnames(newdata)) {
      sample_metadata[[c]] <- newdata[[c]]
    }

    conditional_counts <- estimate_abundances_cond(perturbation_ccm,
      newdata = sample_metadata,
      cond_responses = sel_ccs_counts,
      # cond_responses=Matrix::t(sel_ccs_counts),
      pln_model = "reduced"
    )
    sel_ccs_counts_long <- conditional_counts %>%
      select(embryo = sample, cell_group, log_abund) %>%
      mutate(num_cells = exp(log_abund)) %>%
      select(-log_abund)
  } else {
    sel_ccs_counts <- Matrix::t(Matrix::t(counts(perturbation_ccm@ccs)) / exp(perturbation_ccm@model_aux[["full_model_offsets"]]))
    sel_ccs_counts_long <- tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>%
      as.data.frame(), var = "cell_group") %>%
      pivot_longer(!cell_group, names_to = "embryo", values_to = "num_cells")
  }

  # this is new

  sel_ccs_counts_long <- tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame(), var = "cell_group") %>%
    pivot_longer(!cell_group, names_to = "embryo", values_to = "num_cells")

  cell_group_metadata <- collect_psg_node_metadata(perturbation_ccm@ccs,
    group_nodes_by = group_nodes_by,
    color_nodes_by = group_nodes_by,
    label_nodes_by = group_nodes_by
  ) %>%
    select(id, cell_group = group_nodes_by)

  sel_ccs_counts_long <- left_join(sel_ccs_counts_long,
    cell_group_metadata,
    by = c("cell_group" = "id")
  )

  sel_ccs_counts_long <- left_join(sel_ccs_counts_long,
    colData(perturbation_ccm@ccs) %>% as.data.frame() %>% select(sample, !!sym(interval_col), knockout, expt),
    by = c("embryo" = "sample")
  )

  # INSERT TIME by peak abundances
  perturb_vs_wt_nodes <- left_join(perturb_vs_wt_nodes,
    extant_wt_tbl %>% select(cell_group, !!sym(interval_col), present_above_thresh),
    by = c("cell_group" = "cell_group", "t1" = "timepoint")
  )

  # if (is.null(cell_groups)){
  #   cell_groups = unique(as.character(perturb_vs_wt_nodes$cell_group)) %>% sort()
  # }
  perturb_vs_wt_nodes$cell_group <- factor(as.character(perturb_vs_wt_nodes$cell_group), levels = cell_group_order)
  sel_ccs_counts_long$cell_group <- factor(as.character(sel_ccs_counts_long$cell_group), levels = cell_group_order)

  if (is.null(cell_groups) == FALSE) {
    sel_ccs_counts_long <- sel_ccs_counts_long %>% filter(cell_group %in% cell_groups)
    perturb_vs_wt_nodes <- perturb_vs_wt_nodes %>% filter(cell_group %in% cell_groups)
  }


  perturb_vs_wt_nodes$t1 <- as.numeric(perturb_vs_wt_nodes$t1)
  sel_ccs_counts_long$timepoint <- as.numeric(sel_ccs_counts_long$timepoint)


  res <- list(
    perturb_vs_wt_nodes = perturb_vs_wt_nodes,
    sel_ccs_counts_long = sel_ccs_counts_long
  )

  return(res)
}


plot_cell_type_perturb_kinetics_from_res <- function(perturb_vs_wt_nodes,
                                                     sel_ccs_counts_long,
                                                     cell_groups = NULL,
                                                     interval_step = 1,
                                                     log_abund_detection_thresh = -3,
                                                     delta_log_abund_loss_thresh = 0,
                                                     interval_col = "timepoint",
                                                     q_val = 0.01,
                                                     log_scale = TRUE,
                                                     control_ccm = perturbation_ccm,
                                                     control_start_time = start_time,
                                                     control_stop_time = control_stop_time,
                                                     # group_nodes_by = "cell_type",
                                                     nrow = 1,
                                                     size = 0.75,
                                                     alpha = 0.5,
                                                     raw_counts = FALSE) {
  if (is.null(cell_groups) == FALSE) {
    perturb_vs_wt_nodes <- perturb_vs_wt_nodes %>% filter(cell_group %in% cell_groups)
    sel_ccs_counts_long <- sel_ccs_counts_long %>% filter(cell_group %in% cell_groups)
  }

  kinetic_plot <- ggplot(perturb_vs_wt_nodes, aes(x = !!sym(paste(interval_col, "_x", sep = "")))) +
    geom_point(
      data = perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund > 0 & delta_q_value < q_val),
      aes(
        x = !!sym(paste(interval_col, "_x", sep = "")),
        y = 1.0 * (exp(log_abund_y) + exp(log_abund_detection_thresh))
      ), shape = 8
    ) +
    geom_point(
      data = perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund < 0 & delta_q_value < q_val),
      aes(
        x = !!sym(paste(interval_col, "_x", sep = "")),
        y = 1.0 * (exp(log_abund_y) + exp(log_abund_detection_thresh))
      ), shape = 8
    ) +
    geom_point(
      data = sel_ccs_counts_long %>% filter(knockout == T),
      aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh), shape = knockout),
      color = "black",
      position = "jitter",
      size = size
    ) +
    geom_point(
      data = sel_ccs_counts_long %>% filter(knockout == F),
      aes(x = !!sym(interval_col), y = num_cells + exp(log_abund_detection_thresh), shape = knockout),
      color = "gray",
      position = "jitter",
      size = size
    ) +
    geom_line(aes(y = exp(log_abund_x) + exp(log_abund_detection_thresh), linetype = "Wild-type")) +
    geom_line(aes(y = exp(log_abund_y) + exp(log_abund_detection_thresh), linetype = "Knockout")) +
    ggh4x::stat_difference(aes(
      ymin = exp(log_abund_x) + exp(log_abund_detection_thresh),
      ymax = exp(log_abund_y) + exp(log_abund_detection_thresh)
    ), alpha = alpha) +
    scale_fill_manual(values = c("orangered3", "royalblue3", "lightgray")) +
    facet_wrap(~cell_group, scales = "free_y", nrow = nrow) +
    monocle3:::monocle_theme_opts() +
    geom_hline(yintercept = exp(log_abund_detection_thresh), color = "lightgrey") +
    xlab("timepoint")

  if (log_scale) {
    kinetic_plot <- kinetic_plot + scale_y_log10()
  }

  kinetic_plot <- kinetic_plot + ylab("cells per sample")

  return(kinetic_plot)
}
