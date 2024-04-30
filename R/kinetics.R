#' plots the control 
#' @param control_ccm
#' @param cell_groups
#' @param start_time
#' @param stop_time
#' @param interval_step 
#' @param interval_col 
#' @param batch_col 
#' @param q_val
#' @param newdata
#' @export
plot_cell_type_control_kinetics = function(control_ccm,
                                           cell_groups=NULL,
                                           start_time=NULL,
                                           stop_time=NULL,
                                           interval_step=1,
                                           log_abund_detection_thresh=-3,
                                           delta_log_abund_loss_thresh=0,
                                           interval_col="timepoint",
                                           batch_col=NULL,
                                           q_val=0.01,
                                           min_log_abund = -5,
                                           log_scale=TRUE,
                                           group_nodes_by = "cell_type", 
                                           nrow = 1,
                                           newdata = tibble(), 
                                           reference_batch = NULL, 
                                           raw_counts = FALSE){
  
  
  colData(control_ccm@ccs)[,interval_col] = as.numeric(colData(control_ccm@ccs)[,interval_col])

  if (is.null(start_time))
    start_time = min(colData(control_ccm@ccs)[,interval_col])
  if (is.null(stop_time))
    stop_time = max(colData(control_ccm@ccs)[,interval_col])
  
  if (nrow(newdata) > 0 ){
    newdata = cross_join(tibble(knockout=FALSE) , newdata)
  } else {
    newdata = tibble(knockout=FALSE) 
  }
  
  if (is.null(batch_col)){
    wt_timepoint_pred_df = estimate_abundances_over_interval(control_ccm, start_time, stop_time, 
                                                             interval_col=interval_col, interval_step=interval_step,
                                                             newdata = newdata)
  } else if (is.null(reference_batch)) {

    batches = tibble(batch = unique(colData(control_ccm@ccs)[,batch_col]))
    
    batches = batches %>% mutate(tp_preds = purrr::map(.f = function(expt) {
      newdata[[batch_col]] = expt
      estimate_abundances_over_interval(control_ccm,
                                                start_time,
                                                stop_time,
                                                interval_col=interval_col,
                                                interval_step=interval_step,
                                                min_log_abund = min_log_abund,
                                                newdata = newdata)
    }, .x=batch))
    wt_timepoint_pred_df = batches %>% select(tp_preds) %>% tidyr::unnest(tp_preds)
    
    vibrant.colors =
      c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
        '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')
    my_colors = colorRampPalette(c(vibrant.colors))(nrow(batches))
  } else {
    
    newdata[[batch_col]] = reference_batch
    wt_timepoint_pred_df = estimate_abundances_over_interval(control_ccm, 
                                                             start_time, 
                                                             stop_time, 
                                                             interval_col = interval_col, 
                                                             interval_step=interval_step,
                                                             newdata = newdata)
    batches = tibble(batch = unique(colData(control_ccm@ccs)[,batch_col]))
    vibrant.colors =
      c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
        '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')
    my_colors = colorRampPalette(c(vibrant.colors))(nrow(batches))
  }
  
  
 timepoints = seq(start_time, stop_time, interval_step)

  #print (earliest_loss_tbl)

  peak_wt_abundance = wt_timepoint_pred_df %>% group_by(cell_group) %>% slice_max(log_abund, n=1)

  sel_ccs_counts = normalized_counts(control_ccm@ccs, norm_method="size_only", pseudocount=0)
 
  if (raw_counts == FALSE) {
    sample_metadata = colData(control_ccm@ccs) %>% as_tibble
    
    if (is.null(reference_batch) == FALSE)
      sample_metadata$expt = reference_batch
    
    conditional_counts = estimate_abundances_cond(control_ccm, 
                                                  newdata=sample_metadata, 
                                                  cond_responses= sel_ccs_counts, 
                                                  # cond_responses=Matrix::t(sel_ccs_counts),
                                                  pln_model="reduced")
    sel_ccs_counts_long = conditional_counts %>% 
                          select(embryo=sample, cell_group, log_abund) %>% 
                          mutate(num_cells=exp(log_abund)) %>% select(-log_abund)
    
  } else{
    sel_ccs_counts = Matrix::t(Matrix::t(counts(control_ccm@ccs)) / exp(model.offset(control_ccm@model_aux[["full_model_frame"]])))
    sel_ccs_counts_long = tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% 
                          as.data.frame, var="cell_group") %>%
                          pivot_longer(!cell_group, names_to="embryo", values_to="num_cells")
  }

  cell_group_metadata = collect_psg_node_metadata(control_ccm@ccs,
                                                          group_nodes_by=group_nodes_by,
                                                          color_nodes_by=group_nodes_by,
                                                          label_nodes_by=group_nodes_by) %>%
    select(id, cell_group = group_nodes_by)

  sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                  cell_group_metadata,
                                  by=c("cell_group"="id"))
  
  if (is.null(batch_col)) {
    sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                    colData(control_ccm@ccs) %>% as.data.frame %>%
                                      select(sample, !!sym(interval_col)),
                                    by=c("embryo"="sample"))
  } else {
    sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                    colData(control_ccm@ccs) %>% as.data.frame %>%
                                      select(sample, !!sym(interval_col), !!sym(batch_col)),
                                    by=c("embryo"="sample"))
  }
  

  if (is.null(cell_groups) == FALSE){
    wt_timepoint_pred_df = wt_timepoint_pred_df %>% filter(cell_group %in% cell_groups)
    sel_ccs_counts_long = sel_ccs_counts_long %>% filter(cell_group %in% cell_groups)
  }
  
  if (is.null(cell_groups)){
    cell_groups = unique(as.character(wt_timepoint_pred_df$cell_group)) %>% sort()
  }
  wt_timepoint_pred_df$cell_group = factor(as.character(wt_timepoint_pred_df$cell_group), levels=cell_groups)
  sel_ccs_counts_long$cell_group = factor(as.character(sel_ccs_counts_long$cell_group), levels=cell_groups)
  
  wt_timepoint_pred_df$timepoint = as.numeric(wt_timepoint_pred_df$timepoint)
  sel_ccs_counts_long$timepoint = as.numeric(sel_ccs_counts_long$timepoint)
  
  sel_ccs_counts_long = sel_ccs_counts_long %>% filter(timepoint >= start_time, timepoint <=stop_time)
  
  kinetic_plot = 
    ggplot(wt_timepoint_pred_df, aes(x = timepoint)) +
    geom_point(data=sel_ccs_counts_long,
               aes(x = timepoint, y = num_cells +  exp(log_abund_detection_thresh), color=expt, alpha=0.5),
               position="jitter", size=0.5) + 
    facet_wrap(~cell_group, scales="free_y", nrow = nrow) + monocle3:::monocle_theme_opts()
  
  if (is.null(batch_col)) {
    
    kinetic_plot = ggplot(wt_timepoint_pred_df, aes(x = timepoint)) +
      geom_point(data=sel_ccs_counts_long,
                 aes(x = timepoint, y = num_cells +  exp(log_abund_detection_thresh)), alpha=0.5,
                 position="jitter", size=0.5) + 
      facet_wrap(~cell_group, scales="free_y", nrow = nrow) + monocle3:::monocle_theme_opts()
    kinetic_plot = kinetic_plot +
      # scale_color_manual(values = my_colors) +
      geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh)), linewidth=1)

  } else {
    kinetic_plot = ggplot(wt_timepoint_pred_df, aes(x = timepoint)) +
      geom_point(data=sel_ccs_counts_long,
                 aes(x = timepoint, y = num_cells + exp(log_abund_detection_thresh), color=!!sym(batch_col)), 
                 alpha=0.5,
                 position="jitter", size=0.5) + 
      facet_wrap(~cell_group, scales="free_y", nrow = nrow) + monocle3:::monocle_theme_opts()
    kinetic_plot = kinetic_plot + 
      scale_color_manual(values = my_colors) +
      geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh), color=!!sym(batch_col)), linewidth=1)
      
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

  if (log_scale)
    kinetic_plot = kinetic_plot + scale_y_log10()
  
  kinetic_plot = kinetic_plot + ylab("cells per sample")

  return(kinetic_plot)
}

#' 
#' @param perturbation_ccm
#' @param cell_groups
#' @param start_time
#' @param stop_time
#' @param interval_col 
#' @param q_val
#' @param control_ccm
#' @param control_state_time
#' @param control_stop_time
#' @param newdata 
#' @export
plot_cell_type_perturb_kinetics = function(perturbation_ccm,
                                           cell_groups=NULL,
                                           start_time=NULL,
                                           stop_time=NULL,
                                           interval_step=1,
                                           log_abund_detection_thresh=-3,
                                           delta_log_abund_loss_thresh=0,
                                           interval_col="timepoint",
                                           q_val=0.01,
                                           log_scale=TRUE,
                                           control_ccm=perturbation_ccm,
                                           control_start_time=start_time,
                                           control_stop_time=control_stop_time,
                                           group_nodes_by = "cell_type",
                                           newdata = tibble()
                                           ){
  
  colData(perturbation_ccm@ccs)[,interval_col] = as.numeric(colData(perturbation_ccm@ccs)[,interval_col])

  if (is.null(start_time))
    start_time = min(colData(perturbation_ccm@ccs)[,interval_col])
  if (is.null(stop_time))
    stop_time = max(colData(perturbation_ccm@ccs)[,interval_col])
  
  if (nrow(newdata) > 0 ){
    newdata = cross_join(tibble(knockout=FALSE) , newdata)
  } else {
    newdata = tibble(knockout=FALSE) 
  }

  wt_timepoint_pred_df = hooke:::estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, 
                                                                   interval_col=interval_col, interval_step=interval_step, newdata = newdata)
  ko_timepoint_pred_df = hooke:::estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time,  
                                                                   interval_col=interval_col, interval_step=interval_step, newdata = newdata)

  #print (wt_timepoint_pred_df)
  #print (ko_timepoint_pred_df)
  timepoints = seq(start_time, stop_time, interval_step)

  # Find the pairs of nodes that are both lost in the perturbation at the same time
  perturb_vs_wt_nodes = tibble(t1=timepoints) %>%
    mutate(comp_abund = purrr::map(.f = platt:::compare_ko_to_wt_at_timepoint,
                                   .x = t1,
                                   perturbation_ccm=perturbation_ccm,
                                   interval_col=interval_col,
                                   wt_pred_df = wt_timepoint_pred_df,
                                   ko_pred_df = ko_timepoint_pred_df)) %>% tidyr::unnest(comp_abund)

  extant_wt_tbl = get_extant_cell_types(perturbation_ccm,
                                        start_time,
                                        stop_time,
                                        log_abund_detection_thresh=log_abund_detection_thresh,
                                        knockout=FALSE,
                                        ...)

  sel_ccs_counts = normalized_counts(perturbation_ccm@ccs, norm_method="size_only", pseudocount=0)
  sel_ccs_counts_long = tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame, var="cell_group") %>%
    pivot_longer(!cell_group, names_to="embryo", values_to="num_cells")

  cell_group_metadata = collect_psg_node_metadata(perturbation_ccm@ccs,
                                                          group_nodes_by=group_nodes_by,
                                                          color_nodes_by=group_nodes_by,
                                                          label_nodes_by=group_nodes_by) %>%
    select(id, cell_group = group_nodes_by)

  sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                  cell_group_metadata,
                                  by=c("cell_group"="id"))

  sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                  colData(perturbation_ccm@ccs) %>% as.data.frame %>% select(sample, !!sym(interval_col), knockout, expt, gene_target),
                                  by=c("embryo"="sample"))

  if (is.null(cell_groups) == FALSE){
    sel_ccs_counts_long = sel_ccs_counts_long %>% filter(cell_group %in% cell_groups)
    perturb_vs_wt_nodes = perturb_vs_wt_nodes %>% filter(cell_group %in% cell_groups)
  }


  perturb_vs_wt_nodes = left_join(perturb_vs_wt_nodes,
                                  extant_wt_tbl %>% select(cell_group, !!sym(interval_col), present_above_thresh), by=c("cell_group" = "cell_group", "t1" = "timepoint"))
  
  if (is.null(cell_groups)){
    cell_groups = unique(as.character(perturb_vs_wt_nodes$cell_group)) %>% sort()
  }
  perturb_vs_wt_nodes$cell_group = factor(as.character(perturb_vs_wt_nodes$cell_group), levels=cell_groups)
  sel_ccs_counts_long$cell_group = factor(as.character(sel_ccs_counts_long$cell_group), levels=cell_groups)

  
  perturb_vs_wt_nodes$t1 = as.numeric(perturb_vs_wt_nodes$t1)
  sel_ccs_counts_long$timepoint = as.numeric(sel_ccs_counts_long$timepoint)
  
  kinetic_plot = ggplot(perturb_vs_wt_nodes, aes(x = !!sym(paste(interval_col, "_x", sep=""))))+ 
    
    geom_point(aes(x = !!sym(paste(interval_col, "_x", sep="")), 
                   y =1.0*(exp(log_abund_y) + exp(log_abund_detection_thresh))), 
               shape=8, 
               data=perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund > 0 & delta_q_value < q_val)) +
    geom_point(aes(x = !!sym(paste(interval_col, "_x", sep="")), y =1.0*(exp(log_abund_y) + exp(log_abund_detection_thresh))), 
               shape=8, data=perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund < 0 & delta_q_value < q_val)) +
    geom_jitter(data=sel_ccs_counts_long,
                aes(x = timepoint, y = num_cells+exp(log_abund_detection_thresh), shape=expt),
                height=0,
                width=2, color = "gray", size=0.5) +
    geom_line(aes(y = exp(log_abund_x) + exp(log_abund_detection_thresh), linetype = "Wild-type")) +
    geom_line(aes(y = exp(log_abund_y) + exp(log_abund_detection_thresh), linetype = "Knockout")) +
    ggh4x::stat_difference(aes(ymin = exp(log_abund_x)+exp(log_abund_detection_thresh), 
                               ymax = exp(log_abund_y) +exp(log_abund_detection_thresh)), alpha=0.5) + 
    scale_fill_manual(values = c("orangered3", "royalblue3", "lightgray")) + 
    facet_wrap(~cell_group, scales="free_y", nrow=2) + monocle3:::monocle_theme_opts() + 
    # ggtitle(paste0(perturbation, " kinetics")) + 
    geom_hline(yintercept=exp(log_abund_detection_thresh), color="lightgrey") + xlab("timepoint") + 
    facet_wrap(~cell_group, scales="free_y", nrow=1)
  
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
  
  if (log_scale)
    kinetic_plot = kinetic_plot + scale_y_log10()

  return(kinetic_plot)
}
