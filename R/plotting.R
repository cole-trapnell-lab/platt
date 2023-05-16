
plot_state_graph_annotations_wrapper = function(ccm, state_graph, edge_weights="support") {
  plot_state_graph_annotations(ccm,
                               state_graph,
                               color_nodes_by = "timepoint",
                               group_nodes_by = "cell_type_sub",
                               edge_weights = edge_weights,
                               hide_unlinked_nodes = FALSE)
}



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
                                           log_scale=TRUE,
                                           ...){

  if (is.null(start_time))
    start_time = min(colData(control_ccm@ccs)[,interval_col])
  if (is.null(stop_time))
    stop_time = max(colData(control_ccm@ccs)[,interval_col])

  if (is.null(batch_col)){
    wt_timepoint_pred_df = hooke:::estimate_abundances_over_interval(control_ccm, start_time, stop_time, knockout=FALSE, interval_col=interval_col, interval_step=interval_step, ...)
  }else{

    batches = tibble(batch = unique(colData(control_ccm@ccs)[,batch_col]))
    batches = batches %>% mutate(tp_preds = purrr::map(.f = function(expt) {
      hooke:::estimate_abundances_over_interval(control_ccm,
                                                start_time,
                                                stop_time,
                                                knockout=FALSE,
                                                interval_col=interval_col,
                                                interval_step=interval_step,
                                                expt, # FIXME: this should be generic, not hardcoded as "expt"
                                                ...)
    }, .x=batch))
    wt_timepoint_pred_df = batches %>% select(tp_preds) %>% tidyr::unnest(tp_preds)
  }

  #print (wt_timepoint_pred_df)
  #print (ko_timepoint_pred_df)
  timepoints = seq(start_time, stop_time, interval_step)

  #print (earliest_loss_tbl)


  peak_wt_abundance = wt_timepoint_pred_df %>% group_by(cell_group) %>% slice_max(log_abund, n=1)

  sel_ccs_counts = normalized_counts(control_ccm@ccs, norm_method="size_only", pseudocount=0)
  #sel_ccs_counts = Matrix::t(Matrix::t(counts(control_ccm@ccs)) / exp(model.offset(control_ccm@model_aux[["model_frame"]])))
  sel_ccs_counts_long = tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame, var="cell_group") %>%
    pivot_longer(!cell_group, names_to="embryo", values_to="num_cells")

  #sel_ccs_counts_long = summary(sel_ccs_counts)
  #sel_ccs_counts_long = data.frame(cell_group = rownames(sel_ccs_counts)[sel_ccs_counts_long$i],
  #                                 embryo = colnames(sel_ccs_counts)[sel_ccs_counts_long$j],
  #                                 num_cells      = sel_ccs_counts_long$x)

  cell_group_metadata = hooke:::collect_psg_node_metadata(control_ccm@ccs,
                                                          group_nodes_by="cell_type_broad",
                                                          color_nodes_by="cell_type_broad",
                                                          label_nodes_by="cell_type_broad") %>%
    select(id, cell_type_broad = group_nodes_by)

  sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                  cell_group_metadata,
                                  by=c("cell_group"="id"))

  sel_ccs_counts_long = left_join(sel_ccs_counts_long,
                                  colData(control_ccm@ccs) %>% as.data.frame %>%
                                    select(sample, !!sym(interval_col), !!sym(batch_col)),
                                  by=c("embryo"="sample"))

  if (is.null(cell_groups) == FALSE){
    wt_timepoint_pred_df = wt_timepoint_pred_df %>% filter(cell_group %in% cell_groups)
    sel_ccs_counts_long = sel_ccs_counts_long %>% filter(cell_group %in% cell_groups)
  }

  wt_timepoint_pred_df$cell_group = factor(as.character(wt_timepoint_pred_df$cell_group), levels=cell_groups)
  sel_ccs_counts_long$cell_group = factor(as.character(sel_ccs_counts_long$cell_group), levels=cell_groups)

  kinetic_plot = ggplot(wt_timepoint_pred_df, aes(x = !!sym(interval_col))) +
    geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh), color=!!sym(batch_col))) +
    facet_wrap(~cell_group, scales="free_y", ncol=1)

  kinetic_plot = kinetic_plot +
    geom_point(data=sel_ccs_counts_long,
               aes(x = !!sym(interval_col), y = num_cells +  exp(log_abund_detection_thresh), color=!!sym(batch_col)),
               position="jitter") +
    # facet_wrap(~cell_group, scales="free_y", ncol=1)
    facet_wrap(~cell_group, scales="free_y", nrow=1)
  if (log_scale)
    kinetic_plot = kinetic_plot + scale_y_log10()


  return(kinetic_plot)
}


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
                                           ...){

  if (is.null(start_time))
    start_time = min(colData(perturbation_ccm@ccs)[,interval_col])
  if (is.null(stop_time))
    stop_time = max(colData(perturbation_ccm@ccs)[,interval_col])

  wt_timepoint_pred_df = hooke:::estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=FALSE, interval_col=interval_col, interval_step=interval_step, ...)
  ko_timepoint_pred_df = hooke:::estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=TRUE, interval_col=interval_col, interval_step=interval_step, ...)

  #print (wt_timepoint_pred_df)
  #print (ko_timepoint_pred_df)
  timepoints = seq(start_time, stop_time, interval_step)

  # Find the pairs of nodes that are both lost in the perturbation at the same time
  perturb_vs_wt_nodes = tibble(t1=timepoints) %>%
    mutate(comp_abund = purrr::map(.f = hooke:::compare_ko_to_wt_at_timepoint,
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

  cell_group_metadata = hooke:::collect_psg_node_metadata(perturbation_ccm@ccs,
                                                          group_nodes_by="cell_type_broad",
                                                          color_nodes_by="cell_type_broad",
                                                          label_nodes_by="cell_type_broad") %>%
    select(id, cell_type_broad = group_nodes_by)

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

  perturb_vs_wt_nodes$cell_group = factor(as.character(perturb_vs_wt_nodes$cell_group), levels=cell_groups)
  sel_ccs_counts_long$cell_group = factor(as.character(sel_ccs_counts_long$cell_group), levels=cell_groups)


  kinetic_plot = ggplot(perturb_vs_wt_nodes, aes(x = !!sym(paste(interval_col, "_x", sep="")))) +
    geom_line(aes(y = exp(log_abund_x) +exp(log_abund_detection_thresh), linetype = "Wild-type")) +
    geom_line(aes(y = exp(log_abund_y) +exp(log_abund_detection_thresh), linetype = "Knockout")) +
    ggh4x::stat_difference(aes(ymin = exp(log_abund_x)+exp(log_abund_detection_thresh), ymax = exp(log_abund_y) +exp(log_abund_detection_thresh)), alpha=0.3) +
    geom_point(aes(x = !!sym(paste(interval_col, "_x", sep="")), y =1.0*(exp(log_abund_y) + exp(log_abund_detection_thresh))), shape=8, data=perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund > 0 & delta_q_value < q_val)) +
    geom_point(aes(x = !!sym(paste(interval_col, "_x", sep="")), y =1.0*(exp(log_abund_y) + exp(log_abund_detection_thresh))), shape=8, data=perturb_vs_wt_nodes %>% filter(present_above_thresh & delta_log_abund < 0 & delta_q_value < q_val)) +
    # facet_wrap(~cell_group, scales="free_y", ncol=1)
    facet_wrap(~cell_group, scales="free_y", nrow=1)

  kinetic_plot = kinetic_plot +
    geom_jitter(data=sel_ccs_counts_long,
                aes(x = !!sym(interval_col), y = num_cells+exp(log_abund_detection_thresh), color=gene_target, shape=expt),
                height=0,
                width=2)
  kinetic_plot = kinetic_plot + geom_hline(yintercept=exp(log_abund_detection_thresh), color="lightgrey")
  if (log_scale)
    kinetic_plot = kinetic_plot + scale_y_log10()

  return(kinetic_plot)
}

geom_richnodelabel <- function (mapping = NULL, data = NULL, position = "identity",
                                ..., parse = FALSE, nudge_x = 0, nudge_y = 0, label.padding = unit(0.25,
                                                                                                   "lines"), label.r = unit(0.15, "lines"), label.size = 0.25,
                                na.rm = FALSE, show.legend = NA, inherit.aes = TRUE)
{
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`",
           call. = FALSE)
    }
    position <- ggplot2::position_nudge(nudge_x, nudge_y)
  }
  ggplot2::layer(data = data, mapping = mapping, stat = ggnetwork:::StatNodes,
                 geom = ggtext:::GeomRichtext, position = position, show.legend = show.legend,
                 inherit.aes = inherit.aes, params = list(parse = parse,
                                                          label.padding = label.padding, label.r = label.r,
                                                          label.size = label.size, na.rm = na.rm, ...))
}

plot_state_graph_perturb_effects <- function(ccs,
                                             state_graph,
                                             perturbation_table,
                                             patterns = c("direct", "indirect", "inferred", "predicted"),
                                             num_top_perturbs=3,
                                             num_top_genes=3,
                                             color_nodes_by=NULL,
                                             label_nodes_by=NULL,
                                             group_nodes_by=NULL,
                                             label_edges_by=NULL,
                                             edge_weights=NULL,
                                             arrow.gap=0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             node_size = 2,
                                             num_layers=10,
                                             min_edge_size=0.1,
                                             max_edge_size=2,
                                             fract_expr = 0.0,
                                             mean_expr = 0.0,
                                             unlabeled_groups = c("Unknown"),
                                             label_groups=TRUE,
                                             hide_unlinked_nodes=TRUE,
                                             group_label_font_size=6,
                                             edge_label_font_size=2,
                                             label_conn_linetype="dotted",
                                             legend_position = "none",
                                             con_colour = "darkgrey",
                                             group_outline=FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = hooke:::collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
  }

  node_support_df = igraph::as_data_frame(state_graph, what="vertices") %>%
    tibble::as_tibble() %>% rename(id=name)
  #node_metadata = node_metadata %>% tidyr::unnest(direct_perturb)

  perturbation_table = perturbation_table %>% group_by(id) %>%
    mutate(perturb_display_name = case_when(perturb_effect == "direct" ~ glue::glue("<i style='color:#2B9EB3'>{perturb_name}</i>"),
                                            perturb_effect == "indirect" ~ glue::glue("<i style='color:#FCAB10'>{perturb_name}</i>"),
                                            perturb_effect == "predicted" ~ glue::glue("<i style='color:#F8333C'>{perturb_name}</i>"),
                                            TRUE ~ glue::glue("<i style='color:#44AF69'>{perturb_name}</i>"))) %>%
    summarize(perturb_effect_label = ifelse(n() > num_top_genes,
                                            paste0(c(perturb_display_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "<br>"),
                                            paste0(perturb_display_name, collapse = "<br>")))
  # summarize(top_genes=ifelse(n() > num_top_perturbs,
  #                            paste0(c(gene_short_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "\n"),
  #                            paste0(gene_short_name, collapse = "\n"))) %>%
  #

  node_metadata = node_metadata %>% left_join(perturbation_table %>% select(id, perturb_effect_label), by=c("id"="id"))
  if (is.null(label_nodes_by)){
    label_nodes_by = "perturb_effect_label"
    node_metadata = node_metadata %>% mutate(label_nodes_by=perturb_effect_label)
  }else{
    node_metadata = node_metadata %>% mutate(label_nodes_by=ifelse(is.na(perturb_effect_label), label_nodes_by, paste(label_nodes_by, perturb_effect_label, sep="<br>")))
  }

  if (is.null(label_edges_by)){
    edge_info = edges %>% select(from, to)
  }else{
    if (is(state_graph, "igraph")){
      edge_info = state_graph %>% igraph::as_data_frame() %>% select(from, to, label=!!sym(label_edges_by))
    }else{
      edge_info = state_graph %>% select(from, to, label=!!sym(label_edges_by))
    }

    edges = edges %>% left_join(edge_info)
    #print(edges)
  }

  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE){
    G_df = igraph::as_data_frame(G)
    edge_names =  stringr::str_c(G_df$from, G_df$to, sep="~")
    edge_labels = igraph::E(G)$label
    names(edge_labels) = edge_names
    #print(edge_labels)
    #edge_labels = NULL
  }else{
    edge_labels=NULL
  }

  layout_info = hooke:::layout_state_graph(G, node_metadata, edge_labels, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
                                     unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
                                     edge_thickness = replace_na(edge_thickness, min_edge_size))
  }else{
    bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge = FALSE
  }

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){
    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      label.fontsize=group_label_font_size,
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      data=g)
    } else {
      if (label_groups){
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            label = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        expand = unit(2, "mm"),
                                        label.buffer=unit(1, "mm"),
                                        radius = unit(1.5, "mm"),
                                        label.margin = margin(1, 1, 1, 1, "mm"),
                                        label.fontsize=group_label_font_size,
                                        label.fontface="plain",
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }else{
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        label.fontsize=group_label_font_size,
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }

    }

  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE){
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          geom_richnodelabel(data = g,
                             aes(x, y,
                                 fill = !!sym(color_nodes_by),
                                 label = label_nodes_by),
                             size = node_size) +
          labs(fill = color_nodes_by)
        p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + geom_richnodelabel(data = g,
                                   aes(x, y,
                                       fill = color_nodes_by,
                                       label = label_nodes_by),
                                   size = node_size) +
          labs(fill = color_nodes_by)

      }
    }else{
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(data = g,
                                aes(x, y,
                                    color = !!sym(color_nodes_by)),
                                size = node_size) +
          labs(color = color_nodes_by)
        p = p + scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + ggnetwork::geom_nodes(data = g,
                                      aes(x, y,
                                          color = color_nodes_by),
                                      size = node_size) +
          labs(color = color_nodes_by)
      }
    }

  } else {
    if (is.null(label_nodes_by) == FALSE){
      p = p + geom_richnodelabel(data = g,
                                 aes(x, y, xend = xend, yend = yend,
                                     label = label_nodes_by),
                                 size = node_size)
    }else{
      p = p + ggnetwork::geom_nodes(data = g,
                                    aes(x, y, xend = xend, yend = yend),
                                    size = node_size)
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p = p + geom_text(data = label_df,
                      aes(x,y, label = label),
                      size=edge_label_font_size)
  }

  p = p + annotate(geom='richtext', x=0.15*max(g$xend), y=0.01*max(g$yend),
                   size=2,
                   fill = NA, label.color = NA,
                   label="<i style='color:#2B9EB3'>direct</i> <i style='color:#FCAB10'>indirect</i> <i style='color:#F8333C'>predicted</i> <i style='color:#44AF69'>other</i> ")
  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)
  return(p)
}


plot_state_graph_marker_genes <- function(ccs,
                                             state_graph,
                                             marker_table,
                                             patterns_for_marker_genes = c("Specifically activated",
                                                                           "Selectively activated",
                                                                           "Selectively upregulated",
                                                                           "Specifically upregulated",
                                                                           "Specifically maintained",
                                                                           "Selectively maintained",
                                                                           "Activated",
                                                                           "Upregulated",
                                                                           "Maintained"),
                                             num_top_genes=3,
                                             color_nodes_by=NULL,
                                             label_nodes_by=NULL,
                                             group_nodes_by=NULL,
                                             label_edges_by=NULL,
                                             edge_weights=NULL,
                                             arrow.gap=0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             node_size = 2,
                                             num_layers=10,
                                             min_edge_size=0.1,
                                             max_edge_size=2,
                                             fract_expr = 0.0,
                                             mean_expr = 0.0,
                                             unlabeled_groups = c("Unknown"),
                                             label_groups=TRUE,
                                             hide_unlinked_nodes=TRUE,
                                             group_label_font_size=6,
                                             edge_label_font_size=2,
                                             label_conn_linetype="dotted",
                                             legend_position = "none",
                                             con_colour = "darkgrey",
                                             group_outline=FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = hooke:::collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
  }

  node_support_df = igraph::as_data_frame(state_graph, what="vertices") %>%
    tibble::as_tibble() %>% rename(id=name)
  #node_metadata = node_metadata %>% tidyr::unnest(direct_perturb)

  red_patterns =  c("Specifically activated",
                    "Specifically upregulated",
                    "Specifically maintained")
  orange_patterns = c("Selectively activated",
                      "Selectively upregulated",
                      "Selectively maintained")
  blue_patterns = c("Activated",
                     "Upregulated")
  green_patterns = c("Maintained")

  patterns_for_marker_genes = c("Specifically activated",
                                "Selectively activated",
                                "Selectively upregulated",
                                "Specifically upregulated",
                                "Specifically maintained",
                                "Selectively maintained",
                                "Activated",
                                "Upregulated",
                                "Maintained")


  marker_table = marker_table %>% group_by(id) %>%
    mutate(marker_display_name = case_when(marker_type %in% blue_patterns ~ glue::glue("<i style='color:#2B9EB3'>{marker_name}</i>"),
                                           marker_type %in% orange_patterns ~ glue::glue("<i style='color:#FCAB10'>{marker_name}</i>"),
                                           marker_type %in% red_patterns ~ glue::glue("<i style='color:#F8333C'>{marker_name}</i>"),
                                           marker_type %in% green_patterns ~ glue::glue("<i style='color:#44AF69'>{marker_name}</i>"),
                                            TRUE ~ glue::glue("<i style='color:#D3D3D3'>{marker_name}</i>"))) %>%
    summarize(marker_label = ifelse(n() > num_top_genes,
                                            paste0(c(marker_display_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "<br>"),
                                            paste0(marker_display_name, collapse = "<br>")))
  # summarize(top_genes=ifelse(n() > num_top_perturbs,
  #                            paste0(c(gene_short_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "\n"),
  #                            paste0(gene_short_name, collapse = "\n"))) %>%
  #

  node_metadata = node_metadata %>% left_join(marker_table %>% select(id, marker_label), by=c("id"="id"))
  if (is.null(label_nodes_by)){
    label_nodes_by = "marker_label"
    node_metadata = node_metadata %>% mutate(label_nodes_by=marker_label)
  }else{
    node_metadata = node_metadata %>% mutate(label_nodes_by=ifelse(is.na(marker_label), label_nodes_by, paste(label_nodes_by, marker_label, sep="<br>")))
  }

  if (is.null(label_edges_by)){
    edge_info = edges %>% select(from, to)
  }else{
    if (is(state_graph, "igraph")){
      edge_info = state_graph %>% igraph::as_data_frame() %>% select(from, to, label=!!sym(label_edges_by))
    }else{
      edge_info = state_graph %>% select(from, to, label=!!sym(label_edges_by))
    }

    edges = edges %>% left_join(edge_info)
    #print(edges)
  }

  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE){
    G_df = igraph::as_data_frame(G)
    edge_names =  stringr::str_c(G_df$from, G_df$to, sep="~")
    edge_labels = igraph::E(G)$label
    names(edge_labels) = edge_names
    #print(edge_labels)
    #edge_labels = NULL
  }else{
    edge_labels=NULL
  }

  layout_info = hooke:::layout_state_graph(G, node_metadata, edge_labels, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
                                     unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
                                     edge_thickness = replace_na(edge_thickness, min_edge_size))
  }else{
    bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge = FALSE
  }

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){
    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      label.fontsize=group_label_font_size,
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      data=g)
    } else {
      if (label_groups){
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            label = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        expand = unit(2, "mm"),
                                        label.buffer=unit(1, "mm"),
                                        radius = unit(1.5, "mm"),
                                        label.margin = margin(1, 1, 1, 1, "mm"),
                                        label.fontsize=group_label_font_size,
                                        label.fontface="plain",
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }else{
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        label.fontsize=group_label_font_size,
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }

    }

  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE){
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          geom_richnodelabel(data = g,
                             aes(x, y,
                                 fill = !!sym(color_nodes_by),
                                 label = label_nodes_by),
                             size = node_size) +
          labs(fill = color_nodes_by)
        p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + geom_richnodelabel(data = g,
                                   aes(x, y,
                                       fill = color_nodes_by,
                                       label = label_nodes_by),
                                   size = node_size) +
          labs(fill = color_nodes_by)

      }
    }else{
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(data = g,
                                aes(x, y,
                                    color = !!sym(color_nodes_by)),
                                size = node_size) +
          labs(color = color_nodes_by)
        p = p + scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + ggnetwork::geom_nodes(data = g,
                                      aes(x, y,
                                          color = color_nodes_by),
                                      size = node_size) +
          labs(color = color_nodes_by)
      }
    }

  } else {
    if (is.null(label_nodes_by) == FALSE){
      p = p + geom_richnodelabel(data = g,
                                 aes(x, y, xend = xend, yend = yend,
                                     label = label_nodes_by),
                                 size = node_size)
    }else{
      p = p + ggnetwork::geom_nodes(data = g,
                                    aes(x, y, xend = xend, yend = yend),
                                    size = node_size)
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p = p + geom_text(data = label_df,
                      aes(x,y, label = label),
                      size=edge_label_font_size)
  }

  p = p + annotate(geom='richtext', x=0.15*max(g$xend), y=0.01*max(g$yend),
                   size=2,
                   fill = NA, label.color = NA,
                   label="<i style='color:#2B9EB3'>upregulated</i> <i style='color:#FCAB10'>selective</i> <i style='color:#F8333C'>specific</i> <i style='color:#44AF69'>maintained</i> ")
  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)
  return(p)
}


plot_state_graph_key_genes <- function(ccs,
                                       state_graph,
                                       gene_pattern_table,
                                       num_top_genes=3,
                                       patterns=c("Activated",
                                                  "Upregulated",
                                                  "Specifically activated",
                                                  "Selectively activated",
                                                  "Specifically activated",
                                                  "Selectively upregulated",
                                                  "Specifically upregulated",
                                                  "Specifically maintained",
                                                  "Selectively maintained"
                                       ),
                                       color_nodes_by=NULL,
                                       label_nodes_by=NULL,
                                       group_nodes_by=NULL,
                                       label_edges_by=NULL,
                                       edge_weights=NULL,
                                       arrow.gap=0.03,
                                       arrow_unit = 2,
                                       bar_unit = .075,
                                       node_size = 2,
                                       num_layers=10,
                                       min_edge_size=0.1,
                                       max_edge_size=2,
                                       fract_expr = 0.0,
                                       mean_expr = 0.0,
                                       unlabeled_groups = c("Unknown"),
                                       label_groups=TRUE,
                                       hide_unlinked_nodes=TRUE,
                                       group_label_font_size=6,
                                       edge_label_font_size=2,
                                       label_conn_linetype="dotted",
                                       legend_position = "none",
                                       con_colour = "darkgrey",
                                       group_outline=FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = hooke:::collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
  }

  selected_genes = gene_pattern_table %>% filter(interpretation %in% patterns) %>%
    group_by(cell_state) %>% arrange(desc(pattern_match_score)) %>%
    summarize(top_genes=ifelse(n() > num_top_genes,
                               paste0(c(gene_short_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "\n"),
                               paste0(gene_short_name, collapse = "\n")))

  node_metadata = node_metadata %>% left_join(selected_genes %>% select(cell_state, top_genes), by=c("id"="cell_state"))
  if (is.null(label_nodes_by)){
    label_nodes_by = "top_genes"
    node_metadata = node_metadata %>% mutate(label_nodes_by=top_genes)
  }else{
    node_metadata = node_metadata %>% mutate(label_nodes_by=paste(label_nodes_by, top_genes, sep="\n"))
  }

  if (is.null(label_edges_by)){
    edge_info = edges %>% select(from, to)
  }else{
    if (is(state_graph, "igraph")){
      edge_info = state_graph %>% igraph::as_data_frame() %>% select(from, to, label=!!sym(label_edges_by))
    }else{
      edge_info = state_graph %>% select(from, to, label=!!sym(label_edges_by))
    }

    edges = edges %>% left_join(edge_info)
    #print(edges)
  }

  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE){
    G_df = igraph::as_data_frame(G)
    edge_names =  stringr::str_c(G_df$from, G_df$to, sep="~")
    edge_labels = igraph::E(G)$label
    names(edge_labels) = edge_names
    #print(edge_labels)
    #edge_labels = NULL
  }else{
    edge_labels=NULL
  }

  layout_info = hooke:::layout_state_graph(G, node_metadata, edge_labels, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size)
  }else{
    bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
  }

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){
    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      label.fontsize=group_label_font_size,
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      data=g)
    } else {
      if (label_groups){
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            label = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        expand = unit(2, "mm"),
                                        label.buffer=unit(1, "mm"),
                                        radius = unit(1.5, "mm"),
                                        label.margin = margin(1, 1, 1, 1, "mm"),
                                        label.fontsize=group_label_font_size,
                                        label.fontface="plain",
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }else{
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        label.fontsize=group_label_font_size,
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }

    }

  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE){
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodelabel(data = g,
                                    aes(x, y,
                                        fill = !!sym(color_nodes_by),
                                        label = label_nodes_by),
                                    size = node_size) +
          labs(fill = color_nodes_by)
        p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + ggnetwork::geom_nodelabel(data = g,
                                          aes(x, y,
                                              fill = color_nodes_by,
                                              label = label_nodes_by),
                                          size = node_size) +
          labs(fill = color_nodes_by)

      }
    }else{
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(data = g,
                                aes(x, y,
                                    color = !!sym(color_nodes_by)),
                                size = node_size) +
          labs(color = color_nodes_by)
        p = p + scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + ggnetwork::geom_nodes(data = g,
                                      aes(x, y,
                                          color = color_nodes_by),
                                      size = node_size) +
          labs(color = color_nodes_by)
      }
    }

  } else {
    if (is.null(label_nodes_by) == FALSE){
      p = p + ggnetwork::geom_nodelabel(data = g,
                                        aes(x, y, xend = xend, yend = yend,
                                            label = label_nodes_by),
                                        size = node_size)
    }else{
      p = p + ggnetwork::geom_nodes(data = g,
                                    aes(x, y, xend = xend, yend = yend),
                                    size = node_size)
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p = p + geom_text(data = label_df,
                      aes(x,y, label = label),
                      size=edge_label_font_size)
  }

  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)
  return(p)
}
