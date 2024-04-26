
#' Plot the state graph with annotations.
#'
#' @param ccm cell_count_model A cell_count_model.
#' @param state_graph
#' @param color_nodes_by
#' @param label_nodes_by
#' @param group_nodes_by
#' @param label_edges_by
#' @param edge_weights
#' @param arrow.gap
#' @param arrow_unit
#' @param bar_unit
#' @param node_size
#' @param min_edge_size
#' @param max_edge_size
#' @param unlabeled_groups
#' @param label_groups
#' @param hide_unlinked_nodes
#' @param group_label_font_size
#' @param edge_label_font_size
#' @param label_conn_linetype
#' @param legend_position
#' @param con_colour
#' @param group_outline
#' @return A ggplot2 plot object.
#' @export
plot_state_graph_annotations <- function(ccs,
                                         state_graph,
                                         color_nodes_by=NULL,
                                         label_nodes_by=NULL,
                                         group_nodes_by=NULL,
                                         label_edges_by=NULL,
                                         edge_weights=NULL,
                                         arrow.gap=0.03,
                                         arrow_unit = 2,
                                         bar_unit = .075,
                                         node_size = 2,
                                         min_edge_size=0.1,
                                         max_edge_size=2,
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

  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
    edges$weight = 1
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
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


  if (length(setdiff(edges$to, node_metadata$id)) > 0 ||
      length(setdiff(edges$from, node_metadata$id)) > 0){
    message("Warning: graph edges refers to nodes not in cell_count_set ccs. Dropping edges...")
  }

  edges = edges %>% filter(from %in% node_metadata$id & to %in% node_metadata$id)

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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
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
  
  # add assembly group
  bezier_df$assembly_group = vapply(strsplit(bezier_df$from, "-", fixed=T), "[", "", 1)
  g$assembly_group = vapply(strsplit(g$name, "-", fixed=T), "[", "", 1)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){
    if (group_outline) {
      # if (label_groups) {
      #   label_groups_by = group_nodes_by
      # } else {
      #   label_groups_by = NULL
      # }
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          # label = label_groups_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      expand = unit(2, "mm"),
                                      label.buffer=unit(1, "mm"),
                                      radius = unit(1.5, "mm"),
                                      label.margin = margin(1, 1, 1, 1, "mm"),
                                      label.fontsize=group_label_font_size,
                                      label.fontface="plain",
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
                                        expand = unit(2, "mm"),
                                        label.buffer=unit(1, "mm"),
                                        radius = unit(1.5, "mm"),
                                        label.margin = margin(1, 1, 1, 1, "mm"),
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
    label_df = layout_info$bezier_df %>%
      group_by(edge_name) %>%
      summarize(x = mean(x), y=mean(y))
    label_df$label = edge_labels[label_df$edge_name]
    #label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    #p = p + geom_text(data = label_df,
    #                  aes(x,y, label = label),
    #                  size=edge_label_font_size)
    p = p + ggrepel::geom_text_repel(data = label_df,
                                mapping = aes(x, y, label=label),
                                size=edge_label_font_size)
  }

  p = p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position)
  return(p)
}

#' This function plots when each cell state is lost following a perturbation
#'
#' FIXME: this is really a prototype and needs lots of cleanup and refactoring.
#' It could plot the earliest time a loss is detected, the latest, etc. But
#' right now it only plots the loss at the time of peak abundance in the control
#' see estimate_loss_timing() for more details.
#'
#' @export
plot_state_graph_losses <- function(perturbation_ccm,
                                    state_graph,
                                    start_time,
                                    stop_time,
                                    interval_step,
                                    interval_col,
                                    log_abund_detection_thresh,
                                    q_val,
                                    loss_time = c("largest_loss", "largest_loss_time", "earliest_time", "latest_time", "peak_loss_time", "delta_log_abund_at_peak"),
                                    label_nodes_by=NULL,
                                    group_nodes_by=NULL,
                                    label_edges_by=NULL,
                                    edge_weights=NULL,
                                    arrow.gap=0.03,
                                    arrow_unit = 2,
                                    bar_unit = .075,
                                    node_size = 2,
                                    min_edge_size=0.1,
                                    max_edge_size=2,
                                    unlabeled_groups = c("Unknown"),
                                    label_subset = NULL,
                                    label_groups=TRUE,
                                    hide_unlinked_nodes=TRUE,
                                    group_label_font_size=6,
                                    edge_label_font_size=2,
                                    label_conn_linetype="dotted",
                                    legend_position = "none",
                                    con_colour = "darkgrey",
                                    group_outline=FALSE,
                                    control_ccm=perturbation_ccm,
                                    control_start_time=NULL,
                                    control_stop_time=NULL,
                                    newdata = tibble())
{
  loss_time = match.arg(loss_time)

  if (is.null(control_start_time))
    control_start_time = start_time
  if (is.null(control_stop_time))
    control_stop_time = stop_time

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = collect_psg_node_metadata(perturbation_ccm@ccs, color_nodes_by=NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
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

  earliest_loss_tbl = estimate_loss_timing(perturbation_ccm,
                                           start_time=start_time,
                                           stop_time=stop_time,
                                           interval_step = interval_step,
                                           interval_col=interval_col,
                                           log_abund_detection_thresh=log_abund_detection_thresh,
                                           q_val = q_val,
                                           control_ccm=control_ccm,
                                           control_start_time=control_start_time,
                                           control_stop_time=control_stop_time,
                                           newdata = newdata)

  #earliest_loss_tbl = earliest_loss_tbl %>% mutate(fill_alpha = ifelse(peak_time_in_ctrl_within_perturb_time_range, 1.0, 0.3))
  #print (earliest_loss_tbl)
  node_metadata = node_metadata %>% left_join(earliest_loss_tbl, by=c("id" ="cell_group"))

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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
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
  # add assembly group
  bezier_df$assembly_group = vapply(strsplit(bezier_df$from, "-", fixed=T), "[", "", 1)
  g$assembly_group = vapply(strsplit(g$name, "-", fixed=T), "[", "", 1)
  
  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness), data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  if (is.null(label_subset)) {
    label_subset = unique(node_metadata$group_nodes_by)
  }

  label_subset = label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE){


    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      expand = unit(2, "mm"),
                                      label.buffer=unit(1, "mm"),
                                      radius = unit(1.5, "mm"),
                                      label.margin = margin(1, 1, 1, 1, "mm"),
                                      label.fontsize=group_label_font_size,
                                      label.fontface="plain",
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      data=g)
    } else {
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
    }

  }


  # if numerical

  min_delta_log_abund = -3
  max_delta_log_abund = 3
  g = g %>% mutate(fill_val = !!sym(loss_time))

  if (loss_time %in% c("largest_loss", "delta_log_abund_at_peak")){
    g = g %>% mutate(fill_val = ifelse(fill_val < min_delta_log_abund, min_delta_log_abund, fill_val),
                     fill_val = ifelse(fill_val > max_delta_log_abund, max_delta_log_abund, fill_val),
                     fill_val = ifelse(peak_time_in_ctrl_within_perturb_time_range, fill_val, NA),
                     fill_border_color = ifelse(peak_time_in_ctrl_within_perturb_time_range == FALSE,
                                                "out-of-range",
                                                ifelse(is_lost_at_peak,
                                                  "significant",
                                                  "notsignificant")
                                                ))
  }

  p = p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodelabel(data = g,
                              aes(x, y,
                                  fill = fill_val,
                                  color = fill_border_color,
                                  label = label_nodes_by),
                              size = node_size)
  p = p + scale_color_manual(values=c("out-of-range"="black", "significant"="black", "notsignificant"="lightgrey"))

  if (loss_time %in% c("largest_loss", "delta_log_abund_at_peak")){
    p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3",
                                 limits=c(min_delta_log_abund - abs(min_delta_log_abund)*0.05,
                                          max_delta_log_abund + abs(max_delta_log_abund)*0.05))
  }else{
    p = p + scale_fill_stepsn(n.breaks=5, colours = terrain.colors(5)) #+ scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
  }

  p = p + guides(color = "none")

  #p = p + guides(fill = guide_legend(position=legend_position))

  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$bezier_df %>%
      group_by(edge_name) %>%
      summarize(x = mean(x), y=mean(y))
    label_df$label = edge_labels[label_df$edge_name]
    #label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    #p = p + geom_text(data = label_df,
    #                  aes(x,y, label = label),
    #                  size=edge_label_font_size)
    p = p + ggrepel::geom_text_repel(data = label_df,
                                     mapping = aes(x, y, label=label),
                                     size=edge_label_font_size)
  }

  p = p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() #+
    #theme(legend.position=legend_position)
  return(p)
}


#' @noRd
num_extract <- function(string, as_char = TRUE){

  if (as_char)
    stringr::str_trim(
      format(stringr::str_extract(string, "\\-*\\d+\\.*\\d*"),
             scientific = FALSE,
             trim = TRUE)
    )

  else
    as.numeric(
      stringr::str_trim(
        format(stringr::str_extract(string, "\\-*\\d+\\.*\\d*"),
               scientific = FALSE,
               trim = TRUE)
      )
    )

}

#' @noRd
calc_sig_ind <- function(p_value, html = TRUE) {

  #p_value <- suppressWarnings(
  #  num_extract(p_value, as_char = FALSE)
  #)

  if (html) {
    dplyr::case_when(
      p_value <= 0 ~ "",
      p_value <= 0.001 ~ "\\***",
      p_value <= 0.01 ~ "\\**",
      p_value <= 0.05 ~ "\\*",
      p_value <= 0.1 ~ ".",
      p_value <= 1 ~ "",
      TRUE ~ ""
    )
  } else {
    dplyr::case_when(
      p_value <= 0 ~ "",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      p_value <= 0.1 ~ ".",
      p_value <= 1 ~ "",
      TRUE ~ ""
    )
  }

}

#' Plot the expression of a single gene on the state transition graph
#' @export
plot_state_graph_abundance_changes <- function(ccs,
                                               state_graph,
                                               comp_abund_table,
                                               contrast = "contrast",
                                               label_nodes_by=NULL,
                                               group_nodes_by=NULL,
                                               label_edges_by=NULL,
                                               edge_weights=NULL,
                                               size_nodes_by = NULL,
                                               num_degs = NULL, 
                                               fc_limits=c(-3,3),
                                               arrow.gap=0.03,
                                               arrow_unit = 2,
                                               bar_unit = .075,
                                               node_size = 6,
                                               min_edge_size=0.1,
                                               max_edge_size=2,
                                               unlabeled_groups = c("Unknown"),
                                               label_subset = NULL,
                                               label_groups=TRUE,
                                               hide_unlinked_nodes=TRUE,
                                               group_label_font_size=6,
                                               edge_label_font_size=2,
                                               label_conn_linetype="dotted",
                                               legend_position = "none",
                                               con_colour = "darkgrey",
                                               group_outline=FALSE,
                                               flip_x = FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by=NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)

  comp_abund_table[["contrast"]] = comp_abund_table[[contrast]]
  comp_abund_table$contrast = as.factor(comp_abund_table$contrast)

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

  abund_fc_df = comp_abund_table

  if (is.null(fc_limits)) {
    fc_limits = range(abund_fc_df$delta_log_abund)
  } else {
    min = fc_limits[1]
    max = fc_limits[2]
    abund_fc_df = abund_fc_df %>%
      mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
      mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  }

  abund_fc_df = abund_fc_df %>%
    mutate(
      delta_q_value = pmax(0.0001, delta_q_value),
      q_value_sig_code = calc_sig_ind(delta_q_value, html=FALSE))

  g = left_join(g, abund_fc_df, by=c("name"="cell_group"), relationship = "many-to-many")

  #g$gene_short_name = factor(g$gene_short_name, levels=genes)
  color_nodes_by = "delta_log_abund"
  # group_outline = TRUE

  if (flip_x) {
    g$x = -1*g$x
    bezier_df$x = -1*bezier_df$x
  }

  p <- ggplot(aes(x,y), data=g)
  p = p +
    ggplot2::geom_path(aes(x, y, group=edge_name), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  if (is.null(label_subset)) {
    label_subset = unique(node_metadata$group_nodes_by)
  }

  label_subset = label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE){

    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          fill = contrast,
                                          color = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      expand = unit(2, "mm"),
                                      label.buffer=unit(1, "mm"),
                                      radius = unit(1.5, "mm"),
                                      label.margin = margin(1, 1, 1, 1, "mm"),
                                      label.fontsize=group_label_font_size,
                                      label.fontface="plain",
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      show.legend = F)
        
    } else {
      values = rep("white",length(label_subset))
      names(values) = label_subset
      if (label_groups) {
        p = p +
          ggforce::geom_mark_rect(aes(x, y,
                                      fill = contrast,
                                      color = group_nodes_by,
                                      label = group_nodes_by,
                                      filter = group_nodes_by %in% label_subset),
                                  size=0.5,
                                  expand = unit(2, "mm"),
                                  label.buffer=unit(1, "mm"),
                                  radius = unit(1.5, "mm"),
                                  label.margin = margin(1, 1, 1, 1, "mm"),
                                  label.fontsize=group_label_font_size,
                                  label.fontface="plain",
                                  con.linetype=label_conn_linetype,
                                  con.colour=con_colour,
                                  show.legend = F) +
          scale_color_manual(values = c(values))
        
      } else {
        p = p +
          ggforce::geom_mark_rect(aes(x, y,
                                      fill = contrast,
                                      color = group_nodes_by,
                                      filter = group_nodes_by %in% label_subset),
                                  size=0.5,
                                  expand = unit(2, "mm"),
                                  label.buffer=unit(1, "mm"),
                                  radius = unit(1.5, "mm"),
                                  label.margin = margin(1, 1, 1, 1, "mm"),
                                  label.fontsize=group_label_font_size,
                                  label.fontface="plain",
                                  con.linetype=label_conn_linetype,
                                  con.colour=con_colour,
                                  show.legend = F) +
          scale_color_manual(values = c(values))
      }

      

    }

    p = p + scale_fill_manual(values=rep("white", length(unique(g$contrast))))
    # p = p + scale_color_manual(values = rainbow_colors)

  }
  p = p + guides(fill = "none")
  
  if (!is.null(size_nodes_by)){
    
    
    g = left_join(g, num_degs, by=c("name"="cell_group", "perturb_name"))
    g[["size_nodes_by"]] = as.numeric(g[[size_nodes_by]] )
    p = p + ggnewscale::new_scale_color() +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                size = size_nodes_by*1.2),
                            color=I("black")) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                size = size_nodes_by,
                                color=delta_log_abund)) +
      ggnetwork::geom_nodetext(data = g,
                               aes(x, y,
                                   label = q_value_sig_code),
                               color=I("black"))
    
  } else {
    p = p + ggnewscale::new_scale_color() +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                size = -log10(delta_q_value)*1.2),
                            color=I("black")) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                size = -log10(delta_q_value),
                                color=delta_log_abund)) +
      ggnetwork::geom_nodetext(data = g,
                               aes(x, y,
                                   label = q_value_sig_code),
                               color=I("black"))
  }
  
  
  # ggnetwork::geom_nodetext_repel(data = g,
  #                                aes(x, y,
  #                                    label = q_value_sig_code),
  #                                color=I("black"))
  #labs(fill = color_nodes_by)

  p = p + scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")

  if (is.null(edge_labels) == FALSE) {
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  #p = p +
  #  ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))

  p = p + facet_wrap(~contrast)

  p = p + scale_size(range=c(1, node_size)) +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position) + guides(fill = "none", size="none")
  return(p)
}

#' Plot the expression of a single gene on the state transition graph
#' @export
plot_state_graph_gene_expression <- function(ccs,
                                             state_graph,
                                             genes,
                                             method = "min",
                                             gene_expr = NULL,
                                             fract_expr=0.0,
                                             mean_expr=0.0,
                                             scale_to_range = FALSE,
                                             color_nodes_by=NULL,
                                             label_nodes_by=NULL,
                                             group_nodes_by=NULL,
                                             label_edges_by=NULL,
                                             edge_weights=NULL,
                                             arrow.gap=0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             min_node_size = 0.25,
                                             max_node_size = 2,
                                             min_edge_size=0.1,
                                             max_edge_size=2,
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

  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by=NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
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

  gene_ids = rowData(ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% genes) %>% rownames()

  if (is.null(gene_expr)) {
    gene_expr = hooke:::aggregated_expr_data(ccs@cds[gene_ids,], group_cells_by = ccs@info$cell_group)
  }

  sub_gene_expr = gene_expr %>%
    filter(gene_short_name %in% genes)

  method = get(method)
  # sub_gene_expr = sub_gene_expr %>%
  #   group_by(cell_group) %>%
  #   dplyr::summarise(fraction_expressing = method(fraction_expressing),
  #                    mean_expression = method(mean_expression),
  #                    specificity = method(specificity)) %>%
  #   mutate(gene_short_name =  paste(genes, collapse = "-"))


  #node_metadata = node_metadata %>% left_join(sub_gene_expr, by = c("id" = "cell_group"))

  sub_gene_expr = sub_gene_expr %>%
    group_by(gene_short_name) %>%
    mutate(
      max_expr = max(mean_expression),
      fraction_max = ifelse (max_expr > 0, mean_expression / max_expr, 0),
      gene_expr = case_when(
        fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
        TRUE ~ FALSE))

  color_nodes_by = "mean_expression"
  expression_legend_label = "mean expression"

  if (scale_to_range) {
    sub_gene_expr = sub_gene_expr %>%
      mutate(value = mean_expression) %>%
      group_by(gene_short_name) %>%
      dplyr::mutate(max_val_for_feature = max(value),
                    min_val_for_feature = min(value)) %>%
      dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
    expression_legend_label = "% Max"
    color_nodes_by = "value"
  }

  g = left_join(g, sub_gene_expr, by=c("name"="cell_group"))

  g$gene_short_name = factor(g$gene_short_name, levels=genes)

  # group_outline = TRUE

  p <- ggplot(aes(x,y), data=g)
  p = p +
    ggplot2::geom_path(aes(x, y, group=edge_name), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')


  # if (is.null(label_subset)) {
  #   label_subset = unique(node_metadata$group_nodes_by)
  # }
  #
  # label_subset = label_subset[label_subset != unlabeled_groups]
  group_nodes_by = NULL
  if (is.null(group_nodes_by) == FALSE){
    p = p + ggforce::geom_mark_rect(aes(x, y,
                                        fill = gene_short_name,
                                        color = group_nodes_by,
                                        filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                    size=0.5,
                                    expand = unit(2, "mm"),
                                    label.buffer=unit(1, "mm"),
                                    radius = unit(1.5, "mm"),
                                    label.margin = margin(1, 1, 1, 1, "mm"),
                                    label.fontsize=group_label_font_size,
                                    label.fontface="plain",
                                    con.linetype=label_conn_linetype,
                                    con.colour=con_colour,
                                    show.legend = F) #+
      # ggforce::geom_mark_rect(aes(x, y,
      #                             fill = gene_short_name,
      #                             color = group_nodes_by,
      #                             label = group_nodes_by,
      #                             filter = group_nodes_by %in% label_subset),
      #                         size=0.5,
      #                         expand = unit(2, "mm"),
      #                         label.buffer=unit(1, "mm"),
      #                         radius = unit(1.5, "mm"),
      #                         label.margin = margin(1, 1, 1, 1, "mm"),
      #                         label.fontsize=group_label_font_size,
      #                         label.fontface="plain",
      #                         con.linetype=label_conn_linetype,
      #                         con.colour=con_colour,
      #                         show.legend = F)

    p = p + scale_fill_manual(values=rep("white", length(unique(g$gene_short_name))))

  }
  p = p + guides(fill = "none")
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr),
                          aes(x, y,
                              size = fraction_max,
                              color = I(con_colour))) +
    ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr & fraction_max > 0),
                          aes(x, y,
                              size = fraction_max,
                              color = get(color_nodes_by))) +
    labs(color = expression_legend_label) +
    viridis::scale_color_viridis(option = "viridis")

  if (is.null(edge_labels) == FALSE) {
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  p = p + facet_wrap(~gene_short_name)

  p = p + 
    scale_size_continuous(labels = scales::percent, range=c(min_node_size, max_node_size)) +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position) + guides(fill = "none")
  return(p)
}

#' @noRd
collect_psg_node_metadata <- function(ccs,
                                      color_nodes_by,
                                      label_nodes_by,
                                      group_nodes_by)
{
  cell_groups = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  metadata_cols = c(color_nodes_by,
                    group_nodes_by)
  if (is.null(label_nodes_by) == FALSE && label_nodes_by != "cell_group")
    metadata_cols = c(metadata_cols, label_nodes_by)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = ccs@cds_coldata[,metadata_cols, drop=F] %>%
    as.data.frame
  cell_group_metadata$cell_group = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  if (is.null(color_nodes_by) == FALSE){
    color_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      dplyr::count(cell_group, !!sym(color_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(color_by_metadata) = c("cell_group", "color_nodes_by")
    node_metadata = left_join(node_metadata, color_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(group_nodes_by) == FALSE){
    group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
      as.data.frame %>%
      dplyr::count(cell_group, !!sym(group_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
    node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(label_nodes_by) == FALSE){
    label_by_metadata = cell_group_metadata[,c("cell_group", label_nodes_by), drop=F]
    colnames(label_by_metadata) = c("cell_group", "label_nodes_by")
    label_by_metadata = label_by_metadata %>%
      as.data.frame %>%
      dplyr::count(cell_group, label_nodes_by) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    node_metadata = left_join(node_metadata, label_by_metadata, by=c("id"="cell_group"))
  }else{
    node_metadata$label_nodes_by = node_metadata$id
  }

  node_metadata = node_metadata %>% distinct() %>% as.data.frame(stringsAsFactor=FALSE)
  row.names(node_metadata) = node_metadata$id

  return(node_metadata)
}


#' @noRd
layout_state_graph <- function(G, node_metadata, edge_labels, weighted=FALSE)
{

  if (weighted){
    G_nel = graph::graphAM(igraph::get.adjacency(G, attr="weight") %>% as.matrix(),
                           edgemode = 'directed', values=list(weight=1)) %>% as("graphNEL")
  } else{
    G_nel = graph::graphAM(igraph::get.adjacency(G) %>% as.matrix(),
                           edgemode = 'directed') %>%as("graphNEL")
  }

  edge_weights = unlist(graph::edgeWeights(G_nel))
  names(edge_weights) = stringr::str_replace_all(names(edge_weights), "\\.", "~")

  # make_subgraphs_for_groups <- function(grouping_set, G_nel, node_meta_df){
  #   nodes = node_meta_df %>% filter(group_nodes_by == grouping_set) %>% pull(id) %>% as.character %>% unique()
  #   #sg = list(graph=graph::subGraph(snodes=nodes, graph=G_nel))
  #   sg = graph::subGraph(snodes=nodes, graph=G_nel)
  #   return (sg)
  # }

  make_subgraphs_for_groups <- function(subgraph_ids, G_nel){
    nodes = subgraph_ids %>% pull(id) %>% as.character %>% unique()
    sg = list(graph=graph::subGraph(snodes=nodes, graph=G_nel), cluster=TRUE)
    #sg = graph::subGraph(snodes=nodes, graph=G_nel)
    return (sg)
  }
  
  if (is.null(node_metadata)){
    subgraphs = NULL
  } else {
    subgraph_df = node_metadata %>%
      select(group_nodes_by, id) %>%
      group_by(group_nodes_by) %>%
      tidyr::nest(subgraph_ids = id) %>%
      summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
                                      .x = subgraph_ids,
                                      G_nel))
    subgraphs = subgraph_df$subgraph
    names(subgraphs) = subgraph_df$group_nodes_by
  }
  
    #summarize(subgraph = graph::subGraph(id,graph=G_nel)
    # summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
    #                                 .x = group_nodes_by,
    #                                 G_nel,
    #                                 node_metadata))

  if (is.null(edge_labels)== FALSE) {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraphs, edgeAttrs=list(label=edge_labels), recipEdges="distinct")
    label_df = data.frame("x" = gvizl@renderInfo@edges$labelX, "y" = gvizl@renderInfo@edges$labelY) %>%
      tibble::rownames_to_column("edge_name") %>%
      left_join(tibble("edge_name" = names(edge_labels), label=edge_labels))

  } else {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraphs, recipEdges="distinct")
    label_df=NULL
  }

  gvizl_coords = cbind(gvizl@renderInfo@nodes$nodeX, gvizl@renderInfo@nodes$nodeY)

  beziers = lapply(gvizl@renderInfo@edges$splines, function(bc) {
    bc_segments = lapply(bc, Rgraphviz::bezierPoints)
    bezier_cp_df = do.call(rbind, bc_segments) %>% as.data.frame
    colnames(bezier_cp_df) = c("x", "y")
    #bezier_cp_df$point = "control"
    #bezier_cp_df$point[1] = "end"
    #bezier_cp_df$point[nrow(bezier_cp_df)] = "end"
    bezier_cp_df
    #control_point_coords = lapply(bc, function(cp) Rgraphviz::getPoints)
    #control_point_coords = rbind(control_point_coords)
  })
  bezier_df = do.call(rbind, beziers)
  bezier_df$edge_name = stringr::str_split_fixed(row.names(bezier_df), "\\.", 2)[,1]
  bezier_df$from = stringr::str_split_fixed(bezier_df$edge_name, "~", 2)[,1]
  bezier_df$to = stringr::str_split_fixed(bezier_df$edge_name, "~", 2)[,2]
  #bezier_df = bezier_df %>% mutate(x = ggnetwork:::scale_safely(x),
  #                                 y = ggnetwork:::scale_safely(y))

  bezier_df = left_join(bezier_df, tibble(edge_name=names(gvizl@renderInfo@edges$direction), edge_direction=gvizl@renderInfo@edges$direction))
  bezier_df = bezier_df %>% dplyr::distinct()
  return(list(gvizl_coords=gvizl_coords, bezier_df=bezier_df, label_df=label_df))
}



plot_state_graph_annotations_wrapper = function(ccm, state_graph, edge_weights="support") {
  plot_state_graph_annotations(ccm,
                               state_graph,
                               color_nodes_by = "timepoint",
                               group_nodes_by = "cell_type_sub",
                               edge_weights = edge_weights,
                               hide_unlinked_nodes = FALSE)
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

  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

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

  layout_info = layout_state_graph(G, node_metadata, edge_labels, weighted=FALSE)
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
  # add assembly group
  bezier_df$assembly_group = vapply(strsplit(bezier_df$from, "-", fixed=T), "[", "", 1)
  g$assembly_group = vapply(strsplit(g$name, "-", fixed=T), "[", "", 1)
  
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

  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

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
  # add assembly group
  bezier_df$assembly_group = vapply(strsplit(bezier_df$from, "-", fixed=T), "[", "", 1)
  g$assembly_group = vapply(strsplit(g$name, "-", fixed=T), "[", "", 1)
  
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

  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

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
  # add assembly group
  bezier_df$assembly_group = vapply(strsplit(bezier_df$from, "-", fixed=T), "[", "", 1)
  g$assembly_group = vapply(strsplit(g$name, "-", fixed=T), "[", "", 1)
  
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



# cell_type_ann_df = data.frame(cell_group = unique(colData(cds)[[group_cells_by]]), 
#                               category = rep(1:7, each=2)) %>% 
#   tibble::column_to_rownames("cell_group")

#'
#' @param cds description
#' @param group_cells_by description
#' @param genes description
#' @param cell_types
#' @param gene_ann_df one column must match genes and be called "genes" 
#' @param cell_type_ann_df one column must match cell types and be called "cell_group" 
plot_genes_expr_across_cell_groups = function(cds, 
                                              group_cells_by, 
                                              genes = NULL, 
                                              cell_types = NULL, 
                                              gene_ann_df = NULL, 
                                              cell_type_ann_df = NULL, 
                                              drop_zeroes = TRUE) {
  
  if (is.null(cell_types) == FALSE) {
    cds = cds[,colData(cds)[[group_cells_by]] %in% cell_types]
  }
  
  # maybe warn if this isnt null / longer than a certain length 
  if (is.null(genes) == FALSE) {
    cds = cds[rowData(cds)$gene_short_name %in% genes, ]
  }
  
  agg_expr_data = hooke:::aggregated_expr_data(cds, group_cells_by = group_cells_by)
  
  
  if (is.null(gene_ann_df)==FALSE) {
    gene_ann_df = gene_ann_df %>% tibble::column_to_rownames("genes")
  }
  
  
  if (is.null(cell_type_ann_df)==FALSE) {
    cell_type_ann_df = cell_type_ann_df %>% tibble::column_to_rownames("cell_group")
  }
  
  agg_expr_mat = agg_expr_data %>% 
    filter(gene_short_name %in% genes) %>% 
    select(cell_group, gene_short_name, mean_expression) %>% 
    pivot_wider(names_from = gene_short_name, values_from = mean_expression) %>% 
    tibble::column_to_rownames("cell_group") %>% 
    t()
  
  if (drop_zeroes) {
    agg_expr_mat = agg_expr_mat[rowSums(agg_expr_mat) != 0,
                                colSums(agg_expr_mat) != 0]
  }
  
 
  p = pheatmap::pheatmap(agg_expr_mat, 
                         cluster_rows = T, 
                         cluster_cols = T, 
                         color = viridis::viridis(10), 
                         annotation_row = gene_ann_df, 
                         annotation_col = cell_type_ann_df)
  return(p)
  
  
}


#' selects the fold changes from when a cell type is at peak abundance in the control
#' @param ccm cell_count_model, assumes output of fit_genotype_ccm 
#' @param start_time
#' @param stop_time
#' @param ctrl_abundances you can put in precomputed estimates but 
get_max_abundance_contrast = function(ccm, 
                                      start_time = NULL, 
                                      stop_time = NULL, 
                                      ctrl_abundances = NULL, 
                                      interval_col = "timepoint", 
                                      newdata = tibble()){
  
  timepoints = as.numeric(unique(colData(ccm@ccs@cds)[[interval_col]]))
  timepoints = timepoints[!is.na(timepoints)]
  
  if (is.null(start_time)) {
    start_time = min(timepoints)
  }
  if (is.null(stop_time)) {
    stop_time = max(timepoints)
  }
  
  if (nrow(newdata) > 0 ){
    newdata = cross_join(tibble(knockout=FALSE), newdata)
  } else {
    newdata = tibble(knockout=FALSE) 
  }
  
  # CHECK THIS 
  perturb_effects = platt:::get_perturbation_effects(ccm, ...)
  
  
  if (is.null(ctrl_abundances)) {
    ctrl_abundances = get_extant_cell_types(ccm, start = start_time, stop = stop_time, newdata = newdata) 
  } else {
    ctrl_abundances = ctrl_abundances %>% filter(timepoint %in% unique(perturb_effects$time))
  }
  
  
  peak_abundances = ctrl_abundances %>% 
    filter(!is.na(percent_cell_type_range)) %>% 
    group_by(cell_group) %>% 
    top_n(n = 1, wt = percent_max_abund) %>% 
    select(timepoint, cell_group) 
  
  peak_abundances$timepoint = as.character(peak_abundances$timepoint)
  
  max_abund_tbl = left_join(peak_abundances, perturb_effects, by = c("cell_group", "timepoint" = "time"))
  
  return(max_abund_tbl)
}




plot_graph_simple = function( ccs, 
                              state_graph, 
                              num_degs, 
                              color_nodes_by = NULL, 
                              label_nodes_by= NULL, 
                              group_nodes_by= NULL, 
                              arrow.gap=0.03,
                              arrow_unit = 2,
                              bar_unit = .075,
                              node_size = 2,
                              min_edge_size=0.1,
                              max_edge_size=2, 
                              con_colour = "darkgrey", 
                              legend_position = "none", 
                              hide_unlinked_nodes=TRUE, 
                              flip_x=F) {
  
  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }
  
  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)
  
  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }
  
  edges = edges %>% dplyr::ungroup()
  edges = edges %>% filter(from %in% node_metadata$id & to %in% node_metadata$id)
  
  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)
  
  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  
  

  
  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)
  g = left_join(g, num_degs, by=c("name"="cell_group"))
  
  if (flip_x) {
    g$x = -1*g$x
    bezier_df$x = -1*bezier_df$x
  }
  
  bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
  bezier_df$unsupported_edge = FALSE
  
  g = g %>% filter(!is.na(term))
  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), 
                       colour=con_colour, data=bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')
  
  p <- p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y,
                              color = log10(n),
                              size = node_size)) +
    labs(color = "number degs") + 
    scale_color_viridis_c()
  
  p <- p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position)
  
  p = p + facet_wrap(~term)
  return(p)
  
}
