# to do : add edge stuff and node labels

plot_annotations = function(cell_state_graph, 
                            color_nodes_by = NULL, 
                            label_nodes_by = NULL, 
                            arrow_unit = 7,
                            node_size = 2,
                            con_colour = "darkgrey", 
                            legend_position = "none", 
                            min_edge_size=0.1,
                            max_edge_size=2,
                            edge_weights=NULL, 
                            plot_labels=T) {
  
  if (is.null(color_nodes_by)){
    color_nodes_by = cell_state_graph@ccs@info$cell_group
  } else {
    cell_state_graph@g[["color_nodes_by"]] = color_nodes_by
  }
  
  group_nodes_by = cell_state_graph@metadata$group_nodes_by
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
 
  
  p <- ggplot(aes(x,y), data=g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name, linetype=unsupported_edge), 
                       colour = con_colour, 
                       data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  
  if (is.null(group_nodes_by) == FALSE) {
    
    p =  p + 
      ggnetwork::geom_nodes(data = g,
                            aes(x, y) ,
                            color=I("black"), size=node_size*1.2) +
      # ggnetwork::geom_nodes(data = g,
      #                       aes(x, y,
      #                           color = color_nodes_by),
      #                       size = node_size)  +
      labs(fill = color_nodes_by)
    
    p = p + ggforce::geom_mark_rect(aes(x, y,
                                        fill = group_nodes_by),
                                    size=0,
                                    expand = unit(2, "mm"),
                                    radius = unit(1.5, "mm"),
                                    data=g)
    
  } else {
    p = p + 
      ggnetwork::geom_nodes(data = g,
                            aes(x, y) ,
                            color=I("black"), size=node_size*1.2) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                color = color_nodes_by),
                            size = node_size)  +
      labs(fill = color_nodes_by)
    
  }
  
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black"), 
                                     box.padding = 0.5) 
  }
  
  
  p = p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    scale_color_manual(values = hooke:::get_colors(length(unique(g$name)), type="vibrant")) +
    theme(legend.position=legend_position)
  
  # plot an invisible point to help with nodes not being cut off
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
             mapping = aes(x,y), color="white", alpha=0)
  


  
  
  return(p)
}

#' 
#' @param cell_state_graph
#' @param comp_abund_table
#' @param facet_group

plot_abundance_changes = function(cell_state_graph, 
                                  comp_abund_table,
                                  facet_group = NULL,
                                  scale_node = FALSE,
                                  plot_labels = FALSE, 
                                  arrow_unit = 7,
                                  node_size = 2,
                                  node_scale = 1,
                                  con_colour = "darkgrey", 
                                  legend_position = "none", 
                                  fc_limits = c(-3,3)) {
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  
  if (is.null(facet_group) == FALSE) {
    comp_abund_table[["contrast"]] = comp_abund_table[[facet_group]]
    comp_abund_table$contrast = as.factor(comp_abund_table$contrast)
  }
  min = fc_limits[1]
  max = fc_limits[2]
  comp_abund_table = comp_abund_table %>%
    mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
    mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  
  comp_abund_table = comp_abund_table %>%
    mutate(
      delta_q_value = pmax(0.0001, delta_q_value),
      q_value_sig_code = platt:::calc_sig_ind(delta_q_value, html=FALSE))
  
  g = left_join(g, comp_abund_table, by=c("name"="cell_group"), relationship = "many-to-many")
  
  p <- ggplot(aes(x,y), data = g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour = con_colour, data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  
  if (scale_node) {
    
    p = p + ggnewscale::new_scale_color() +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y, 
                                size = -log10(delta_q_value)*1.2 * node_scale), 
                            color=I("black")) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                size = -log10(delta_q_value) * node_scale,
                                color=delta_log_abund)) +
      ggnetwork::geom_nodetext(data = g,
                               aes(x, y,
                                   label = q_value_sig_code),
                               color=I("black")) +  
      scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3") + 
      scale_size_identity() +
      ggnetwork::theme_blank() +
      hooke_theme_opts() +
      theme(legend.position=legend_position)
    
    
    
  } else {
    
    p = p + ggnewscale::new_scale_color() +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y) ,
                            color=I("black"), size=node_size*1.2) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                color=delta_log_abund), size=node_size) +
      ggnetwork::geom_nodetext(data = g,
                               aes(x, y,
                                   label = q_value_sig_code),
                               color=I("black")) +
      scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3") + 
      # scale_size(range=c(2, 6)) +
      scale_size_identity() +
      ggnetwork::theme_blank() +
      hooke_theme_opts() +
      theme(legend.position=legend_position) 
    
  }
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                             aes(x, y, label=name),
                             color=I("black")) 
  }
  
  
  if (is.null(facet_group) == FALSE) {
    p = p + facet_wrap(~contrast)
  }
  
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
}


# given a list of genes will plot the sum, mean, min or max expression of those genes
# on top of a cell_state_graph
#' @param cell_state_graph
#' @param genes
#' @param color_nodes_by
#' 
plot_gene_expr = function(cell_state_graph, 
                          genes, 
                          color_nodes_by = c("sum_expr", "mean_expr", "min_expr", "max_expr"),
                          arrow_unit = 7,
                          node_size = 2,
                          con_colour = "darkgrey",
                          fract_expr = 0.0,
                          mean_expr = 0.0,
                          legend_position = "none", 
                          plot_labels = F) {
  
  g = cell_state_graph@g
  ccs = cell_state_graph@ccs
  bezier_df = cell_state_graph@layout_info$bezier_df
  color_nodes_by = match.arg(color_nodes_by)
  
  gene_ids = rowData(ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% genes) %>% rownames()
  gene_expr = hooke:::aggregated_expr_data(ccs@cds[gene_ids,], group_cells_by = ccs@info$cell_group)
  sub_gene_expr = gene_expr %>%
    group_by(gene_short_name) %>%
    mutate(
      max_expr = max(mean_expression),
      fraction_max = ifelse (max_expr > 0, mean_expression / max_expr, 0),
      gene_expr = case_when(
        fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
        TRUE ~ FALSE)) 
  scale_to_range = T
  if (scale_to_range) {
    sub_gene_expr = sub_gene_expr %>%
      mutate(value = mean_expression) %>%
      group_by(gene_short_name) %>%
      dplyr::mutate(max_val_for_feature = max(value),
                    min_val_for_feature = min(value)) %>%
      dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
  }
  
  gene_expr_summary =  sub_gene_expr %>% 
    group_by(cell_group) %>% 
    summarize( sum_expr = sum(mean_expression), 
               mean_expr = mean(mean_expression), 
               min_expr = min(mean_expression), 
               max_expr = max(mean_expression), 
               fraction_max = sum(fraction_max),
               gene_expr = (min(gene_expr) == 1))
  
  g = left_join(g, gene_expr_summary, by = c("name" = "cell_group"), relationship = "many-to-many")
  
  p <- ggplot(aes(x,y), data=g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour=con_colour, data=bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr),
                          aes(x, y,
                              size = fraction_max * node_size * 1.2), 
                          color=I("black")) +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr),
                          aes(x, y,
                              size = fraction_max * node_size,
                              color = I(con_colour))) +
    ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr & fraction_max > 0),
                          aes(x, y,
                              size = fraction_max * node_size,
                              color = sum_expr)) +
    labs(color = color_nodes_by) +
    viridis::scale_color_viridis(option = "viridis") + 
    ggnetwork::theme_blank() +
    scale_size_identity() +
    scale_size(range=c(1, 5)) + 
    hooke_theme_opts() + theme(legend.position = "none") 
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black")) 
  }
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  return(p)
}


plot_deviation_plot = function(cell_state_graph, 
                               deviation_table, 
                               facet_group = NULL, 
                               arrow_unit = 7,
                               node_size = 2,
                               con_colour = "darkgrey",
                               fract_expr = 0.0,
                               mean_expr = 0.0,
                               legend_position = "none", 
                               plot_labels=F){
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  
  g = left_join(g, deviation_table, by=c("name"="cell_group"), relationship = "many-to-many")
  
  if (is.null(facet_group) == FALSE) {
    g[["contrast"]] = g[[facet_group]]  
    g$contrast = as.factor(g$contrast)
  }
 
  
  p <- ggplot(aes(x,y), data = g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour = con_colour, data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  # deviation plot 
  p = p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y,
                              size = node_size * 1.2), 
                          color=I("black")) +
    ggnetwork::geom_nodes(data = g %>% filter(!is.na(deviation)),
                          aes(x, y,
                              color = log10(deviation+1),
                              size = node_size)) +
    ggnetwork::geom_nodes(data = g %>% filter(is.na(deviation)),
                          aes(x, y,
                              size = node_size), color="white") +
    labs(color = "number degs") + 
    scale_color_viridis_c() + 
    scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position='none') 
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black")) 
  }
  
  if (is.null(facet_group) == FALSE) {
    p = p + facet_wrap(~contrast)
  }
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
  
}


plot_deg_change = function(cell_state_graph, 
                          deg_table, 
                          facet_group = "term", 
                          arrow_unit = 7,
                          node_size = 2,
                          con_colour = "darkgrey",
                          fract_expr = 0.0,
                          mean_expr = 0.0,
                          legend_position = "none", 
                          fc_limits = c(-3,3)){ 
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  
  # g = left_join(g, deg_table, by = c("name"="cell_group", "term"), relationship = "many-to-many") 
  
  g = left_join(g, deg_table, by = c("name"="cell_group"), relationship = "many-to-many") 
  
  # myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  # sc <- scale_colour_gradientn(colours = myPalette(100), limits=fc_limits)
  
  p <- ggplot(aes(x,y), data = g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour = con_colour, data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  p =  p + 
    ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y) ,
                          color=I("black"), size=node_size*1.2) +
    ggnetwork::geom_nodes(data = g ,
                          mapping = aes(x, y, color = n),
                          size = node_size) + 
    ggnetwork::theme_blank() + 
    # sc + 
    scale_color_viridis_c(option="B")+
    hooke_theme_opts() +
    theme(legend.position='none')
  
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
  
  
    
}

#' this plotting function 
#' 
plot_degs = function(cell_state_graph, 
                     deg_table, 
                     facet_group = "term", 
                     arrow_unit = 7,
                     node_size = 2,
                     con_colour = "darkgrey",
                     fract_expr = 0.0,
                     mean_expr = 0.0,
                     legend_position = "none", 
                     fc_limits = c(-3,3)){ 
  
  # deg_table = deg_table %>%
  #   mutate(
  #     delta_q_value = pmax(0.0001, perturb_to_ctrl_p_value),
  #     p_value_sig_code = calc_sig_ind(perturb_to_ctrl_p_value, html=FALSE)) 
  # 
  # deg_table[["contrast"]] = deg_table[[facet_group]]
  # deg_table$contrast = as.factor(deg_table$contrast)
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  
  # g = left_join(g, deg_table, by = c("name"="cell_group", "term"), relationship = "many-to-many") 
  
  g = left_join(g, deg_table, by = c("name"="cell_group"), relationship = "many-to-many") 
  
  # myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  # sc <- scale_colour_gradientn(colours = myPalette(100), limits=fc_limits)
  
  p <- ggplot(aes(x,y), data = g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour = con_colour, data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  p =  p + 
    ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y) ,
                          color=I("black"), size=node_size*1.2) +
    ggnetwork::geom_nodes(data = g ,
                          mapping = aes(x, y, color = n),
                          size = node_size) + 
    ggnetwork::theme_blank() + 
    # sc + 
    scale_color_viridis_c(option="B")+
    hooke_theme_opts() +
    theme(legend.position='none')
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
  
}
