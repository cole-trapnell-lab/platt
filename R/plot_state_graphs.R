#' plots a cell state graph
#' @param cell_state_graph
#' @param color_cells_by
#' @export
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

#' plots the cell abundance changes over a platt graph
#' @param cell_state_graph
#' @param comp_abund_table
#' @param facet_group
#' @export
plot_abundance_changes = function(cell_state_graph, 
                                  comp_abund_table,
                                  facet_group = NULL,
                                  scale_node = FALSE,
                                  plot_labels = TRUE, 
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
      scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3", limits = fc_limits) + 
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
      scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3", limits = fc_limits) +  
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
  
  p = p + guides(color=guide_colourbar(title="log(\u0394 Abundance)"))
  
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
#' @param node_size Does not actually control the sizes of nodes, but is used to specify the offset to use for invisible points that help prevent the plot from getting clipped upon saving
#' @export
plot_gene_expr = function(cell_state_graph, 
                          genes, 
                          color_nodes_by = c("sum_expr", "mean_expr", "min_expr", "max_expr"),
                          arrow_unit = 7,
                          node_size = 2,
                          con_colour = "lightgrey",
                          fract_expr = 0.0,
                          mean_expr = 0.0,
                          legend_position = "right", 
                          plot_labels = T, 
                          aggregate = F, 
                          scale_to_range = T, 
                          expr_limits = NULL) {
  
  if (scale_to_range && aggregate) {
    message("Warning: scale_to_range is not compatible with aggregate. Setting scale_to_range to FALSE.")
    scale_to_range = FALSE
  }
  
  g = cell_state_graph@g
  ccs = cell_state_graph@ccs
  bezier_df = cell_state_graph@layout_info$bezier_df
  color_nodes_by = match.arg(color_nodes_by)
  
  gene_info = rowData(ccs@cds) %>%
    as.data.frame() %>%
    filter(gene_short_name %in% genes) %>%
    select(gene_short_name)
  duplicated_genes = gene_info %>%
    summarize(n = n(), .by = gene_short_name) %>%
    filter(n > 1) %>%
    pull(gene_short_name)
  if (aggregate == FALSE && length(duplicated_genes) > 0) {
    message(
      "The following gene names are represented by multiple ENSEMBL IDs: ",
      paste(duplicated_genes, collapse = ", "),
      ". Each ENSEMBL ID will be represented separately in the plot."
    )
    gene_info = gene_info %>%
      mutate(gene_short_name = ave(gene_short_name, gene_short_name, FUN = function(x) {
        if (length(x) > 1) {
          paste0(x, "_", seq_along(x))
        } else {
          x
        }
      }))
  }
  gene_expr = hooke:::aggregated_expr_data(ccs@cds[rownames(gene_info), ], group_cells_by = ccs@info$cell_group)
  sub_gene_expr = gene_expr %>%
    mutate(
      max_expr = max(mean_expression),
      min_expr = min(mean_expression),
      fraction_max = ifelse(max_expr > 0, mean_expression / max_expr, 0),
      gene_expr = case_when(
        fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
        TRUE ~ FALSE
      ),
      .by = gene_id
    )
  
  if (aggregate == FALSE) {
    gene_expr_summary = sub_gene_expr %>%
      mutate(gene_short_name = gene_info[gene_id, "gene_short_name"]) %>%
      dplyr::rename(sum_expr = mean_expression)
    if (scale_to_range) {
      gene_expr_summary = gene_expr_summary %>%
        mutate(
          sum_expr = (sum_expr - min_expr) / (max_expr - min_expr),
          .by = gene_id
        )
    }
  } else {
    gene_expr_summary =  sub_gene_expr %>% 
      group_by(cell_group) %>% 
      summarize( sum_expr = sum(mean_expression), 
                #  mean_expr = mean(mean_expression), 
                #  min_expr = min(mean_expression), 
                #  max_expr = max(mean_expression), 
                 fraction_max = mean(fraction_max),
                 gene_expr = (min(gene_expr) == 1)) %>% 
      mutate(gene_short_name = paste(sort(genes), collapse =";"))
  }
  
  
  
  if (is.null(expr_limits) == FALSE){
    min = expr_limits[1]
    max = expr_limits[2]
    gene_expr_summary = gene_expr_summary %>%
      mutate(sum_expr = ifelse(sum_expr > max, max, sum_expr)) %>%
      mutate(sum_expr = ifelse(sum_expr < min, min, sum_expr))
  }
  
  g = left_join(g, gene_expr_summary, by = c("name" = "cell_group"), relationship = "many-to-many")
  
  
  
  p <- ggplot(aes(x,y), data=g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour=con_colour, data=bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')
  
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr),
                          aes(x, y,
                              size = fraction_max,
                          ),
                          shape = "circle filled",
                          fill = I(con_colour),
                          color = I("black")
      ) +
    ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr & fraction_max > 0),
                          aes(x, y,
                              size = fraction_max,
                              fill = sum_expr
                            ),
                            shape = "circle filled",
                            color = I("black")
      ) +
    labs(color = color_nodes_by) +
    scale_fill_viridis_c(limits = expr_limits) +
    ggnetwork::theme_blank() +
    scale_size_identity() +
    scale_size(range=c(1, 5)) + 
    hooke_theme_opts() + theme(legend.position = legend_position) 
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black")) 
  }
  
  p = p +
    guides(fill = guide_colourbar(title = ifelse(scale_to_range, "Rel. Gene Expr.", "Gene Expr."))) +
    labs(size = "Fract. of Embryos") +
    facet_wrap(~gene_short_name)
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
    scale_color_viridis_c(option="plasma") + 
    scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position) 
  
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
    theme(legend.position=legend_position)
  
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
  
  
    
}

#' this plotting function plots the fold change of a given gene in a deg table
#' over the platt graph
#' @param cell_state_graph 
#' @param deg_table
#' @export
plot_degs = function(cell_state_graph, 
                     deg_table, 
                     facet_group = "term", 
                     arrow_unit = 7,
                     node_size = 2,
                     con_colour = "darkgrey",
                     fract_expr = 0.0,
                     mean_expr = 0.0,
                     legend_position = "none", 
                     fc_limits = c(-3,3), 
                     plot_labels=T){ 
  
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
  
  if (is.null(fc_limits)==FALSE) {
    min = fc_limits[1]
    max = fc_limits[2]
    deg_table = deg_table %>%
      mutate(perturb_to_ctrl_raw_lfc = ifelse(perturb_to_ctrl_raw_lfc > max, max, perturb_to_ctrl_raw_lfc)) %>%
      mutate(perturb_to_ctrl_raw_lfc = ifelse(perturb_to_ctrl_raw_lfc < min, min, perturb_to_ctrl_raw_lfc))
  }
  
  
  deg_table = deg_table %>% mutate(
    delta_q_value = pmax(0.0001, perturb_to_ctrl_p_value),
    q_value_sig_code = platt:::calc_sig_ind(perturb_to_ctrl_p_value, html=FALSE))
  
  g = left_join(g, deg_table, by = c("name"="cell_group"), relationship = "many-to-many") 
  
  
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
    ggnetwork::geom_nodes(data = g %>% filter(!is.na(perturb_to_ctrl_raw_lfc)),
                          mapping = aes(x, y, color = perturb_to_ctrl_raw_lfc),
                          size = node_size) + 
    ggnetwork::geom_nodes(data = g %>% filter(is.na(perturb_to_ctrl_raw_lfc)),
                          mapping = aes(x, y), color="lightgray",
                          size = node_size) +
    ggnetwork::theme_blank() + 
    scale_color_gradient2(low = "#006600",  mid = "white", high = "#800080", limits = fc_limits) + 
    # scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3", limits = fc_limits) + 
    # scale_color_viridis_c(option="B", limits = fc_limits)+
    ggnetwork::geom_nodetext(data = g,
                             aes(x, y,
                                 label = q_value_sig_code),
                             color=I("black")) + 
    hooke_theme_opts() +
    scale_size_identity() +
    theme(legend.position=legend_position)
  
  p = p + guides(color=guide_colourbar(title="log2(fc)"))
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black"), 
                                     box.padding = 0.5) 
  }
  
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
  
}

#' plot a larger graph
#' @export
plot_by_annotations = function(cell_state_graph, 
                               colname = "projection_group",
                               arrow_unit = 4,
                               node_size = 1,
                               node_scale = 1, 
                               node_label_width = 50,
                               legend_position = "right", 
                               group_label_font_size=1) {
  
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  
  
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
  
  # hidden_node_labels = cell_state_graph@layout_info$hidden_gvizl_coords
  # hidden_node_labels = hidden_node_labels %>% filter(grepl("head", name)) %>% select(x,y,group) %>% distinct()
  # hidden_node_labels = hidden_node_labels %>% group_by(group) %>% summarize(x = mean(x), y = mean(y))
  
  
  p <- ggplot(aes(x,y), data = g) + 
    ggplot2::geom_path(aes(x, y, 
                           group = edge_name),
                       # colour = con_colour,
                       data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=20, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre') 
  
  if (is.null(grouping_df) == FALSE && identical(grouping_df$group_nodes_by, grouping_df$id) == FALSE){
    
    p = p + ggforce::geom_mark_rect(aes(x, y, group=group_nodes_by, color=I("lightgrey")),
                                    size=0.25,
                                    radius = unit(0.5, "mm"),
                                    expand = unit(1, "mm"),
                                    #con.linetype="dotted",
                                    con.type="straight",
                                    con.colour="lightgrey",
                                    con.size=0.25,
                                    con.border="one",
                                    na.rm=TRUE,
                                    data=g)
    
    p = p + geom_text(data = group_label_position_df, aes(x, y, label = group_nodes_by), size = group_label_font_size)
    
  }
  
  
  vibrant.colors =
    c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
      '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')
  
  bright.colors =
    c('#4477AA',
      '#EE6677',
      '#228833',
      '#CCBB44',
      '#66CCEE',
      '#AA3377',
      '#BBBBBB')
  
  
  if (colname == "cell_type") {
    
    num.colors.tissue = igraph::V(cell_state_graph@graph)$name %>% unique() %>% length()
    
    tissue.colors =
      colorRampPalette(c(vibrant.colors,bright.colors))(num.colors.tissue)
    
    names(tissue.colors) = igraph::V(cell_state_graph@graph)$name %>% unique()
    
  }else {
    
    num.colors.tissue = ref_ccs@cds_coldata[[colname]] %>% unique() %>% length()
    
    tissue.colors =
      colorRampPalette(c(vibrant.colors,bright.colors))(num.colors.tissue)
    
    names(tissue.colors) = ref_ccs@colData$projection_group %>% unique()
    
  }
  
  
  
  # p = p + geom_text(data = hidden_node_labels, aes(x, y, label = group), size = 2)
  
  p = p + ggnewscale::new_scale_color() + 
    ggnetwork::geom_nodes(data = g,
                          aes(x, y) ,
                          color=I("black"), size=node_size*1.2) +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y,
                              color = color_nodes_by)) + 
    scale_color_manual(values = tissue.colors) + 
    monocle3:::monocle_theme_opts() 
  
  p = p + 
    scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position) 
  
  return(p)
  
  
}


plot_perturb_effects <- function(cell_state_graph,
                                 num_top_perturbs=3,
                                 num_top_genes=3, 
                                 arrow_unit = 7,
                                 node_size = 2,
                                 con_colour = "darkgrey",
                                 fract_expr = 0.0,
                                 mean_expr = 0.0,
                                 legend_position = "none", 
                                 plot_labels=T) {
  
  
  num_top_genes=3
  node_support_df = igraph::as_data_frame(cell_state_graph@graph, what="vertices") %>%
    tibble::as_tibble() %>% rename(id=name)
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  
  perturbation_table = cell_state_graph@genetic_requirements
  perturbation_table = perturbation_table %>% 
    mutate(name = sapply(strsplit(id, "-"), `[`, 2))
  
  perturbation_table = perturbation_table %>% group_by(id) %>%
    mutate(perturb_display_name = case_when(perturb_effect == "direct" ~ glue::glue("<i style='color:#2B9EB3'>{perturb_name}</i>"),
                                            perturb_effect == "indirect" ~ glue::glue("<i style='color:#FCAB10'>{perturb_name}</i>"),
                                            perturb_effect == "predicted" ~ glue::glue("<i style='color:#F8333C'>{perturb_name}</i>"),
                                            TRUE ~ glue::glue("<i style='color:#44AF69'>{perturb_name}</i>"))) %>%
    summarize(perturb_effect_label = ifelse(n() > num_top_genes,
                                            paste0(c(perturb_display_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "<br>"),
                                            paste0(perturb_display_name, collapse = "<br>")))

  
  g = left_join(g, perturbation_table, by = c("name"="id"), relationship = "many-to-many") 

  
  g = g %>% mutate(label_nodes_by=perturb_effect_label)
  
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
    ggnetwork::geom_nodes(data = g,
                          mapping = aes(x, y, color = color_nodes_by),
                          size = node_size) +
    ggnetwork::theme_blank() + 
    hooke_theme_opts() +
    scale_size_identity() +
    theme(legend.position=legend_position)
  
  p = p + annotate(geom='richtext', x=0.15*max(g$xend), y=0.01*max(g$yend),
                   size=2,
                   fill = NA, label.color = NA,
                   label="<i style='color:#2B9EB3'>direct</i> <i style='color:#FCAB10'>indirect</i> <i style='color:#F8333C'>predicted</i> <i style='color:#44AF69'>other</i> ")
  
  p = p + guides(color=guide_colourbar(title="log2(fc)"))
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black"), 
                                     box.padding = 0.5) 
  }
  
  
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  return(p)
  
  
}