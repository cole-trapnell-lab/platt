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
                            plot_labels=T,
                            group_label_font_size=1, 
                            node_label_width = 50) {
  
  if (is.null(color_nodes_by)){
    color_nodes_by = cell_state_graph@ccs@info$cell_group
  } else {
    cell_state_graph@g[["color_nodes_by"]] = color_nodes_by
  }
  
  group_nodes_by = cell_state_graph@metadata$group_nodes_by
  
  g = cell_state_graph@g

  bezier_df = cell_state_graph@layout_info$bezier_df
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
 
  
  p <- ggplot(aes(x,y), data=g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name, linetype=unsupported_edge), 
                       colour = con_colour, 
                       data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
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
    plot_labels = FALSE
    color_nodes_by = group_nodes_by
  } 
  
  p = p + 
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                fill = color_nodes_by),
                            size = node_size, 
                            shape = "circle filled",
                            color = I("black"))  +
      labs(fill = color_nodes_by)
    
  if (plot_labels) {
      p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                       aes(x, y, label=name),
                                       color=I("black"), 
                                       box.padding = 0.5) 
  }
    

  num.colors = cell_state_graph@g[["color_nodes_by"]] %>% unique() %>% length()
  
  p = p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    scale_fill_manual(values = hooke:::get_colors(num.colors, type="vibrant")) +
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
                                  node_label_width = 50,
                                  group_label_font_size = 1,
                                  fc_limits = c(-3,3)) {
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
  
  
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
    plot_labels = F
    scale_node = T
    
  }
  
  
  
  if (scale_node) {
    
    p = p + ggnewscale::new_scale_color() +
      # ggnetwork::geom_nodes(data = g,
      #                       aes(x, y, 
      #                           size = -log10(delta_q_value)*1.2 * node_scale), 
      #                       color=I("black")) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                size = -log10(delta_q_value) * node_scale,
                                fill=delta_log_abund), 
                            shape = "circle filled",
                            color = I("black")) +
      ggnetwork::geom_nodetext(data = g,
                               aes(x, y,
                                   label = q_value_sig_code),
                               color=I("black")) +  
      scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3", limits = fc_limits) + 
      scale_size_identity() +
      ggnetwork::theme_blank() +
      hooke_theme_opts() +
      theme(legend.position=legend_position)
    
    
    
  } else {
    
    p = p + ggnewscale::new_scale_color() +
      # ggnetwork::geom_nodes(data = g,
      #                       aes(x, y) ,
      #                       color=I("black"), size=node_size*1.2) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                fill=delta_log_abund), 
                            shape = "circle filled",
                            color = I("black"),
                            size=node_size) +
      ggnetwork::geom_nodetext(data = g,
                               aes(x, y,
                                   label = q_value_sig_code),
                               color=I("black")) +
      scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3", limits = fc_limits) +  
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
  
  p = p + guides(fill=guide_colourbar(title="log(\u0394 Abundance)"))
  
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
#' @param node_size Does not actually control the sizes of nodes, but is used to specify the offset to use for invisible points that help prevent the plot from getting clipped upon saving
#' @param log_expr If TRUE, the expression values will be log-transformed
#' @export
plot_gene_expr = function(cell_state_graph, 
                          genes, 
                          arrow_unit = 7,
                          node_size = 2,
                          con_colour = "lightgrey",
                          fract_expr = 0.0,
                          mean_expr = 0.0,
                          legend_position = "right", 
                          plot_labels = TRUE, 
                          aggregate = FALSE, 
                          scale_to_range = FALSE,
                          log_expr = FALSE,
                          pseudocount = 1e-5,
                          expr_limits = NULL, 
                          node_label_width = 50,
                          group_label_font_size=1) {
  
  if (scale_to_range && aggregate) {
    message("Warning: scale_to_range is not compatible with aggregate. Setting scale_to_range to FALSE.")
    scale_to_range = FALSE
  }
  
  g = cell_state_graph@g
  ccs = cell_state_graph@ccs
  bezier_df = cell_state_graph@layout_info$bezier_df
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
  
  
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
      is_expr = (fraction_expressing >= fract_expr & mean_expression >= mean_expr),
      .by = gene_id
    )
  if (aggregate) {
    gene_expr_summary = sub_gene_expr %>%
      group_by(cell_group) %>%
      summarize(
        sum_expr = sum(mean_expression),
        fraction_expressing = mean(fraction_expressing),
        is_expr = all(is_expr)
      ) %>%
      mutate(gene_short_name = paste(sort(genes), collapse = ";"))
  } else {
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
  }
  
  if (log_expr) {
    gene_expr_summary = gene_expr_summary %>%
      mutate(sum_expr = log10(sum_expr + pseudocount))
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
    plot_labels = FALSE
  }
  
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(is_expr),
                          aes(x, y,
                              size = fraction_expressing,
                              # size = fraction_max,
                          ),
                          shape = "circle filled",
                          fill = I(con_colour),
                          color = I("black")
      ) +
    ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g %>% filter(is_expr & fraction_expressing > 0),
                          aes(x, y,
                              size = fraction_expressing,
                              fill = sum_expr
                            ),
                            shape = "circle filled",
                            color = I("black")
      ) +
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
  
  fill_title = "Gene Expr."
  if (scale_to_range) fill_title = paste("Rel.", fill_title)
  if (log_expr) fill_title = paste("log10(", fill_title, "+", pseudocount, ")")
  p = p +
    guides(fill = guide_colourbar(title = fill_title)) +
    labs(size = "Fract. of Cells") +
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
#' @param perturb_table output of Hooke compare_abundances, if provide will scale the node size by delta_log_abundance
#' @export
plot_degs = function(cell_state_graph, 
                     deg_table, 
                     perturb_table = NULL,
                     facet_group = "term", 
                     arrow_unit = 7,
                     node_size = 1,
                     con_colour = "darkgrey",
                     fract_expr = 0.0,
                     mean_expr = 0.0,
                     legend_position = "none", 
                     fc_limits = c(-3,3), 
                     plot_labels=T, 
                     node_label_width = 50,
                     group_label_font_size=1){ 
  
  
  g = cell_state_graph@g
  bezier_df = cell_state_graph@layout_info$bezier_df
  
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
  
  # g = left_join(g, deg_table, by = c("name"="cell_group", "term"), relationship = "many-to-many") 
  
  if (is.null(fc_limits)==FALSE) {
    min = fc_limits[1]
    max = fc_limits[2]
    deg_table = deg_table %>%
      mutate(perturb_to_ctrl_shrunken_lfc = ifelse(perturb_to_ctrl_shrunken_lfc > max, max, perturb_to_ctrl_shrunken_lfc)) %>%
      mutate(perturb_to_ctrl_shrunken_lfc = ifelse(perturb_to_ctrl_shrunken_lfc < min, min, perturb_to_ctrl_shrunken_lfc))
  }
  
  
  deg_table = deg_table %>% mutate(
    delta_q_value = pmax(0.0001, perturb_to_ctrl_p_value),
    q_value_sig_code = platt:::calc_sig_ind(perturb_to_ctrl_p_value, html=FALSE))
  
  g = left_join(g, deg_table, by = c("name"="cell_group"), relationship = "many-to-many") 
  
  
  if (is.null(perturb_table) == FALSE) {
    # make non sig dacts change
    perturb_table = perturb_table %>% mutate(delta_log_abund = ifelse(delta_q_value < 0.05, delta_log_abund, 0))
    
    g = left_join(g, perturb_table, by = c("name"="cell_group"), relationship = "many-to-many")
    g = g %>% mutate(size = exp(delta_log_abund)*node_size)
    
  } else {
    g = g %>% mutate(size = node_size)
  }
  
  
  p <- ggplot(aes(x,y), data = g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       colour = con_colour, data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
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
    plot_labels=FALSE 
  }
  
  p =  p + 
    ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(data = g %>% filter(!is.na(perturb_to_ctrl_shrunken_lfc)),
                          mapping = aes(x, y, fill = perturb_to_ctrl_shrunken_lfc, size = size), 
                          shape = "circle filled",
                          color = I("black")) + 
    ggnetwork::geom_nodes(data = g %>% filter(is.na(perturb_to_ctrl_shrunken_lfc)),
                          mapping = aes(x, y, size = size), 
                          fill="lightgray", 
                          shape = "circle filled",
                          color = I("black")) +
    ggnetwork::theme_blank() + 
    scale_fill_gradient2(low = "#006600",  mid = "white", high = "#800080", limits = fc_limits) + 
    # ggnetwork::geom_nodetext(data = g,
    #                          aes(x, y,
    #                              label = q_value_sig_code),
    #                          color=I("black")) + 
    hooke_theme_opts() +
    theme(legend.position=legend_position)
  
  
  if (is.null(perturb_table) == FALSE) {
    p = p + scale_size_identity(guide = "legend") 
    ggb <- ggplot_build(p)
    size_breaks <- ggb$plot$scales$get_scales("size")$get_breaks()
    size_breaks_scale = size_breaks/node_size
    
    new_breaks = round(size_breaks_scale, 1)*node_size
    new_breaks_labels = round(size_breaks_scale, 1)
    p = p + scale_size_continuous(breaks = new_breaks,
                                  labels = new_breaks_labels)
    
  } else {
    p = p + guides(size="none") + scale_size_identity() 
  }
  
  
  p = p + guides(fill=guide_colourbar(title="log(FC)"))
  
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
                     mapping = aes(x,y), color="white", alpha=0) #+ facet_wrap(gene_short_name)
  
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
    # ggnetwork::geom_nodes(data = g,
    #                       aes(x, y) ,
    #                       color=I("black"), size=node_size*1.2) +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y,
                              fill = color_nodes_by), 
                          size=node_size,
                          shape = "circle filled",
                          color = I("black"),) + 
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
    # ggnetwork::geom_nodes(data = g,
    #                       aes(x, y) ,
    #                       color=I("black"), size=node_size*1.2) +
    ggnetwork::geom_nodes(data = g,
                          mapping = aes(x, y, fill = color_nodes_by),
                          shape = "circle filled",
                          color = I("black"),
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



#' plots a cell_state_graph by any column name in the provide table
#' @param cell_state_graph
#' @param table with attributes for each node in the graph
plot_by_table = function(cell_state_graph, 
                         node_table, 
                            color_nodes_by = NULL, 
                            label_nodes_by = NULL, 
                            arrow_unit = 7,
                            node_size = 2,
                            con_colour = "darkgrey", 
                            legend_position = "none", 
                            min_edge_size=0.1,
                            max_edge_size=2,
                            edge_weights=NULL, 
                            plot_labels=T,
                           node_label_width = 50,
                           group_label_font_size=1) {
  
  
  group_nodes_by = cell_state_graph@metadata$group_nodes_by
  
  g = cell_state_graph@g
  
  node_table[["cell_group"]] = node_table[[cell_state_graph@ccs@info$cell_group]]
  g = left_join(g, node_table, by = c("name"="cell_group"))
  g[["color_nodes_by_col"]] = g[[color_nodes_by]]
  
  bezier_df = cell_state_graph@layout_info$bezier_df
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
  
  
  p <- ggplot(aes(x,y), data=g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name, linetype=unsupported_edge), 
                       colour = con_colour, 
                       data = bezier_df %>% distinct(), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
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
    plot_labels=FALSE 
  }
  
  p = p + 
      # ggnetwork::geom_nodes(data = g,
      #                       aes(x, y) ,
      #                       color=I("black"), size=node_size*1.2) +
      ggnetwork::geom_nodes(data = g,
                            aes(x, y,
                                fill = color_nodes_by_col),
                            shape = "circle filled",
                            color = I("black"),
                            size = node_size)  +
      labs(fill = color_nodes_by)
    
  
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black"), 
                                     box.padding = 0.5) 
  }
  
  
  p = p + 
    scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() + 
    theme(legend.position=legend_position)
  
  if (!is.numeric(g[["color_nodes_by_col"]])) {
    p = p + scale_color_manual(values = hooke:::get_colors(length(unique(g$name)), type="vibrant"))
  }
  
  # plot an invisible point to help with nodes not being cut off
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  
  return(p)
}


#' plots a cell_state_graph by any column name in the provide table
#' @param cell_state_graph
#' @param table with attributes for each node in the graph
plot_by_support = function(cell_state_graph, 
                           edge_table, 
                           color_edges_by, 
                           label_nodes_by = NULL, 
                           arrow_unit = 7,
                           node_size = 2,
                           con_colour = "lightgrey", 
                           legend_position = "none", 
                           min_edge_size=0.1,
                           max_edge_size=2,
                           edge_weights=NULL, 
                           plot_labels=T,
                           node_label_width = 50,
                           group_label_font_size=1) {
  
  
  group_nodes_by = cell_state_graph@metadata$group_nodes_by
  
  g = cell_state_graph@g
  
  bezier_df = cell_state_graph@layout_info$bezier_df
  grouping_df = cell_state_graph@layout_info$grouping_df
  
  bezier_df = left_join(bezier_df, edge_table, by = c("from","to"))
  bezier_df[["color_edges_by_col"]] = bezier_df[[color_edges_by]]
  y_plot_range = max(g$y)
  group_label_position_df = g %>%
    select(x, y, group_nodes_by) %>% distinct() %>%
    group_by(group_nodes_by) %>% summarize(x = mean(x), y = max(y) + y_plot_range * 0.02)
  
  
  p <- ggplot(aes(x,y), data=g) + 
    ggplot2::geom_path(aes(x, y, group = edge_name, 
                           color = color_edges_by_col), 
                       data = bezier_df %>% distinct() %>% filter(!is.na(color_edges_by_col)), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre') +
    ggplot2::geom_path(aes(x, y, group = edge_name), 
                       color = con_colour,
                       data = bezier_df %>% distinct() %>% filter(is.na(color_edges_by_col)), 
                       arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                       linejoin='mitre')+
    labs(color = color_edges_by)
    
  
  
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
    plot_labels=FALSE 
  }
  
  p = p + 
    ggnetwork::geom_nodes(data = g,
                          aes(x, y),
                          fill = "black",
                          shape = "circle filled",
                          color = I("black"),
                          size = node_size)  
  
  
  
  if (plot_labels) {
    p = p + ggrepel::geom_text_repel(data= g %>% select(x, y, name) %>% distinct(), 
                                     aes(x, y, label=name),
                                     color=I("black"), 
                                     box.padding = 0.5) 
  }
  
  
  p = p + 
    scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() + 
    theme(legend.position=legend_position)

  
  # plot an invisible point to help with nodes not being cut off
  x_range = range(g$x) + c(-node_size*1.2, node_size*1.2)
  y_range = range(g$y) + c(-node_size*1.2, node_size*1.2)
  point_df = expand.grid("x" = x_range, "y" = y_range)
  p = p + geom_point(point_df,
                     mapping = aes(x,y), color="white", alpha=0)
  
  
  return(p)
}
