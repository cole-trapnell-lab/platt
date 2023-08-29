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
      count(cell_group, !!sym(color_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(color_by_metadata) = c("cell_group", "color_nodes_by")
    node_metadata = left_join(node_metadata, color_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(group_nodes_by) == FALSE){
    group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(group_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
    node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(label_nodes_by) == FALSE){
    label_by_metadata = cell_group_metadata[,c("cell_group", label_nodes_by), drop=F]
    colnames(label_by_metadata) = c("cell_group", "label_nodes_by")
    label_by_metadata = label_by_metadata %>%
      as.data.frame %>%
      count(cell_group, label_nodes_by) %>%
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

  subgraph_df = node_metadata %>%
    select(group_nodes_by, id) %>%
    group_by(group_nodes_by) %>%
    tidyr::nest(subgraph_ids = id) %>%
    summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
                                    .x = subgraph_ids,
                                    G_nel))
  subgraphs = subgraph_df$subgraph
  names(subgraphs) = subgraph_df$group_nodes_by

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

# # try out rewriting the plotting functions


# # plot_state_graph_annotations
# # plot_state_graph_abundance_changes
#     comp_abund_table,
#     contrast = "contrast",
# # plot_state_graph_gene_expression
#     genes,
#     method = "min",
#     gene_expr = NULL,
#     fract_expr=0.0,
#     mean_expr=0.0,
#     scale_to_range = FALSE


# # plot_state_graph_perturb_effects
#     perturbation_table,
#     patterns = c("direct", "indirect", "inferred", "predicted"),
#     num_top_perturbs=3,
#     num_top_genes=3,
# # plot_state_graph_marker_genes
#     marker_table,
#     patterns_for_marker_genes = c("Specifically activated",
#                             "Selectively activated",
#                             "Selectively upregulated",
#                             "Specifically upregulated",
#                             "Specifically maintained",
#                             "Selectively maintained",
#                             "Activated",
#                             "Upregulated",
#                             "Maintained"),
#     num_top_genes=3
# # plot_state_graph_key_genes
# # plot_state_graph_losses



plot_state_graph <- function(state_graph, 
                            color_nodes_by=NULL,
                            label_nodes_by=NULL,
                            group_nodes_by=NULL,
                            label_edges_by=NULL,
                            edge_weights=NULL,
                            fc_limits=c(-3,3),
                            contrast = "contrast",
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
                            group_outline = FALSE, 
                            comp_abund_table = NULL, 
                            genes = NULL, 
                            scale_to_range=FALSE, 
                            gene_expr = NULL) {

    edges = state_graph@graph %>% igraph::as_data_frame()

    edges = edges %>% dplyr::ungroup()
    
    if (color_nodes_by == "gene_expression") {
      color_nodes_by = NULL
      color_cells_by = "gene_expression" 
    }

    node_metadata = collect_psg_node_metadata(state_graph@ccs, color_nodes_by, label_nodes_by, group_nodes_by)

    if (hide_unlinked_nodes){
        node_metadata = node_metadata %>% dplyr::filter(id %in% edges$from | id %in% edges$to)
    }

    if (is.null(edge_weights)){
        edges = edges %>% select(from, to)
        edges$weight = 1
    }else{
        edges = edges %>% select(from, to, weight=!!sym(edge_weights))
    }

    if (is.null(label_edges_by)){
        edge_info = edges %>% select(from, to)
    } else {
        if (is(state_graph, "igraph")){
            edge_info = state_graph %>% igraph::as_data_frame() %>% select(from, to, label=!!sym(label_edges_by))
    } else {
        edge_info = state_graph %>% select(from, to, label=!!sym(label_edges_by))
    }

        edges = edges %>% left_join(edge_info)
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
    }else{
        edge_labels=NULL
    }   


    # if () {
    # 
    #     node_metadata <- perturbation_nodes(node_metadata, perturbation_table)
    # }



    layout_info = layout_state_graph(G, node_metadata, edge_labels =edge_labels, weighted=FALSE)
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
    bezier_df$assembly_group = vapply(strsplit(bezier_df$from, "-", fixed=T), "[", "", 1)
    g$assembly_group = vapply(strsplit(g$name, "-", fixed=T), "[", "", 1)
    

    if (is.null(comp_abund_table) == FALSE) {
        g <- add_abundance_changes(g, comp_abund_table, facet_group = contrast, fc_limits = fc_limits)
        color_nodes_by = "delta_log_abund"
        group_outline = TRUE # for some reason cant get the geom_mark_rect to facet
    }
    
    if (is.null(genes) == FALSE) {
      g <- add_gene_expression(g, state_graph@ccs, genes=genes, gene_expr=gene_expr, scale_to_range=scale_to_range)
      color_nodes_by = "mean_expression"
      expression_legend_label = "mean expression"
    }
    
    
    # if (color_cells_by == "") {
    # 
    #     if (scale_to_range) {
    #         color_nodes_by = "mean_expression"
    #         expression_legend_label = "mean expression"
    #     }
    #     
    #     expression_legend_label = "% Max"
    #     color_nodes_by = "value"
    # 
    #     g <- add_gene_expression(ccs, 
    #                             genes, 
    #                             method = "min",
    #                             gene_expr = NULL,
    #                             fract_expr=0.0,
    #                             mean_expr=0.0,
    #                             scale_to_range = FALSE)
    # }
    
    
    
    p <- ggplot() +
        ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), 
                           colour=con_colour, data=bezier_df %>% distinct(), 
                           arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
                           linejoin='mitre')

    if (is.null(edge_labels) == FALSE) {
        label_df = layout_info$label_df
        p = p + geom_text(data = label_df, aes(x,y, label = label), size=edge_label_font_size)
    }
    
    p = p + ggforce::geom_mark_rect(data=g, 
                                    aes(x, y,
                                        col = switch(group_outline, group_nodes_by,NULL), 
                                        fill = switch(!group_outline, group_nodes_by,NULL), 
                                        label = switch(label_groups, group_nodes_by, NULL), 
                                        filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                    size=ifelse(group_outline, 0.5, 0),
                                    label.fontsize=group_label_font_size,
                                    con.linetype=label_conn_linetype,
                                    con.colour=con_colour,
                                    expand = unit(2, "mm"),
                                    label.buffer=unit(1, "mm"),
                                    radius = unit(1.5, "mm"),
                                    label.margin = margin(1, 1, 1, 1, "mm"),
                                    label.fontface="plain")
    
    # if it is numeric, use a color scale
    color_grad = switch(is.numeric(g[[color_nodes_by]]), 
                        scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3"), NULL)

    # color_nodes_by = ifelse(is.null(color_nodes_by))
    
    # either use geom_nodel label 
    if (is.null(label_nodes_by) == FALSE) {
        p = p + ggnewscale::new_scale_fill() +
            ggnetwork::geom_nodelabel(data = g,
                                        aes(x, y,
                                            fill = color_nodes_by, 
                                            label = switch(is.null(label_nodes_by) == FALSE, label_nodes_by, NULL)),
                                        size = node_size) +
            labs(fill = color_nodes_by) + theme(legend.position = "none") + 
            color_grad 
    # or plot dots
    } else {
      
      if (is.null(comp_abund_table) == FALSE){
        
        p = p + guides(fill = "none")
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
                                   color=I("black")) + 
          scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")
        p = p + facet_wrap(~contrast) 
        
      } else if (is.null(genes)== FALSE) {
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
          viridis::scale_color_viridis(option = "viridis") + facet_wrap(~gene_short_name)
      }else {
        p = p + 
          ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(data = g,
                                aes(x, y,
                                    color = color_nodes_by),
                                size = node_size) +
          labs(fill = color_nodes_by) + theme(legend.position = "none") + 
          color_grad
      }
      

    }
    
    if (is.null(edge_labels) == FALSE) {
        label_df = layout_info$bezier_df %>%
        group_by(edge_name) %>%
        summarize(x = mean(x), y=mean(y))
        label_df$label = edge_labels[label_df$edge_name]
        p = p + ggrepel::geom_text_repel(data = label_df,
                                    mapping = aes(x, y, label=label),
                                    size=edge_label_font_size)
    }

    p = p + 
        # scale_size(range=c(1, node_size))+ 
      scale_size_identity() +
        monocle3:::monocle_theme_opts() +
        ggnetwork::theme_blank() +
        theme(legend.position=legend_position)

    return(p)

}


# plot_state_graph <- function(state_graph) {
#   
#   if (type == "abundance") {
#     
#   } else if (type == "expression") {
#     
#   } else if (type == "perturb effects") {
#     
#   }
# }


# plot_state_graph(state_graph, type = "abundance_changes")





add_abundance_changes <- function(g, 
                                  comp_abund_table, 
                                  facet_group = "", 
                                  fc_limits = c(-3,3)) {

    if (is.null(facet_group)) {
        comp_abund_table["contrast"] = ""
    } else {
        comp_abund_table[["contrast"]] = comp_abund_table[[facet_group]]
        comp_abund_table$contrast = as.factor(comp_abund_table$contrast)
    }

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
    # color_nodes_by = "delta_log_abund"
    return(g)
}


add_gene_expression <- function(g, 
                                ccs, 
                                genes, 
                                method = "min",
                                gene_expr = NULL,
                                fract_expr=0.0,
                                mean_expr=0.0,
                                scale_to_range = FALSE) {

    gene_ids = rowData(ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% genes) %>% rownames()
    
    if (is.null(gene_expr)) {
        gene_expr = hooke:::aggregated_expr_data(ccs@cds[gene_ids,], group_cells_by = ccs@info$cell_group)
    }

    sub_gene_expr = gene_expr %>%
        filter(gene_short_name %in% genes)

    method = get(method)        
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
    
    return(g)
}

# earliest_loss_nodes <- function(node_metadata, 
#                                 perturbation_ccm,
#                                 start_time=start_time,
#                                 stop_time=stop_time,
#                                 interval_step = interval_step,
#                                 interval_col=interval_col,
#                                 log_abund_detection_thresh=log_abund_detection_thresh,
#                                 q_val = q_val,
#                                 control_ccm=control_ccm,
#                                 control_start_time=control_start_time,
#                                 control_stop_time=control_stop_time,
#                                 ...) {
#     earliest_loss_tbl = hooke:::estimate_loss_timing(perturbation_ccm,
#                                                    start_time=start_time,
#                                                    stop_time=stop_time,
#                                                    interval_step = interval_step,
#                                                    interval_col=interval_col,
#                                                    log_abund_detection_thresh=log_abund_detection_thresh,
#                                                    q_val = q_val,
#                                                    control_ccm=control_ccm,
#                                                    control_start_time=control_start_time,
#                                                    control_stop_time=control_stop_time,
#                                                    ...)

#   #earliest_loss_tbl = earliest_loss_tbl %>% mutate(fill_alpha = ifelse(peak_time_in_ctrl_within_perturb_time_range, 1.0, 0.3))
#   #print (earliest_loss_tbl)
#   node_metadata = node_metadata %>% left_join(earliest_loss_tbl, by=c("id" ="cell_group"))
#   return(node_metadata)
# }



# plot_everything_together <- function(g, bezier_df) {
#     p <- ggplot(aes(x,y), data=g)
    
#     p <- ggplot() +
#         ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), 
#                            colour=con_colour, data=bezier_df %>% distinct(), 
#                            arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), 
#                            linejoin='mitre')

# }

# # nodes by annotation

# plot_annotation_nodes <- function(node_metadata) {

# }



perturbation_nodes <- function(node_metadata, perturbation_table) {

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


    node_metadata = node_metadata %>% left_join(perturbation_table %>% select(id, perturb_effect_label), by=c("id"="id"))
    if (is.null(label_nodes_by)){
        label_nodes_by = "perturb_effect_label"
        node_metadata = node_metadata %>% mutate(label_nodes_by=perturb_effect_label)
    }else{
        node_metadata = node_metadata %>% mutate(label_nodes_by=ifelse(is.na(perturb_effect_label), label_nodes_by, paste(label_nodes_by, perturb_effect_label, sep="<br>")))
    }

    return(node_metadata)
}

