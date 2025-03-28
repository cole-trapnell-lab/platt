
#' Plot State Graph Annotations
#'
#' This function generates a state graph visualization with customizable node and edge attributes.
#'
#' @param ccs A cell count set object containing metadata for nodes.
#' @param state_graph An igraph object or data frame representing the state graph.
#' @param color_nodes_by Column name in `ccs` to color nodes by. Default is `NULL`.
#' @param label_nodes_by Column name in `ccs` to label nodes by. Default is `NULL`.
#' @param group_nodes_by Column name in `ccs` to group nodes by. Default is `NULL`.
#' @param label_edges_by Column name in `state_graph` to label edges by. Default is `NULL`.
#' @param edge_weights Column name in `state_graph` for edge weights. Default is `NULL`.
#' @param arrow.gap Numeric value for the gap between arrows and nodes. Default is `0.03`.
#' @param arrow_unit Numeric value for arrow size in points. Default is `2`.
#' @param bar_unit Numeric value for bar size. Default is `0.075`.
#' @param node_size Numeric value for node size. Default is `2`.
#' @param min_edge_size Minimum edge thickness. Default is `0.1`.
#' @param max_edge_size Maximum edge thickness. Default is `2`.
#' @param unlabeled_groups Character vector of group names to exclude from labeling. Default is `c("Unknown")`.
#' @param label_groups Logical indicating whether to label groups. Default is `TRUE`.
#' @param hide_unlinked_nodes Logical indicating whether to hide nodes not linked by edges. Default is `TRUE`.
#' @param group_label_font_size Font size for group labels. Default is `6`.
#' @param edge_label_font_size Font size for edge labels. Default is `2`.
#' @param label_conn_linetype Line type for group label connectors. Default is `"dotted"`.
#' @param legend_position Position of the legend in the plot. Default is `"none"`.
#' @param con_colour Color for connectors and edges. Default is `"darkgrey"`.
#' @param group_outline Logical indicating whether to draw an outline around groups. Default is `FALSE`.
#'
#' @return A ggplot object representing the state graph visualization.
#'
#' @details
#' The function allows for flexible customization of node and edge attributes, including coloring, labeling, and grouping. 
#' It supports both igraph objects and data frames as input for the state graph. Nodes and edges can be filtered, styled, 
#' and annotated based on the provided metadata.
#'
#' @examples
#' # Example usage:
#' plot_state_graph_annotations(
#'   ccs = my_ccs,
#'   state_graph = my_graph,
#'   color_nodes_by = "group",
#'   label_nodes_by = "name",
#'   edge_weights = "weight"
#' )
#'
#' @import ggplot2
#' @importFrom dplyr filter select mutate left_join ungroup group_by summarize
#' @importFrom igraph as_data_frame graph_from_data_frame E
#' @importFrom ggnetwork ggnetwork geom_nodes geom_nodelabel theme_blank
#' @importFrom ggforce geom_mark_rect
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_fill
#' @importFrom stringr str_c
#' @importFrom tidyr replace_na
#' @importFrom grid unit
#' @export 
plot_state_graph_annotations <- function(ccs,
                                         state_graph,
                                         color_nodes_by = NULL,
                                         label_nodes_by = NULL,
                                         group_nodes_by = NULL,
                                         label_edges_by = NULL,
                                         edge_weights = NULL,
                                         arrow.gap = 0.03,
                                         arrow_unit = 2,
                                         bar_unit = .075,
                                         node_size = 2,
                                         min_edge_size = 0.1,
                                         max_edge_size = 2,
                                         unlabeled_groups = c("Unknown"),
                                         label_groups = TRUE,
                                         hide_unlinked_nodes = TRUE,
                                         group_label_font_size = 6,
                                         edge_label_font_size = 2,
                                         label_conn_linetype = "dotted",
                                         legend_position = "none",
                                         con_colour = "darkgrey",
                                         group_outline = FALSE) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
    edges$weight <- 1
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }


  if (length(setdiff(edges$to, node_metadata$id)) > 0 ||
    length(setdiff(edges$from, node_metadata$id)) > 0) {
    message("Warning: graph edges refers to nodes not in cell_count_set ccs. Dropping edges...")
  }

  edges <- edges %>% filter(from %in% node_metadata$id & to %in% node_metadata$id)

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, NULL, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
      unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
      edge_thickness = replace_na(edge_thickness, min_edge_size)
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge <- FALSE
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)

  # add assembly group
  bezier_df$assembly_group <- vapply(strsplit(bezier_df$from, "-", fixed = T), "[", "", 1)
  g$assembly_group <- vapply(strsplit(g$name, "-", fixed = T), "[", "", 1)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group = edge_name, size = edge_thickness, linetype = unsupported_edge), colour = con_colour, data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")

  # draw activator edges
  # ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE) {
    if (group_outline) {
      # if (label_groups) {
      #   label_groups_by = group_nodes_by
      # } else {
      #   label_groups_by = NULL
      # }
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          col = group_nodes_by,
          # label = label_groups_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0.5,
        expand = unit(2, "mm"),
        label.buffer = unit(1, "mm"),
        radius = unit(1.5, "mm"),
        label.margin = margin(1, 1, 1, 1, "mm"),
        label.fontsize = group_label_font_size,
        label.fontface = "plain",
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        data = g
      )
    } else {
      if (label_groups) {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            label = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          expand = unit(2, "mm"),
          label.buffer = unit(1, "mm"),
          radius = unit(1.5, "mm"),
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = group_label_font_size,
          label.fontface = "plain",
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      } else {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            col = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          expand = unit(2, "mm"),
          label.buffer = unit(1, "mm"),
          radius = unit(1.5, "mm"),
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = group_label_font_size,
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      }
    }
  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE) {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodelabel(
            data = g,
            aes(x, y,
              fill = !!sym(color_nodes_by),
              label = label_nodes_by
            ),
            size = node_size
          ) +
          labs(fill = color_nodes_by)
        p <- p + scale_fill_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + ggnetwork::geom_nodelabel(
          data = g,
          aes(x, y,
            fill = color_nodes_by,
            label = label_nodes_by
          ),
          size = node_size
        ) +
          labs(fill = color_nodes_by)
      }
    } else {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(
            data = g,
            aes(x, y,
              color = !!sym(color_nodes_by)
            ),
            size = node_size
          ) +
          labs(color = color_nodes_by)
        p <- p + scale_color_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + ggnetwork::geom_nodes(
          data = g,
          aes(x, y,
            color = color_nodes_by
          ),
          size = node_size
        ) +
          labs(color = color_nodes_by)
      }
    }
  } else {
    if (is.null(label_nodes_by) == FALSE) {
      p <- p + ggnetwork::geom_nodelabel(
        data = g,
        aes(x, y,
          xend = xend, yend = yend,
          label = label_nodes_by
        ),
        size = node_size
      )
    } else {
      p <- p + ggnetwork::geom_nodes(
        data = g,
        aes(x, y, xend = xend, yend = yend),
        size = node_size
      )
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df <- layout_info$bezier_df %>%
      group_by(edge_name) %>%
      summarize(x = mean(x), y = mean(y))
    label_df$label <- edge_labels[label_df$edge_name]
    # label_df = layout_info$label_df
    # p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    # p = p + geom_text(data = label_df,
    #                  aes(x,y, label = label),
    #                  size=edge_label_font_size)
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      mapping = aes(x, y, label = label),
      size = edge_label_font_size
    )
  }

  p <- p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position = legend_position)
  return(p)
}

#' This function plots when each cell state is lost following a perturbation
#'
#' FIXME: this is really a prototype and needs lots of cleanup and refactoring.
#' It could plot the earliest time a loss is detected, the latest, etc. But
#' right now it only plots the loss at the time of peak abundance in the control
#' see estimate_loss_timing() for more details.
#'
plot_state_graph_losses <- function(perturbation_ccm,
                                    state_graph,
                                    start_time,
                                    stop_time,
                                    interval_step,
                                    interval_col,
                                    log_abund_detection_thresh,
                                    q_val,
                                    loss_time = c("largest_loss", "largest_loss_time", "earliest_time", "latest_time", "peak_loss_time", "delta_log_abund_at_peak"),
                                    label_nodes_by = NULL,
                                    group_nodes_by = NULL,
                                    label_edges_by = NULL,
                                    edge_weights = NULL,
                                    arrow.gap = 0.03,
                                    arrow_unit = 2,
                                    bar_unit = .075,
                                    node_size = 2,
                                    min_edge_size = 0.1,
                                    max_edge_size = 2,
                                    unlabeled_groups = c("Unknown"),
                                    label_subset = NULL,
                                    label_groups = TRUE,
                                    hide_unlinked_nodes = TRUE,
                                    group_label_font_size = 6,
                                    edge_label_font_size = 2,
                                    label_conn_linetype = "dotted",
                                    legend_position = "none",
                                    con_colour = "darkgrey",
                                    group_outline = FALSE,
                                    control_ccm = perturbation_ccm,
                                    control_start_time = NULL,
                                    control_stop_time = NULL,
                                    newdata = tibble()) {
  loss_time <- match.arg(loss_time)

  if (is.null(control_start_time)) {
    control_start_time <- start_time
  }
  if (is.null(control_stop_time)) {
    control_stop_time <- stop_time
  }

  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(perturbation_ccm@ccs, color_nodes_by = NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }

  earliest_loss_tbl <- estimate_loss_timing(perturbation_ccm,
    start_time = start_time,
    stop_time = stop_time,
    interval_step = interval_step,
    interval_col = interval_col,
    log_abund_detection_thresh = log_abund_detection_thresh,
    q_val = q_val,
    control_ccm = control_ccm,
    control_start_time = control_start_time,
    control_stop_time = control_stop_time,
    newdata = newdata
  )

  # earliest_loss_tbl = earliest_loss_tbl %>% mutate(fill_alpha = ifelse(peak_time_in_ctrl_within_perturb_time_range, 1.0, 0.3))
  # print (earliest_loss_tbl)
  node_metadata <- node_metadata %>% left_join(earliest_loss_tbl, by = c("id" = "cell_group"))

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, NULL, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)
  # add assembly group
  bezier_df$assembly_group <- vapply(strsplit(bezier_df$from, "-", fixed = T), "[", "", 1)
  g$assembly_group <- vapply(strsplit(g$name, "-", fixed = T), "[", "", 1)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group = edge_name, size = edge_thickness), data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")

  if (is.null(label_subset)) {
    label_subset <- unique(node_metadata$group_nodes_by)
  }

  label_subset <- label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE) {
    if (group_outline) {
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          col = group_nodes_by,
          label = group_nodes_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0.5,
        expand = unit(2, "mm"),
        label.buffer = unit(1, "mm"),
        radius = unit(1.5, "mm"),
        label.margin = margin(1, 1, 1, 1, "mm"),
        label.fontsize = group_label_font_size,
        label.fontface = "plain",
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        data = g
      )
    } else {
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          fill = group_nodes_by,
          label = group_nodes_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0,
        expand = unit(2, "mm"),
        label.buffer = unit(1, "mm"),
        radius = unit(1.5, "mm"),
        label.margin = margin(1, 1, 1, 1, "mm"),
        label.fontsize = group_label_font_size,
        label.fontface = "plain",
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        data = g
      )
    }
  }


  # if numerical

  min_delta_log_abund <- -3
  max_delta_log_abund <- 3
  g <- g %>% mutate(fill_val = !!sym(loss_time))

  if (loss_time %in% c("largest_loss", "delta_log_abund_at_peak")) {
    g <- g %>% mutate(
      fill_val = ifelse(fill_val < min_delta_log_abund, min_delta_log_abund, fill_val),
      fill_val = ifelse(fill_val > max_delta_log_abund, max_delta_log_abund, fill_val),
      fill_val = ifelse(peak_time_in_ctrl_within_perturb_time_range, fill_val, NA),
      fill_border_color = ifelse(peak_time_in_ctrl_within_perturb_time_range == FALSE,
        "out-of-range",
        ifelse(is_lost_at_peak,
          "significant",
          "notsignificant"
        )
      )
    )
  }

  p <- p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodelabel(
      data = g,
      aes(x, y,
        fill = fill_val,
        color = fill_border_color,
        label = label_nodes_by
      ),
      size = node_size
    )
  p <- p + scale_color_manual(values = c("out-of-range" = "black", "significant" = "black", "notsignificant" = "lightgrey"))

  if (loss_time %in% c("largest_loss", "delta_log_abund_at_peak")) {
    p <- p + scale_fill_gradient2(
      low = "royalblue3", mid = "white", high = "orangered3",
      limits = c(
        min_delta_log_abund - abs(min_delta_log_abund) * 0.05,
        max_delta_log_abund + abs(max_delta_log_abund) * 0.05
      )
    )
  } else {
    p <- p + scale_fill_stepsn(n.breaks = 5, colours = terrain.colors(5)) #+ scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
  }

  p <- p + guides(color = "none")

  # p = p + guides(fill = guide_legend(position=legend_position))

  if (is.null(edge_labels) == FALSE) {
    label_df <- layout_info$bezier_df %>%
      group_by(edge_name) %>%
      summarize(x = mean(x), y = mean(y))
    label_df$label <- edge_labels[label_df$edge_name]
    # label_df = layout_info$label_df
    # p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    # p = p + geom_text(data = label_df,
    #                  aes(x,y, label = label),
    #                  size=edge_label_font_size)
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      mapping = aes(x, y, label = label),
      size = edge_label_font_size
    )
  }

  p <- p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() #+
  # theme(legend.position=legend_position)
  return(p)
}


#' @noRd
num_extract <- function(string, as_char = TRUE) {
  if (as_char) {
    stringr::str_trim(
      format(stringr::str_extract(string, "\\-*\\d+\\.*\\d*"),
        scientific = FALSE,
        trim = TRUE
      )
    )
  } else {
    as.numeric(
      stringr::str_trim(
        format(stringr::str_extract(string, "\\-*\\d+\\.*\\d*"),
          scientific = FALSE,
          trim = TRUE
        )
      )
    )
  }
}

#' @noRd
calc_sig_ind <- function(p_value, html = TRUE) {
  # p_value <- suppressWarnings(
  #  num_extract(p_value, as_char = FALSE)
  # )

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
calc_sig_rank <- function(p_value) {
  # p_value <- suppressWarnings(
  #  num_extract(p_value, as_char = FALSE)
  # )
  dplyr::case_when(
    p_value <= 0 ~ 2,
    p_value <= 0.001 ~ 2,
    p_value <= 0.01 ~ 1,
    p_value <= 0.05 ~ 1,
    p_value <= 0.1 ~ 0.5,
    p_value <= 1 ~ 0.5,
    TRUE ~ 1
  )
}

#' Plot State Graph Abundance Changes
#'
#' This function plots the changes in abundance across different states in a state graph.
#'
#' @param ccs A data frame containing cell cluster information.
#' @param state_graph An igraph object or data frame representing the state graph.
#' @param comp_abund_table A data frame containing the comparison abundance table.
#' @param contrast A string specifying the contrast column in the comparison abundance table. Default is "contrast".
#' @param label_nodes_by A string specifying the column to label nodes by. Default is NULL.
#' @param group_nodes_by A string specifying the column to group nodes by. Default is NULL.
#' @param label_edges_by A string specifying the column to label edges by. Default is NULL.
#' @param edge_weights A string specifying the column to use for edge weights. Default is NULL.
#' @param size_nodes_by A string specifying the column to size nodes by. Default is NULL.
#' @param num_degs A data frame containing the number of differentially expressed genes. Default is NULL.
#' @param fc_limits A numeric vector specifying the fold change limits. Default is c(-3, 3).
#' @param arrow.gap A numeric value specifying the gap for arrows. Default is 0.03.
#' @param arrow_unit A numeric value specifying the unit for arrows. Default is 2.
#' @param bar_unit A numeric value specifying the unit for bars. Default is 0.075.
#' @param node_size A numeric value specifying the size of nodes. Default is 6.
#' @param min_edge_size A numeric value specifying the minimum edge size. Default is 0.1.
#' @param max_edge_size A numeric value specifying the maximum edge size. Default is 2.
#' @param unlabeled_groups A character vector specifying the groups to be unlabeled. Default is c("Unknown").
#' @param label_subset A character vector specifying the subset of labels. Default is NULL.
#' @param label_groups A logical value specifying whether to label groups. Default is TRUE.
#' @param hide_unlinked_nodes A logical value specifying whether to hide unlinked nodes. Default is TRUE.
#' @param group_label_font_size A numeric value specifying the font size for group labels. Default is 6.
#' @param edge_label_font_size A numeric value specifying the font size for edge labels. Default is 2.
#' @param label_conn_linetype A string specifying the line type for label connections. Default is "dotted".
#' @param legend_position A string specifying the position of the legend. Default is "none".
#' @param con_colour A string specifying the color for connections. Default is "darkgrey".
#' @param group_outline A logical value specifying whether to outline groups. Default is FALSE.
#' @param flip_x A logical value specifying whether to flip the x-axis. Default is FALSE.
#' @return A ggplot object representing the state graph with abundance changes.
#' @import dplyr
#' @import ggplot2
#' @import ggnetwork
#' @import ggforce
#' @import ggnewscale
#' @import igraph
#' @import stringr
plot_state_graph_abundance_changes <- function(ccs,
                                               state_graph,
                                               comp_abund_table,
                                               contrast = "contrast",
                                               label_nodes_by = NULL,
                                               group_nodes_by = NULL,
                                               label_edges_by = NULL,
                                               edge_weights = NULL,
                                               size_nodes_by = NULL,
                                               num_degs = NULL,
                                               fc_limits = c(-3, 3),
                                               arrow.gap = 0.03,
                                               arrow_unit = 2,
                                               bar_unit = .075,
                                               node_size = 6,
                                               min_edge_size = 0.1,
                                               max_edge_size = 2,
                                               unlabeled_groups = c("Unknown"),
                                               label_subset = NULL,
                                               label_groups = TRUE,
                                               hide_unlinked_nodes = TRUE,
                                               group_label_font_size = 6,
                                               edge_label_font_size = 2,
                                               label_conn_linetype = "dotted",
                                               legend_position = "none",
                                               con_colour = "darkgrey",
                                               group_outline = FALSE,
                                               flip_x = FALSE) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by = NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, NULL, weighted = FALSE)

  comp_abund_table[["contrast"]] <- comp_abund_table[[contrast]]
  comp_abund_table$contrast <- as.factor(comp_abund_table$contrast)

  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)

  abund_fc_df <- comp_abund_table

  if (is.null(fc_limits)) {
    fc_limits <- range(abund_fc_df$delta_log_abund)
  } else {
    min <- fc_limits[1]
    max <- fc_limits[2]
    abund_fc_df <- abund_fc_df %>%
      mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
      mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  }

  abund_fc_df <- abund_fc_df %>%
    mutate(
      delta_q_value = pmax(0.0001, delta_q_value),
      q_value_sig_code = calc_sig_ind(delta_q_value, html = FALSE)
    )

  g <- left_join(g, abund_fc_df, by = c("name" = "cell_group"), relationship = "many-to-many")

  # g$gene_short_name = factor(g$gene_short_name, levels=genes)
  color_nodes_by <- "delta_log_abund"
  # group_outline = TRUE

  if (flip_x) {
    g$x <- -1 * g$x
    bezier_df$x <- -1 * bezier_df$x
  }

  p <- ggplot(aes(x, y), data = g)
  p <- p +
    ggplot2::geom_path(aes(x, y, group = edge_name), colour = con_colour, data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")

  if (is.null(label_subset)) {
    label_subset <- unique(node_metadata$group_nodes_by)
  }

  label_subset <- label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE) {
    if (group_outline) {
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          fill = contrast,
          color = group_nodes_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0.5,
        expand = unit(2, "mm"),
        label.buffer = unit(1, "mm"),
        radius = unit(1.5, "mm"),
        label.margin = margin(1, 1, 1, 1, "mm"),
        label.fontsize = group_label_font_size,
        label.fontface = "plain",
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        show.legend = F
      )
    } else {
      values <- rep("white", length(label_subset))
      names(values) <- label_subset
      if (label_groups) {
        p <- p +
          ggforce::geom_mark_rect(
            aes(x, y,
              fill = contrast,
              color = group_nodes_by,
              label = group_nodes_by,
              filter = group_nodes_by %in% label_subset
            ),
            size = 0.5,
            expand = unit(2, "mm"),
            label.buffer = unit(1, "mm"),
            radius = unit(1.5, "mm"),
            label.margin = margin(1, 1, 1, 1, "mm"),
            label.fontsize = group_label_font_size,
            label.fontface = "plain",
            con.linetype = label_conn_linetype,
            con.colour = con_colour,
            show.legend = F
          ) +
          scale_color_manual(values = c(values))
      } else {
        p <- p +
          ggforce::geom_mark_rect(
            aes(x, y,
              fill = contrast,
              color = group_nodes_by,
              filter = group_nodes_by %in% label_subset
            ),
            size = 0.5,
            expand = unit(2, "mm"),
            label.buffer = unit(1, "mm"),
            radius = unit(1.5, "mm"),
            label.margin = margin(1, 1, 1, 1, "mm"),
            label.fontsize = group_label_font_size,
            label.fontface = "plain",
            con.linetype = label_conn_linetype,
            con.colour = con_colour,
            show.legend = F
          ) +
          scale_color_manual(values = c(values))
      }
    }

    p <- p + scale_fill_manual(values = rep("white", length(unique(g$contrast))))
    # p = p + scale_color_manual(values = rainbow_colors)
  }
  p <- p + guides(fill = "none")

  if (!is.null(size_nodes_by)) {
    g <- left_join(g, num_degs, by = c("name" = "cell_group", "perturb_name"))
    g[["size_nodes_by"]] <- as.numeric(g[[size_nodes_by]])
    p <- p + ggnewscale::new_scale_color() +
      ggnetwork::geom_nodes(
        data = g,
        aes(x, y,
          size = size_nodes_by * 1.2
        ),
        color = I("black")
      ) +
      ggnetwork::geom_nodes(
        data = g,
        aes(x, y,
          size = size_nodes_by,
          color = delta_log_abund
        )
      ) +
      ggnetwork::geom_nodetext(
        data = g,
        aes(x, y,
          label = q_value_sig_code
        ),
        color = I("black")
      )
  } else {
    p <- p + ggnewscale::new_scale_color() +
      ggnetwork::geom_nodes(
        data = g,
        aes(x, y,
          size = -log10(delta_q_value) * 1.2
        ),
        color = I("black")
      ) +
      ggnetwork::geom_nodes(
        data = g,
        aes(x, y,
          size = -log10(delta_q_value),
          color = delta_log_abund
        )
      ) +
      ggnetwork::geom_nodetext(
        data = g,
        aes(x, y,
          label = q_value_sig_code
        ),
        color = I("black")
      )
  }


  # ggnetwork::geom_nodetext_repel(data = g,
  #                                aes(x, y,
  #                                    label = q_value_sig_code),
  #                                color=I("black"))
  # labs(fill = color_nodes_by)

  p <- p + scale_color_gradient2(low = "royalblue3", mid = "white", high = "orangered3")

  if (is.null(edge_labels) == FALSE) {
    p <- p + ggnetwork::geom_nodetext(
      data = label_df,
      aes(x, y, label = label), size = 3
    )
  }

  # p = p +
  #  ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))

  p <- p + facet_wrap(~contrast)

  p <- p + scale_size(range = c(1, node_size)) +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position = legend_position) + guides(fill = "none", size = "none")
  return(p)
}

#' Plot State Graph Gene Expression
#'
#' This function plots a state graph with gene expression data.
#'
#' @param ccs Cell cluster data.
#' @param state_graph State graph data, either as an igraph object or a data frame.
#' @param genes A vector of gene names to plot.
#' @param method Method to aggregate gene expression data, default is "min".
#' @param gene_expr Gene expression data, default is NULL.
#' @param fract_expr Minimum fraction of cells expressing the gene, default is 0.0.
#' @param mean_expr Minimum mean expression level, default is 0.0.
#' @param scale_to_range Logical, whether to scale expression data to range, default is FALSE.
#' @param color_nodes_by Column name to color nodes by, default is NULL.
#' @param label_nodes_by Column name to label nodes by, default is NULL.
#' @param group_nodes_by Column name to group nodes by, default is NULL.
#' @param label_edges_by Column name to label edges by, default is NULL.
#' @param edge_weights Column name for edge weights, default is NULL.
#' @param arrow.gap Gap size for arrows, default is 0.03.
#' @param arrow_unit Unit size for arrows, default is 2.
#' @param bar_unit Unit size for bars, default is 0.075.
#' @param min_node_size Minimum node size, default is 0.25.
#' @param max_node_size Maximum node size, default is 2.
#' @param min_edge_size Minimum edge size, default is 0.1.
#' @param max_edge_size Maximum edge size, default is 2.
#' @param unlabeled_groups Vector of unlabeled groups, default is c("Unknown").
#' @param label_groups Logical, whether to label groups, default is TRUE.
#' @param hide_unlinked_nodes Logical, whether to hide unlinked nodes, default is TRUE.
#' @param group_label_font_size Font size for group labels, default is 6.
#' @param edge_label_font_size Font size for edge labels, default is 2.
#' @param label_conn_linetype Line type for label connections, default is "dotted".
#' @param legend_position Position of the legend, default is "none".
#' @param con_colour Color for connections, default is "darkgrey".
#' @param group_outline Logical, whether to outline groups, default is FALSE.
#' @return A ggplot object representing the state graph with gene expression data.
#' @import dplyr
#' @import igraph
#' @import ggplot2
#' @import ggnetwork
#' @import ggforce
#' @import ggnewscale
#' @import viridis
#' @importFrom stringr str_c
#' @importFrom scales percent
#' @importFrom tidyr replace_na
plot_state_graph_gene_expression <- function(ccs,
                                             state_graph,
                                             genes,
                                             method = "min",
                                             gene_expr = NULL,
                                             fract_expr = 0.0,
                                             mean_expr = 0.0,
                                             scale_to_range = FALSE,
                                             color_nodes_by = NULL,
                                             label_nodes_by = NULL,
                                             group_nodes_by = NULL,
                                             label_edges_by = NULL,
                                             edge_weights = NULL,
                                             arrow.gap = 0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             min_node_size = 0.25,
                                             max_node_size = 2,
                                             min_edge_size = 0.1,
                                             max_edge_size = 2,
                                             unlabeled_groups = c("Unknown"),
                                             label_groups = TRUE,
                                             hide_unlinked_nodes = TRUE,
                                             group_label_font_size = 6,
                                             edge_label_font_size = 2,
                                             label_conn_linetype = "dotted",
                                             legend_position = "none",
                                             con_colour = "darkgrey",
                                             group_outline = FALSE) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by = NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, NULL, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
      unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
      edge_thickness = replace_na(edge_thickness, min_edge_size)
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge <- FALSE
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)

  gene_ids <- rowData(ccs@cds) %>%
    as.data.frame() %>%
    filter(gene_short_name %in% genes) %>%
    rownames()

  if (is.null(gene_expr)) {
    gene_expr <- hooke:::aggregated_expr_data(ccs@cds[gene_ids, ], group_cells_by = ccs@info$cell_group)
  }

  sub_gene_expr <- gene_expr %>%
    filter(gene_short_name %in% genes)

  method <- get(method)
  # sub_gene_expr = sub_gene_expr %>%
  #   group_by(cell_group) %>%
  #   dplyr::summarise(fraction_expressing = method(fraction_expressing),
  #                    mean_expression = method(mean_expression),
  #                    specificity = method(specificity)) %>%
  #   mutate(gene_short_name =  paste(genes, collapse = "-"))


  # node_metadata = node_metadata %>% left_join(sub_gene_expr, by = c("id" = "cell_group"))

  sub_gene_expr <- sub_gene_expr %>%
    group_by(gene_short_name) %>%
    mutate(
      max_expr = max(mean_expression),
      fraction_max = ifelse(max_expr > 0, mean_expression / max_expr, 0),
      gene_expr = case_when(
        fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
        TRUE ~ FALSE
      )
    )

  color_nodes_by <- "mean_expression"
  expression_legend_label <- "mean expression"

  if (scale_to_range) {
    sub_gene_expr <- sub_gene_expr %>%
      mutate(value = mean_expression) %>%
      group_by(gene_short_name) %>%
      dplyr::mutate(
        max_val_for_feature = max(value),
        min_val_for_feature = min(value)
      ) %>%
      dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
    expression_legend_label <- "% Max"
    color_nodes_by <- "value"
  }

  g <- left_join(g, sub_gene_expr, by = c("name" = "cell_group"))

  g$gene_short_name <- factor(g$gene_short_name, levels = genes)

  # group_outline = TRUE

  p <- ggplot(aes(x, y), data = g)
  p <- p +
    ggplot2::geom_path(aes(x, y, group = edge_name), colour = con_colour, data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")


  # if (is.null(label_subset)) {
  #   label_subset = unique(node_metadata$group_nodes_by)
  # }
  #
  # label_subset = label_subset[label_subset != unlabeled_groups]
  group_nodes_by <- NULL
  if (is.null(group_nodes_by) == FALSE) {
    p <- p + ggforce::geom_mark_rect(
      aes(x, y,
        fill = gene_short_name,
        color = group_nodes_by,
        filter = group_nodes_by %in% unlabeled_groups == FALSE
      ),
      size = 0.5,
      expand = unit(2, "mm"),
      label.buffer = unit(1, "mm"),
      radius = unit(1.5, "mm"),
      label.margin = margin(1, 1, 1, 1, "mm"),
      label.fontsize = group_label_font_size,
      label.fontface = "plain",
      con.linetype = label_conn_linetype,
      con.colour = con_colour,
      show.legend = F
    ) #+
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

    p <- p + scale_fill_manual(values = rep("white", length(unique(g$gene_short_name))))
  }
  p <- p + guides(fill = "none")
  p <- p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(
      data = g %>% filter(gene_expr),
      aes(x, y,
        size = fraction_max,
        color = I(con_colour)
      )
    ) +
    ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(
      data = g %>% filter(gene_expr & fraction_max > 0),
      aes(x, y,
        size = fraction_max,
        color = get(color_nodes_by)
      )
    ) +
    labs(color = expression_legend_label) +
    viridis::scale_color_viridis(option = "viridis")

  if (is.null(edge_labels) == FALSE) {
    p <- p + ggnetwork::geom_nodetext(
      data = label_df,
      aes(x, y, label = label), size = 3
    )
  }

  p <- p + facet_wrap(~gene_short_name)

  p <- p +
    scale_size_continuous(labels = scales::percent, range = c(min_node_size, max_node_size)) +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position = legend_position) + guides(fill = "none")
  return(p)
}

#' Collect PSG Node Metadata
#'
#' This function collects metadata for nodes in a cell clustering structure (CCS).
#'
#' @param ccs A cell clustering structure object containing metadata and column data.
#' @param color_nodes_by A string specifying the metadata column to color nodes by.
#' @param label_nodes_by A string specifying the metadata column to label nodes by.
#' @param group_nodes_by A string specifying the metadata column to group nodes by.
#'
#' @return A data frame containing node metadata with columns for node ID, color, group, and label.
#'
#' @details
#' The function extracts unique cell groups from the CCS metadata and constructs a data frame
#' with node IDs. It then adds metadata columns for coloring, grouping, and labeling nodes based
#' on the specified parameters. If a parameter is not provided, the corresponding metadata column
#' is not added. The resulting data frame is returned with distinct rows and node IDs as row names.
#'
#' @import dplyr
#' @import tibble
#' @importFrom rlang sym
#'
#' @examples
#' \dontrun{
#' ccs <- load_ccs_data() # Assuming a function to load CCS data
#' node_metadata <- collect_psg_node_metadata(ccs, "color_column", "label_column", "group_column")
#' }
#'
#' @noRd
collect_psg_node_metadata <- function(ccs,
                                      color_nodes_by,
                                      label_nodes_by,
                                      group_nodes_by) {
  cell_groups <- ccs@metadata[["cell_group_assignments"]] %>%
    pull(cell_group) %>%
    unique()
  node_metadata <- tibble(id = cell_groups)

  metadata_cols <- c(
    color_nodes_by,
    group_nodes_by
  )
  if (is.null(label_nodes_by) == FALSE && label_nodes_by != "cell_group") {
    metadata_cols <- c(metadata_cols, label_nodes_by)
  }

  # G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata <- ccs@cds_coldata[, metadata_cols, drop = F] %>%
    as.data.frame()
  cell_group_metadata$cell_group <- ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  if (is.null(color_nodes_by) == FALSE) {
    color_by_metadata <- cell_group_metadata[, c("cell_group", color_nodes_by)] %>%
      as.data.frame() %>%
      dplyr::count(cell_group, !!sym(color_nodes_by)) %>%
      group_by(cell_group) %>%
      slice_max(n, with_ties = FALSE) %>%
      dplyr::select(-n)
    colnames(color_by_metadata) <- c("cell_group", "color_nodes_by")
    node_metadata <- left_join(node_metadata, color_by_metadata, by = c("id" = "cell_group"))
  }
  if (is.null(group_nodes_by) == FALSE) {
    group_by_metadata <- cell_group_metadata[, c("cell_group", group_nodes_by)] %>%
      as.data.frame() %>%
      dplyr::count(cell_group, !!sym(group_nodes_by)) %>%
      group_by(cell_group) %>%
      slice_max(n, with_ties = FALSE) %>%
      dplyr::select(-n)
    colnames(group_by_metadata) <- c("cell_group", "group_nodes_by")
    node_metadata <- left_join(node_metadata, group_by_metadata, by = c("id" = "cell_group"))
  }
  if (is.null(label_nodes_by) == FALSE) {
    label_by_metadata <- cell_group_metadata[, c("cell_group", label_nodes_by), drop = F]
    colnames(label_by_metadata) <- c("cell_group", "label_nodes_by")
    label_by_metadata <- label_by_metadata %>%
      as.data.frame() %>%
      dplyr::count(cell_group, label_nodes_by) %>%
      group_by(cell_group) %>%
      slice_max(n, with_ties = FALSE) %>%
      dplyr::select(-n)
    node_metadata <- left_join(node_metadata, label_by_metadata, by = c("id" = "cell_group"))
  } else {
    node_metadata$label_nodes_by <- node_metadata$id
  }

  node_metadata <- node_metadata %>%
    distinct() %>%
    as.data.frame(stringsAsFactor = FALSE)
  row.names(node_metadata) <- node_metadata$id

  return(node_metadata)
}

#' Connect Isolated Nodes in a Graph
#'
#' This function takes a graph and node metadata, and connects isolated nodes within each group specified in the metadata.
#'
#' @param G An igraph object representing the graph.
#' @param node_metadata A data frame containing node metadata. It must have columns "id" and "group_nodes_by".
#'
#' @return An igraph object with isolated nodes connected within each group.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks if the node metadata contains the required columns "id" and "group_nodes_by".
#'   \item Creates a mapping from node id to group.
#'   \item Iterates over each group and identifies isolated nodes within the group.
#'   \item Connects isolated nodes in a grid-like structure within each group.
#' }
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' # Create a sample graph
#' G <- make_empty_graph(n = 10, directed = FALSE)
#'
#' # Create sample node metadata
#' node_metadata <- data.frame(
#'   id = 1:10,
#'   group_nodes_by = rep(1:2, each = 5)
#' )
#'
#' # Connect isolated nodes
#' G_aug <- connect_isolated_nodes(G, node_metadata)
#'
#' # Plot the augmented graph
#' plot(G_aug)
#' }
#'
#' @import igraph
#' @import dplyr
connect_isolated_nodes <- function(G, node_metadata) {
  # Ensure node_metadata has columns "id" and "group_nodes_by"
  if (!all(c("id", "group_nodes_by") %in% colnames(node_metadata))) {
    stop("node_metadata must have 'id' and 'group_nodes_by' columns.")
  }

  # Create a mapping from node id to group
  node_group_map <- node_metadata %>%
    setNames(c("node", "group")) %>%
    as.list()

  # Get the list of groups
  groups <- unique(node_metadata$group_nodes_by)

  G_aug <- G

  for (group in groups) {
    # Extract nodes belonging to the current group
    group_nodes <- node_metadata %>%
      filter(group_nodes_by == group) %>%
      pull(id)

    # Induce subgraph for the current group
    G_sub <- igraph::induced_subgraph(G_aug, v = group_nodes)

    # Identify isolated nodes in the subgraph
    isolated_nodes <- igraph::V(G_sub)[igraph::degree(G_sub) == 0]$name

    if (length(isolated_nodes) > 1) {
      # Find connected components in the subgraph
      components <- igraph::decompose.graph(G_sub)

      # Calculate the number of rows and columns for the grid
      num_nodes <- length(isolated_nodes)

      # Determine the height of the grid based on the tallest connected component
      max_component_size <- max(sapply(components, igraph::diameter))
      if (max_component_size == 0) {
        num_rows <- floor(sqrt(num_nodes))
      } else {
        num_rows <- ceiling(num_nodes / max_component_size)
      }

      num_cols <- max(1, num_nodes %/% num_rows)

      # Create a grid structure to connect isolated nodes
      for (i in seq_len(num_nodes - 1)) {
        row_i <- (i %/% num_cols) + 1
        col_i <- (i %% num_cols) + 1

        row_j <- ((i + 1) %/% num_cols) + 1
        col_j <- ((i + 1) %% num_cols) + 1

        if (row_i == row_j && col_j == col_i + 1) {
          G_aug <- igraph::add_edges(G_aug, edges = c(isolated_nodes[i], isolated_nodes[i + 1]))
        } else if (row_i == row_j + 1 && col_j == col_i) {
          G_aug <- igraph::add_edges(G_aug, edges = c(isolated_nodes[i], isolated_nodes[i + 1]))
        }
      }
    }
  }

  return(G_aug)
}

assign_nodes_to_layers <- function(G, num_layers, node_metadata) {
  # Ensure node_metadata has columns "id" and "group_nodes_by"
  if (!all(c("id", "group_nodes_by") %in% colnames(node_metadata))) {
    stop("node_metadata must have 'id' and 'group_nodes_by' columns.")
  }

  # Create a mapping from node id to group
  node_group_map <- setNames(node_metadata$group_nodes_by, node_metadata$id)

  # Get the list of groups
  groups <- unique(node_metadata$group_nodes_by)

  # Create a list to hold nodes for each group
  group_nodes <- vector("list", length(groups))
  names(group_nodes) <- groups

  # Populate the group_nodes list with nodes from node_metadata
  for (group in groups) {
    group_nodes[[group]] <- node_metadata %>%
      filter(group_nodes_by == group) %>%
      pull(id)
  }

  # Calculate the total number of nodes
  total_nodes <- vcount(G)

  # Calculate the target number of nodes per layer
  target_nodes_per_layer <- ceiling(total_nodes / num_layers)

  # Create a list to hold layers
  layers <- vector("list", num_layers)

  current_layer <- 1
  node_count_in_current_layer <- 0

  for (group in groups) {
    group_node_ids <- unlist(group_nodes[[group]])

    # Check if adding this group would exceed the target number of nodes per layer
    if (node_count_in_current_layer + length(group_node_ids) > target_nodes_per_layer && current_layer < num_layers) {
      current_layer <- current_layer + 1
      node_count_in_current_layer <- 0
    }

    # Add the group to the current layer
    layers[[current_layer]] <- c(layers[[current_layer]], group_node_ids)
    node_count_in_current_layer <- node_count_in_current_layer + length(group_node_ids)
  }

  return(layers)
}

#'
#' @examples
#' # Example usage:
#' # G_with_hidden <- igraph::make_graph(...)
#' ##' Connect Hidden Nodes for Layers
#'
#' This function connects hidden nodes between layers in a given graph.
#'
#' @param G_with_hidden An igraph object representing the graph with hidden nodes.
#' @param layers A list of vectors, where each vector contains the node IDs for a specific layer.
#'
#' @return An igraph object with added connections between hidden nodes of consecutive layers.
#'
#' @details
#' The function iterates through the layers and connects the tail nodes of the current layer to the head nodes of the next layer.
#' It uses the `make_ego_graph` function to get the subgraph for each layer and its hidden nodes.
#' The tail nodes are identified by the prefix "tail_" and the head nodes by the prefix "head_".
#' The function ensures that each tail node is connected to a subset of head nodes in a round-robin fashion.
#'
#' @import igraph layers <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
#' # connected_graph <- connect_hidden_nodes_for_layers(G_with_hidden, layers)
#'
connect_hidden_nodes_for_layers <- function(G_with_hidden, layers) {
  connected_G_with_hidden <- G_with_hidden
  num_layers <- length(layers)

  # Add connections between tail nodes of one layer and head nodes of the next layer
  for (i in seq_len(num_layers - 1)) {
    # current_layer_subgraph <- igraph::induced_subgraph(graph = connected_G_with_hidden, vids = layers[[i]])
    current_layer_subgraph <- igraph::make_ego_graph(connected_G_with_hidden, order = 1, nodes = layers[[i]], mode = "out") # get current layer plus its hidden tail nodes
    current_layer_subgraph <- do.call(igraph::union, current_layer_subgraph)
    current_layer_tail_nodes <- igraph::V(current_layer_subgraph)[igraph::degree(current_layer_subgraph, mode = "out") == 0]$name
    current_layer_tail_nodes <- current_layer_tail_nodes[grepl("^tail_", current_layer_tail_nodes)] # select only the hidden tail nodes

    # next_layer_subgraph <- igraph::induced_subgraph(graph = connected_G_with_hidden, vids = layers[[i + 1]])
    next_layer_subgraph <- igraph::make_ego_graph(connected_G_with_hidden, order = 1, nodes = layers[[i + 1]], mode = "in") # get current layer plus its hidden tail nodes
    next_layer_subgraph <- do.call(igraph::union, next_layer_subgraph)
    next_layer_head_nodes <- igraph::V(next_layer_subgraph)[igraph::degree(next_layer_subgraph, mode = "in") == 0]$name
    next_layer_head_nodes <- next_layer_head_nodes[grepl("^head_", next_layer_head_nodes)] # select only the hidden head nodes

    # Convert vertex sequences to numeric IDs
    # current_layer_tail_ids <- as.numeric(current_layer_tail_nodes)
    # next_layer_head_ids <- as.numeric(next_layer_head_nodes)

    num_heads <- length(next_layer_head_nodes)

    for (j in seq_along(current_layer_tail_nodes)) {
      tail_n <- current_layer_tail_nodes[j]

      # Determine the indices of head nodes to connect
      start_idx <- ((j - 1) %% num_heads) + 1
      end_idx <- min(start_idx + 5, num_heads)

      selected_heads <- next_layer_head_nodes[start_idx:end_idx]

      for (head_n in selected_heads) {
        connected_G_with_hidden <- add_edges(connected_G_with_hidden, edges = c(tail_n, head_n))
      }
    }
  }

  return(connected_G_with_hidden)
}


#' Layout State Graph
#'
#' This function layouts a directed graph with optional node metadata and edge labels.
#' It ensures the graph is directed, connects isolated nodes, assigns nodes to layers,
#' and adds hidden head and tail nodes to each component. The function then performs
#' layout using Rgraphviz and extracts Bezier curve control points for edges.
#'
#' @param G An igraph object representing the directed graph.
#' @param node_metadata A data frame containing node metadata. Must include columns 'group_nodes_by' and 'id'.
#' @param edge_labels A named vector of edge labels (optional).
#' @param num_layers An integer specifying the number of layers for node assignment (default is 1).
#' @param weighted A logical value indicating whether the graph is weighted (default is FALSE).
#'
#' @return A list containing:
#' \item{gvizl_coords}{A matrix of layout coordinates for the nodes.}
#' \item{bezier_df}{A data frame of Bezier curve control points for edges.}
#' \item{label_df}{A data frame of edge labels (if provided).}
#' \item{grouping_df}{A data frame of node groupings based on metadata.}
#' \item{hidden_gvizl_coords}{A matrix of layout coordinates for hidden nodes.}
#' \item{hidden_bezier_df}{A data frame of Bezier curve control points for edges involving hidden nodes.}
#' \item{hidden_label_df}{A data frame of edge labels for edges involving hidden nodes (if provided).}
#'
#' @import igraph
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import Rgraphviz
#' @importFrom graph subGraph graphAM
#' @importFrom stringr str_split_fixed
#' @importFrom tibble tibble
#' @importFrom stats setNames
#' @importFrom utils head
#'
#' @examples
#' # Example usage:
#' G <- igraph::make_ring(10, directed = TRUE)
#' node_metadata <- data.frame(id = 1:10, group_nodes_by = rep(1:2, each = 5))
#' result <- layout_state_graph(G, node_metadata)
#'
#' @export
layout_state_graph <- function(G, node_metadata, edge_labels = NULL, num_layers = 1, weighted = FALSE) {
  make_subgraphs_for_groups <- function(subgraph_ids, G_nel) {
    sg_nodes <- subgraph_ids %>%
      pull(id) %>%
      as.character() %>%
      unique()
    sg <- list(graph = graph::subGraph(snodes = sg_nodes, graph = G_nel), cluster = TRUE)
    return(sg)
  }

  if (!igraph::is_directed(G)) {
    stop("The graph must be directed.")
  }

  G_orig <- G
  G <- connect_isolated_nodes(G, node_metadata)
  components <- igraph::decompose(G)
  num_components <- length(components)
  layers <- assign_nodes_to_layers(G, num_layers, node_metadata)

  G_with_hidden <- G
  head_nodes <- paste0("head_", seq_len(num_components))
  tail_nodes <- paste0("tail_", seq_len(num_components))

  hidden_nodes_tibble <- tibble(node_id = character(), component = integer(), group = character())

  for (i in seq_len(num_components)) {
    comp <- components[[i]]
    root_nodes <- V(comp)[igraph::degree(comp, mode = "in") == 0]$name
    if (length(root_nodes) == 0) {
      root_nodes <- V(comp)[1]$name
    }
    leaf_nodes <- V(comp)[igraph::degree(comp, mode = "out") == 0]$name
    if (length(leaf_nodes) == 0) {
      leaf_nodes <- V(comp)[1]$name
    }

    G_v_names <- igraph::V(G_with_hidden)$name
    G_with_hidden <- add_vertices(G_with_hidden, 2)
    igraph::V(G_with_hidden)$name <- c(G_v_names, head_nodes[i], tail_nodes[i])

    for (root in root_nodes) {
      G_with_hidden <- add_edges(G_with_hidden, edges = c(head_nodes[i], as.character(root)))
    }
    for (leaf in leaf_nodes) {
      G_with_hidden <- add_edges(G_with_hidden, edges = c(as.character(leaf), tail_nodes[i]))
    }

    # Determine the group for the hidden nodes
    root_groups <- node_metadata %>%
      filter(id %in% root_nodes) %>%
      pull(group_nodes_by)
    leaf_groups <- node_metadata %>%
      filter(id %in% leaf_nodes) %>%
      pull(group_nodes_by)

    head_group <- names(sort(table(root_groups), decreasing = TRUE))[1]
    tail_group <- names(sort(table(leaf_groups), decreasing = TRUE))[1]

    hidden_nodes_tibble <- hidden_nodes_tibble %>%
      add_row(node_id = head_nodes[i], component = i, group = head_group) %>%
      add_row(node_id = tail_nodes[i], component = i, group = tail_group)
  }

  G_with_hidden <- connect_hidden_nodes_for_layers(G_with_hidden, layers)

  if (weighted) {
    G_with_hidden_nel <- graph::graphAM(asMatrix = as_adjacency_matrix(G_with_hidden, sparse = FALSE, attr = "weight"), edgemode = "directed") %>% as("graphNEL")
  } else {
    G_with_hidden_nel <- graph::graphAM(adjMat = as_adjacency_matrix(G_with_hidden, sparse = FALSE), edgemode = "directed") %>% as("graphNEL")
  }

  if (is.null(node_metadata)) {
    subgraphs <- NULL
  } else {
    subgraph_df <- node_metadata %>%
      select(group_nodes_by, id) %>%
      group_by(group_nodes_by) %>%
      tidyr::nest(subgraph_ids = id) %>%
      summarize(subgraph = purrr::map(
        .f = purrr::possibly(make_subgraphs_for_groups, NULL),
        .x = subgraph_ids,
        G_with_hidden_nel
      ))
    subgraphs <- subgraph_df$subgraph
    names(subgraphs) <- subgraph_df$group_nodes_by
  }

  gvizl <- Rgraphviz::layoutGraph(G_with_hidden_nel, layoutType = "dot", subGList = subgraphs, recipEdges = "distinct")
  gvizl_coords <- cbind(gvizl@renderInfo@nodes$nodeX, gvizl@renderInfo@nodes$nodeY)

  beziers <- lapply(gvizl@renderInfo@edges$splines, function(bc) {
    bc_segments <- lapply(bc, Rgraphviz::bezierPoints)
    bezier_cp_df <- do.call(rbind, bc_segments) %>% as.data.frame()
    colnames(bezier_cp_df) <- c("x", "y")
    bezier_cp_df
  })
  bezier_df <- do.call(rbind, beziers)
  # bezier_df$edge_name = stringr::str_split_fixed(names(gvizl@renderInfo@edges$splines), "\\.", 2)[,1]
  bezier_df$edge_name <- str_split(row.names(bezier_df), "\\.[0-9]+$", simplify = TRUE)[, 1]
  # bezier_df$edge_name = stringr::str_split_fixed(row.names(bezier_df), "\\.", 2)[,1]
  bezier_df$from <- stringr::str_split_fixed(bezier_df$edge_name, "~", 2)[, 1]
  bezier_df$to <- stringr::str_split_fixed(bezier_df$edge_name, "~", 2)[, 2]
  bezier_df <- left_join(bezier_df, tibble(edge_name = names(gvizl@renderInfo@edges$direction), edge_direction = gvizl@renderInfo@edges$direction))
  bezier_df <- bezier_df %>% dplyr::distinct()
  
  
  bezier_df = left_join(bezier_df, igraph::as_data_frame(G) %>% select(-weight), by = c("from", "to"))

  if (!is.null(edge_labels)) {
    label_df <- data.frame(
      edge_id = gsub("\\.", "~", names(gvizl@renderInfo@edges$splines)),
      label = unname(edge_labels[edge_names])
    )
  } else {
    label_df <- NULL
  }

  hidden_nodes <- c(head_nodes, tail_nodes)
  gvizl_coords_clean <- gvizl_coords[!rownames(gvizl_coords) %in% hidden_nodes, ]
  gvizl_coords_hidden <- gvizl_coords[rownames(gvizl_coords) %in% hidden_nodes, ]

  hidden_gvizl_coords_tibble <- tibble(
    x = gvizl_coords_hidden[, 1],
    y = gvizl_coords_hidden[, 2],
    name = rownames(gvizl_coords_hidden)
  )
  hidden_gvizl_coords_tibble <- hidden_gvizl_coords_tibble %>% left_join(hidden_nodes_tibble, by = c("name" = "node_id"))

  G_orig_edgelist <- igraph::get.edgelist(G_orig)
  colnames(G_orig_edgelist) <- c("from", "to")
  G_orig_edgelist <- G_orig_edgelist %>% as_tibble()

  bezier_df_clean <- bezier_df %>%
    filter(from %in% hidden_nodes == FALSE & to %in% hidden_nodes == FALSE)
  bezier_df_clean <- bezier_df_clean %>% inner_join(G_orig_edgelist, by = c("from" = "from", "to" = "to"))

  bezier_df_hidden <- bezier_df %>%
    filter(from %in% hidden_nodes | to %in% hidden_nodes)

  if (!is.null(label_df)) {
    label_df_clean <- label_df %>%
      filter(edge_id %in% bezier_df_clean$edge_id)
    label_df_hidden <- label_df %>%
      filter(edge_id %in% bezier_df_hidden$edge_id)
  } else {
    label_df_clean <- NULL
    label_df_hidden <- NULL
  }

  grouping_df <- node_metadata %>%
    select(group_nodes_by, id)

  return(list(
    gvizl_coords = gvizl_coords_clean,
    bezier_df = bezier_df_clean,
    label_df = label_df_clean,
    grouping_df = grouping_df,
    hidden_gvizl_coords = hidden_gvizl_coords_tibble,
    hidden_bezier_df = bezier_df_hidden,
    hidden_label_df = label_df_hidden
  ))
}



#' geom_richnodelabel
#'
#' This function creates a custom ggplot2 layer for rich text node labels in network plots.
#'
#' @param mapping Set of aesthetic mappings created by `aes()` or `aes_()`. If specified and `inherit.aes = TRUE` (the default), it is combined with the default mapping at the top level of the plot. You must supply `mapping` if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three options:
#'   - If `NULL`, the default, the data is inherited from the plot data as specified in the call to `ggplot()`.
#'   - A `data.frame`, or other object, will override the plot data. All objects will be fortified to produce a data frame. See `fortify()` for which variables will be created.
#'   - A `function` will be called with a single argument, the plot data. The return value must be a `data.frame`, and will be used as the layer data.
#' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param ... Other arguments passed on to `layer()`. These are often aesthetics, used to set an aesthetic to a fixed value, like `color = "red"` or `size = 3`. They may also be parameters to the paired geom/stat.
#' @param parse If `TRUE`, the labels will be parsed into expressions and displayed as described in `?plotmath`.
#' @param nudge_x Horizontal adjustment to nudge labels by. Useful for offsetting text from points, particularly on discrete scales.
#' @param nudge_y Vertical adjustment to nudge labels by. Useful for offsetting text from points, particularly on discrete scales.
#' @param label.padding Amount of padding around label. Defaults to `unit(0.25, "lines")`.
#' @param label.r Radius of rounded corners. Defaults to `unit(0.15, "lines")`.
#' @param label.size Size of label border. Defaults to `0.25`.
#' @param na.rm If `FALSE`, the default, missing values are removed with a warning. If `TRUE`, missing values are silently removed.
#' @param show.legend Logical. Should this layer be included in the legends? `NA`, the default, includes if any aesthetics are mapped. `FALSE` never includes, and `TRUE` always includes.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics, rather than combining with them. This is most useful for helper functions that define both data and aesthetics and shouldn't inherit behaviour from the default plot specification, e.g. `borders()`.
#'
#' @return A ggplot2 layer that can be added to a ggplot object.
geom_richnodelabel <- function(mapping = NULL, data = NULL, position = "identity",
                               ..., parse = FALSE, nudge_x = 0, nudge_y = 0, label.padding = unit(
                                 0.25,
                                 "lines"
                               ), label.r = unit(0.15, "lines"), label.size = 0.25,
                               na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`",
        call. = FALSE
      )
    }
    position <- ggplot2::position_nudge(nudge_x, nudge_y)
  }
  ggplot2::layer(
    data = data, mapping = mapping, stat = ggnetwork:::StatNodes,
    geom = ggtext:::GeomRichtext, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, params = list(
      parse = parse,
      label.padding = label.padding, label.r = label.r,
      label.size = label.size, na.rm = na.rm, ...
    )
  )
}

plot_state_graph_perturb_effects <- function(ccs,
                                             state_graph,
                                             perturbation_table,
                                             patterns = c("direct", "indirect", "inferred", "predicted"),
                                             num_top_perturbs = 3,
                                             num_top_genes = 3,
                                             color_nodes_by = NULL,
                                             label_nodes_by = NULL,
                                             group_nodes_by = NULL,
                                             label_edges_by = NULL,
                                             edge_weights = NULL,
                                             arrow.gap = 0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             node_size = 2,
                                             num_layers = 10,
                                             min_edge_size = 0.1,
                                             max_edge_size = 2,
                                             fract_expr = 0.0,
                                             mean_expr = 0.0,
                                             unlabeled_groups = c("Unknown"),
                                             label_groups = TRUE,
                                             hide_unlinked_nodes = TRUE,
                                             group_label_font_size = 6,
                                             edge_label_font_size = 2,
                                             label_conn_linetype = "dotted",
                                             legend_position = "none",
                                             con_colour = "darkgrey",
                                             group_outline = FALSE) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  node_support_df <- igraph::as_data_frame(state_graph, what = "vertices") %>%
    tibble::as_tibble() %>%
    rename(id = name)
  # node_metadata = node_metadata %>% tidyr::unnest(direct_perturb)

  perturbation_table <- perturbation_table %>%
    group_by(id) %>%
    mutate(perturb_display_name = case_when(
      perturb_effect == "direct" ~ glue::glue("<i style='color:#2B9EB3'>{perturb_name}</i>"),
      perturb_effect == "indirect" ~ glue::glue("<i style='color:#FCAB10'>{perturb_name}</i>"),
      perturb_effect == "predicted" ~ glue::glue("<i style='color:#F8333C'>{perturb_name}</i>"),
      TRUE ~ glue::glue("<i style='color:#44AF69'>{perturb_name}</i>")
    )) %>%
    summarize(perturb_effect_label = ifelse(n() > num_top_genes,
      paste0(c(perturb_display_name[1:num_top_genes], paste("+", n() - num_top_genes, " more", sep = "")), collapse = "<br>"),
      paste0(perturb_display_name, collapse = "<br>")
    ))
  # summarize(top_genes=ifelse(n() > num_top_perturbs,
  #                            paste0(c(gene_short_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "\n"),
  #                            paste0(gene_short_name, collapse = "\n"))) %>%
  #

  node_metadata <- node_metadata %>% left_join(perturbation_table %>% select(id, perturb_effect_label), by = c("id" = "id"))
  if (is.null(label_nodes_by)) {
    label_nodes_by <- "perturb_effect_label"
    node_metadata <- node_metadata %>% mutate(label_nodes_by = perturb_effect_label)
  } else {
    node_metadata <- node_metadata %>% mutate(label_nodes_by = ifelse(is.na(perturb_effect_label), label_nodes_by, paste(label_nodes_by, perturb_effect_label, sep = "<br>")))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, edge_labels, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
      unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
      edge_thickness = replace_na(edge_thickness, min_edge_size)
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge <- FALSE
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)
  # add assembly group
  bezier_df$assembly_group <- vapply(strsplit(bezier_df$from, "-", fixed = T), "[", "", 1)
  g$assembly_group <- vapply(strsplit(g$name, "-", fixed = T), "[", "", 1)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group = edge_name, size = edge_thickness, linetype = unsupported_edge), colour = con_colour, data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")

  # draw activator edges
  # ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE) {
    if (group_outline) {
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          col = group_nodes_by,
          label = group_nodes_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0.5,
        label.fontsize = group_label_font_size,
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        data = g
      )
    } else {
      if (label_groups) {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            label = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          expand = unit(2, "mm"),
          label.buffer = unit(1, "mm"),
          radius = unit(1.5, "mm"),
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = group_label_font_size,
          label.fontface = "plain",
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      } else {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          label.fontsize = group_label_font_size,
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      }
    }
  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE) {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          geom_richnodelabel(
            data = g,
            aes(x, y,
              fill = !!sym(color_nodes_by),
              label = label_nodes_by
            ),
            size = node_size
          ) +
          labs(fill = color_nodes_by)
        p <- p + scale_fill_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + geom_richnodelabel(
          data = g,
          aes(x, y,
            fill = color_nodes_by,
            label = label_nodes_by
          ),
          size = node_size
        ) +
          labs(fill = color_nodes_by)
      }
    } else {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(
            data = g,
            aes(x, y,
              color = !!sym(color_nodes_by)
            ),
            size = node_size
          ) +
          labs(color = color_nodes_by)
        p <- p + scale_color_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + ggnetwork::geom_nodes(
          data = g,
          aes(x, y,
            color = color_nodes_by
          ),
          size = node_size
        ) +
          labs(color = color_nodes_by)
      }
    }
  } else {
    if (is.null(label_nodes_by) == FALSE) {
      p <- p + geom_richnodelabel(
        data = g,
        aes(x, y,
          xend = xend, yend = yend,
          label = label_nodes_by
        ),
        size = node_size
      )
    } else {
      p <- p + ggnetwork::geom_nodes(
        data = g,
        aes(x, y, xend = xend, yend = yend),
        size = node_size
      )
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df <- layout_info$label_df
    # p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p <- p + geom_text(
      data = label_df,
      aes(x, y, label = label),
      size = edge_label_font_size
    )
  }

  p <- p + annotate(
    geom = "richtext", x = 0.15 * max(g$xend), y = 0.01 * max(g$yend),
    size = 2,
    fill = NA, label.color = NA,
    label = "<i style='color:#2B9EB3'>direct</i> <i style='color:#FCAB10'>indirect</i> <i style='color:#F8333C'>predicted</i> <i style='color:#44AF69'>other</i> "
  )
  p <- p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position = legend_position)
  return(p)
}


plot_state_graph_marker_genes <- function(ccs,
                                          state_graph,
                                          marker_table,
                                          patterns_for_marker_genes = c(
                                            "Specifically activated",
                                            "Selectively activated",
                                            "Selectively upregulated",
                                            "Specifically upregulated",
                                            "Specifically maintained",
                                            "Selectively maintained",
                                            "Activated",
                                            "Upregulated",
                                            "Maintained"
                                          ),
                                          num_top_genes = 3,
                                          color_nodes_by = NULL,
                                          label_nodes_by = NULL,
                                          group_nodes_by = NULL,
                                          label_edges_by = NULL,
                                          edge_weights = NULL,
                                          arrow.gap = 0.03,
                                          arrow_unit = 2,
                                          bar_unit = .075,
                                          node_size = 2,
                                          num_layers = 10,
                                          min_edge_size = 0.1,
                                          max_edge_size = 2,
                                          fract_expr = 0.0,
                                          mean_expr = 0.0,
                                          unlabeled_groups = c("Unknown"),
                                          label_groups = TRUE,
                                          hide_unlinked_nodes = TRUE,
                                          group_label_font_size = 6,
                                          edge_label_font_size = 2,
                                          label_conn_linetype = "dotted",
                                          legend_position = "none",
                                          con_colour = "darkgrey",
                                          group_outline = FALSE) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  node_support_df <- igraph::as_data_frame(state_graph, what = "vertices") %>%
    tibble::as_tibble() %>%
    rename(id = name)
  # node_metadata = node_metadata %>% tidyr::unnest(direct_perturb)

  red_patterns <- c(
    "Specifically activated",
    "Specifically upregulated",
    "Specifically maintained"
  )
  orange_patterns <- c(
    "Selectively activated",
    "Selectively upregulated",
    "Selectively maintained"
  )
  blue_patterns <- c(
    "Activated",
    "Upregulated"
  )
  green_patterns <- c("Maintained")

  patterns_for_marker_genes <- c(
    "Specifically activated",
    "Selectively activated",
    "Selectively upregulated",
    "Specifically upregulated",
    "Specifically maintained",
    "Selectively maintained",
    "Activated",
    "Upregulated",
    "Maintained"
  )


  marker_table <- marker_table %>%
    group_by(id) %>%
    mutate(marker_display_name = case_when(
      marker_type %in% blue_patterns ~ glue::glue("<i style='color:#2B9EB3'>{marker_name}</i>"),
      marker_type %in% orange_patterns ~ glue::glue("<i style='color:#FCAB10'>{marker_name}</i>"),
      marker_type %in% red_patterns ~ glue::glue("<i style='color:#F8333C'>{marker_name}</i>"),
      marker_type %in% green_patterns ~ glue::glue("<i style='color:#44AF69'>{marker_name}</i>"),
      TRUE ~ glue::glue("<i style='color:#D3D3D3'>{marker_name}</i>")
    )) %>%
    summarize(marker_label = ifelse(n() > num_top_genes,
      paste0(c(marker_display_name[1:num_top_genes], paste("+", n() - num_top_genes, " more", sep = "")), collapse = "<br>"),
      paste0(marker_display_name, collapse = "<br>")
    ))
  # summarize(top_genes=ifelse(n() > num_top_perturbs,
  #                            paste0(c(gene_short_name[1:num_top_genes], paste("+", n()-num_top_genes, " more", sep="")), collapse = "\n"),
  #                            paste0(gene_short_name, collapse = "\n"))) %>%
  #

  node_metadata <- node_metadata %>% left_join(marker_table %>% select(id, marker_label), by = c("id" = "id"))
  if (is.null(label_nodes_by)) {
    label_nodes_by <- "marker_label"
    node_metadata <- node_metadata %>% mutate(label_nodes_by = marker_label)
  } else {
    node_metadata <- node_metadata %>% mutate(label_nodes_by = ifelse(is.na(marker_label), label_nodes_by, paste(label_nodes_by, marker_label, sep = "<br>")))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, edge_labels, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
      unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
      edge_thickness = replace_na(edge_thickness, min_edge_size)
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge <- FALSE
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)
  # add assembly group
  bezier_df$assembly_group <- vapply(strsplit(bezier_df$from, "-", fixed = T), "[", "", 1)
  g$assembly_group <- vapply(strsplit(g$name, "-", fixed = T), "[", "", 1)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group = edge_name, size = edge_thickness, linetype = unsupported_edge), colour = con_colour, data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")

  # draw activator edges
  # ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE) {
    if (group_outline) {
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          col = group_nodes_by,
          label = group_nodes_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0.5,
        label.fontsize = group_label_font_size,
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        data = g
      )
    } else {
      if (label_groups) {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            label = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          expand = unit(2, "mm"),
          label.buffer = unit(1, "mm"),
          radius = unit(1.5, "mm"),
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = group_label_font_size,
          label.fontface = "plain",
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      } else {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          label.fontsize = group_label_font_size,
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      }
    }
  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE) {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          geom_richnodelabel(
            data = g,
            aes(x, y,
              fill = !!sym(color_nodes_by),
              label = label_nodes_by
            ),
            size = node_size
          ) +
          labs(fill = color_nodes_by)
        p <- p + scale_fill_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + geom_richnodelabel(
          data = g,
          aes(x, y,
            fill = color_nodes_by,
            label = label_nodes_by
          ),
          size = node_size
        ) +
          labs(fill = color_nodes_by)
      }
    } else {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(
            data = g,
            aes(x, y,
              color = !!sym(color_nodes_by)
            ),
            size = node_size
          ) +
          labs(color = color_nodes_by)
        p <- p + scale_color_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + ggnetwork::geom_nodes(
          data = g,
          aes(x, y,
            color = color_nodes_by
          ),
          size = node_size
        ) +
          labs(color = color_nodes_by)
      }
    }
  } else {
    if (is.null(label_nodes_by) == FALSE) {
      p <- p + geom_richnodelabel(
        data = g,
        aes(x, y,
          xend = xend, yend = yend,
          label = label_nodes_by
        ),
        size = node_size
      )
    } else {
      p <- p + ggnetwork::geom_nodes(
        data = g,
        aes(x, y, xend = xend, yend = yend),
        size = node_size
      )
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df <- layout_info$label_df
    # p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p <- p + geom_text(
      data = label_df,
      aes(x, y, label = label),
      size = edge_label_font_size
    )
  }

  p <- p + annotate(
    geom = "richtext", x = 0.15 * max(g$xend), y = 0.01 * max(g$yend),
    size = 2,
    fill = NA, label.color = NA,
    label = "<i style='color:#2B9EB3'>upregulated</i> <i style='color:#FCAB10'>selective</i> <i style='color:#F8333C'>specific</i> <i style='color:#44AF69'>maintained</i> "
  )
  p <- p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position = legend_position)
  return(p)
}


plot_state_graph_key_genes <- function(ccs,
                                       state_graph,
                                       gene_pattern_table,
                                       num_top_genes = 3,
                                       patterns = c(
                                         "Activated",
                                         "Upregulated",
                                         "Specifically activated",
                                         "Selectively activated",
                                         "Specifically activated",
                                         "Selectively upregulated",
                                         "Specifically upregulated",
                                         "Specifically maintained",
                                         "Selectively maintained"
                                       ),
                                       color_nodes_by = NULL,
                                       label_nodes_by = NULL,
                                       group_nodes_by = NULL,
                                       label_edges_by = NULL,
                                       edge_weights = NULL,
                                       arrow.gap = 0.03,
                                       arrow_unit = 2,
                                       bar_unit = .075,
                                       node_size = 2,
                                       num_layers = 10,
                                       min_edge_size = 0.1,
                                       max_edge_size = 2,
                                       fract_expr = 0.0,
                                       mean_expr = 0.0,
                                       unlabeled_groups = c("Unknown"),
                                       label_groups = TRUE,
                                       hide_unlinked_nodes = TRUE,
                                       group_label_font_size = 6,
                                       edge_label_font_size = 2,
                                       label_conn_linetype = "dotted",
                                       legend_position = "none",
                                       con_colour = "darkgrey",
                                       group_outline = FALSE) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  # edges = hooke:::distance_to_root(edges)
  edges <- edges %>% dplyr::ungroup()

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)) {
    edges <- edges %>% select(from, to)
  } else {
    edges <- edges %>% select(from, to, weight = !!sym(edge_weights))
  }

  selected_genes <- gene_pattern_table %>%
    filter(interpretation %in% patterns) %>%
    group_by(cell_state) %>%
    arrange(desc(pattern_match_score)) %>%
    summarize(top_genes = ifelse(n() > num_top_genes,
      paste0(c(gene_short_name[1:num_top_genes], paste("+", n() - num_top_genes, " more", sep = "")), collapse = "\n"),
      paste0(gene_short_name, collapse = "\n")
    ))

  node_metadata <- node_metadata %>% left_join(selected_genes %>% select(cell_state, top_genes), by = c("id" = "cell_state"))
  if (is.null(label_nodes_by)) {
    label_nodes_by <- "top_genes"
    node_metadata <- node_metadata %>% mutate(label_nodes_by = top_genes)
  } else {
    node_metadata <- node_metadata %>% mutate(label_nodes_by = paste(label_nodes_by, top_genes, sep = "\n"))
  }

  if (is.null(label_edges_by)) {
    edge_info <- edges %>% select(from, to)
  } else {
    if (is(state_graph, "igraph")) {
      edge_info <- state_graph %>%
        igraph::as_data_frame() %>%
        select(from, to, label = !!sym(label_edges_by))
    } else {
      edge_info <- state_graph %>% select(from, to, label = !!sym(label_edges_by))
    }

    edges <- edges %>% left_join(edge_info)
    # print(edges)
  }

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE) {
    G_df <- igraph::as_data_frame(G)
    edge_names <- stringr::str_c(G_df$from, G_df$to, sep = "~")
    edge_labels <- igraph::E(G)$label
    names(edge_labels) <- edge_names
    # print(edge_labels)
    # edge_labels = NULL
  } else {
    edge_labels <- NULL
  }

  layout_info <- layout_state_graph(G, node_metadata, edge_labels, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df
  if (is.null(edge_weights) == FALSE) {
    bezier_df <- left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
  }

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)
  # add assembly group
  bezier_df$assembly_group <- vapply(strsplit(bezier_df$from, "-", fixed = T), "[", "", 1)
  g$assembly_group <- vapply(strsplit(g$name, "-", fixed = T), "[", "", 1)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group = edge_name, size = edge_thickness), colour = con_colour, data = bezier_df %>% distinct(), arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre")

  # draw activator edges
  # ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE) {
    if (group_outline) {
      p <- p + ggforce::geom_mark_rect(
        aes(x, y,
          col = group_nodes_by,
          label = group_nodes_by,
          filter = group_nodes_by %in% unlabeled_groups == FALSE
        ),
        size = 0.5,
        label.fontsize = group_label_font_size,
        con.linetype = label_conn_linetype,
        con.colour = con_colour,
        data = g
      )
    } else {
      if (label_groups) {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            label = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          expand = unit(2, "mm"),
          label.buffer = unit(1, "mm"),
          radius = unit(1.5, "mm"),
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = group_label_font_size,
          label.fontface = "plain",
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      } else {
        p <- p + ggforce::geom_mark_rect(
          aes(x, y,
            fill = group_nodes_by,
            filter = group_nodes_by %in% unlabeled_groups == FALSE
          ),
          size = 0,
          label.fontsize = group_label_font_size,
          con.linetype = label_conn_linetype,
          con.colour = con_colour,
          data = g
        )
      }
    }
  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE) {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodelabel(
            data = g,
            aes(x, y,
              fill = !!sym(color_nodes_by),
              label = label_nodes_by
            ),
            size = node_size
          ) +
          labs(fill = color_nodes_by)
        p <- p + scale_fill_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + ggnetwork::geom_nodelabel(
          data = g,
          aes(x, y,
            fill = color_nodes_by,
            label = label_nodes_by
          ),
          size = node_size
        ) +
          labs(fill = color_nodes_by)
      }
    } else {
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p <- p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(
            data = g,
            aes(x, y,
              color = !!sym(color_nodes_by)
            ),
            size = node_size
          ) +
          labs(color = color_nodes_by)
        p <- p + scale_color_gradient2(low = "royalblue3", mid = "white", high = "orangered3")
      } else {
        # if categorical
        p <- p + ggnetwork::geom_nodes(
          data = g,
          aes(x, y,
            color = color_nodes_by
          ),
          size = node_size
        ) +
          labs(color = color_nodes_by)
      }
    }
  } else {
    if (is.null(label_nodes_by) == FALSE) {
      p <- p + ggnetwork::geom_nodelabel(
        data = g,
        aes(x, y,
          xend = xend, yend = yend,
          label = label_nodes_by
        ),
        size = node_size
      )
    } else {
      p <- p + ggnetwork::geom_nodes(
        data = g,
        aes(x, y, xend = xend, yend = yend),
        size = node_size
      )
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df <- layout_info$label_df
    # p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p <- p + geom_text(
      data = label_df,
      aes(x, y, label = label),
      size = edge_label_font_size
    )
  }

  p <- p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position = legend_position)
  return(p)
}



# cell_type_ann_df = data.frame(cell_group = unique(colData(cds)[[group_cells_by]]),
#                               category = rep(1:7, each=2)) %>%
#   tibble::column_to_rownames("cell_group")


#' Plot Gene Expression Across Cell Groups
#'
#' This function generates a heatmap to visualize the expression of specified genes across 
#' different cell groups. It allows filtering of cells and genes, and provides options for 
#' annotating rows and columns of the heatmap.
#'
#' @param cds A cell_data_set object containing single-cell data.
#' @param group_cells_by A string specifying the column in `colData(cds)` to group cells by.
#' @param genes A character vector of gene names to include in the plot. If NULL, all genes are included.
#' @param cell_types A character vector of cell types to include in the plot. If NULL, all cell types are included.
#' @param gene_ann_df A data frame containing annotations for genes. Must have a column named "genes".
#' @param cell_type_ann_df A data frame containing annotations for cell groups. Must have a column named "cell_group".
#' @param drop_zeroes A logical value indicating whether to drop rows and columns with all zero values. Default is TRUE.
#'
#' @return A heatmap object generated by `pheatmap::pheatmap`.
#'
#' @importFrom dplyr filter select group_by top_n left_join
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#' @importFrom viridis viridis
#'
#' @examples
#' # Example usage:
#' plot_genes_expr_across_cell_groups(cds, group_cells_by = "cell_type", genes = c("GeneA", "GeneB"))
plot_genes_expr_across_cell_groups <- function(cds,
                                               group_cells_by,
                                               genes = NULL,
                                               cell_types = NULL,
                                               gene_ann_df = NULL,
                                               cell_type_ann_df = NULL,
                                               drop_zeroes = TRUE) {
  if (is.null(cell_types) == FALSE) {
    cds <- cds[, colData(cds)[[group_cells_by]] %in% cell_types]
  }

  # maybe warn if this isnt null / longer than a certain length
  if (is.null(genes) == FALSE) {
    cds <- cds[rowData(cds)$gene_short_name %in% genes, ]
  }

  agg_expr_data <- hooke:::aggregated_expr_data(cds, group_cells_by = group_cells_by)


  if (is.null(gene_ann_df) == FALSE) {
    gene_ann_df <- gene_ann_df %>% tibble::column_to_rownames("genes")
  }


  if (is.null(cell_type_ann_df) == FALSE) {
    cell_type_ann_df <- cell_type_ann_df %>% tibble::column_to_rownames("cell_group")
  }

  agg_expr_mat <- agg_expr_data %>%
    filter(gene_short_name %in% genes) %>%
    select(cell_group, gene_short_name, mean_expression) %>%
    pivot_wider(names_from = gene_short_name, values_from = mean_expression) %>%
    tibble::column_to_rownames("cell_group") %>%
    t()

  if (drop_zeroes) {
    agg_expr_mat <- agg_expr_mat[
      rowSums(agg_expr_mat) != 0,
      colSums(agg_expr_mat) != 0
    ]
  }


  p <- pheatmap::pheatmap(agg_expr_mat,
    cluster_rows = T,
    cluster_cols = T,
    color = viridis::viridis(10),
    annotation_row = gene_ann_df,
    annotation_col = cell_type_ann_df
  )
  return(p)
}



#' @title Get Maximum Abundance Contrast
#' 
#' @description This function calculates the maximum abundance contrast for a given 
#' cell clustering model (ccm) over a specified time interval. It allows for 
#' customization of the start and stop times, control abundances, and additional 
#' data inputs.
#' 
#' @param ccm A cell clustering model object containing cell clustering data.
#' @param start_time (Optional) The starting timepoint for the analysis. Defaults 
#' to the minimum timepoint in the data if not provided.
#' @param stop_time (Optional) The ending timepoint for the analysis. Defaults 
#' to the maximum timepoint in the data if not provided.
#' @param ctrl_abundances (Optional) A data frame of control abundances. If not 
#' provided, control abundances will be calculated based on the specified time 
#' interval.
#' @param interval_col A string specifying the column name in the cell clustering 
#' model that contains the timepoint information. Defaults to "timepoint".
#' @param newdata A tibble containing additional data to be used in the analysis. 
#' Defaults to an empty tibble.
#' 
#' @details The function first determines the range of timepoints to analyze 
#' based on the provided or default start and stop times. It then calculates 
#' perturbation effects and filters control abundances to match the timepoints 
#' present in the perturbation effects.
#' 
#' @return A data frame or tibble containing the maximum abundance contrast 
#' results.
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#' 
#' @examples
#' # Example usage:
#' # result <- get_max_abundance_contrast(ccm, start_time = 0, stop_time = 10)
#' 
get_max_abundance_contrast <- function(ccm,
                                       start_time = NULL,
                                       stop_time = NULL,
                                       ctrl_abundances = NULL,
                                       interval_col = "timepoint",
                                       newdata = tibble()) {
  timepoints <- as.numeric(unique(colData(ccm@ccs@cds)[[interval_col]]))
  timepoints <- timepoints[!is.na(timepoints)]

  if (is.null(start_time)) {
    start_time <- min(timepoints)
  }
  if (is.null(stop_time)) {
    stop_time <- max(timepoints)
  }

  if (nrow(newdata) > 0) {
    newdata <- cross_join(tibble(knockout = FALSE), newdata)
  } else {
    newdata <- tibble(knockout = FALSE)
  }

  # CHECK THIS
  perturb_effects <- platt:::get_perturbation_effects(ccm, ...)


  if (is.null(ctrl_abundances)) {
    ctrl_abundances <- get_extant_cell_types(ccm, start = start_time, stop = stop_time, newdata = newdata)
  } else {
    ctrl_abundances <- ctrl_abundances %>% filter(timepoint %in% unique(perturb_effects$time))
  }


  peak_abundances <- ctrl_abundances %>%
    filter(!is.na(percent_cell_type_range)) %>%
    group_by(cell_group) %>%
    top_n(n = 1, wt = percent_max_abund) %>%
    select(timepoint, cell_group)

  peak_abundances$timepoint <- as.character(peak_abundances$timepoint)

  max_abund_tbl <- left_join(peak_abundances, perturb_effects, by = c("cell_group", "timepoint" = "time"))

  return(max_abund_tbl)
}




plot_graph_simple <- function(ccs,
                              state_graph,
                              num_degs,
                              color_nodes_by = NULL,
                              label_nodes_by = NULL,
                              group_nodes_by = NULL,
                              arrow.gap = 0.03,
                              arrow_unit = 2,
                              bar_unit = .075,
                              node_size = 2,
                              min_edge_size = 0.1,
                              max_edge_size = 2,
                              con_colour = "darkgrey",
                              legend_position = "none",
                              hide_unlinked_nodes = TRUE,
                              flip_x = F) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  edges <- edges %>% dplyr::ungroup()
  edges <- edges %>% filter(from %in% node_metadata$id & to %in% node_metadata$id)

  G <- edges %>%
    distinct() %>%
    igraph::graph_from_data_frame(directed = T, vertices = node_metadata)

  layout_info <- layout_state_graph(G, node_metadata, NULL, weighted = FALSE)
  gvizl_coords <- layout_info$gvizl_coords
  bezier_df <- layout_info$bezier_df




  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale = F)
  g <- left_join(g, num_degs, by = c("name" = "cell_group"))

  if (flip_x) {
    g$x <- -1 * g$x
    bezier_df$x <- -1 * bezier_df$x
  }

  bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
  bezier_df$unsupported_edge <- FALSE

  g <- g %>% filter(!is.na(term))
  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group = edge_name, size = edge_thickness, linetype = unsupported_edge),
      colour = con_colour, data = bezier_df %>% distinct(),
      arrow = arrow(angle = 30, length = unit(arrow_unit, "pt"), type = "closed"), linejoin = "mitre"
    )

  p <- p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodes(
      data = g,
      aes(x, y,
        color = log10(n),
        size = node_size
      )
    ) +
    labs(color = "number degs") +
    scale_color_viridis_c()

  p <- p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position = legend_position)

  p <- p + facet_wrap(~term)
  return(p)
}
