#' Plot a Clustered Heatmap for a DEG Gene List
#'
#' Build a heatmap from a DEG table for a supplied set of genes, typically a
#' set of common DEGs shared across perturbations or cell groups. Cell groups
#' are retained only if they contain at least one gene from `genes` that is
#' detected, significant, and above the specified mean-expression threshold.
#' The retained gene-by-cell-group log fold-change matrix is then hierarchically
#' clustered along both axes and plotted with significant tiles outlined.
#'
#' @param deg_table A data frame containing DEG results. By default it must
#'   include gene identifiers, gene labels, cell-group labels, shrunken log fold
#'   changes, p-values, a detection flag, and log mean expression values.
#' @param genes Character vector of genes or IDs to include in the heatmap.
#'   Matching is performed against `gene_col`.
#' @param cell_groups Optional character vector of cell groups to include. If
#'   `NULL`, cell groups are selected automatically by requiring at least one
#'   gene from `genes` to satisfy the detection, significance, and mean
#'   expression thresholds.
#' @param sig_p_value_thresh Numeric significance threshold applied to
#'   `p_value_col` when selecting cell groups and when outlining significant
#'   tiles in the plot.
#' @param min_log_mean_expression Numeric minimum value in
#'   `log_mean_expression_col` required for a DEG to count toward retaining a
#'   cell group.
#' @param gene_col Name of the column containing the gene identifiers used to
#'   match `genes`.
#' @param gene_label_col Name of the column used for y-axis labels in the
#'   heatmap.
#' @param cell_group_col Name of the column containing cell-group labels.
#' @param lfc_col Name of the column containing log fold-change values used for
#'   tile fill.
#' @param p_value_col Name of the column containing DEG p-values.
#' @param present_col Name of the logical column indicating whether the gene is
#'   detected above the expression threshold for that cell group.
#' @param log_mean_expression_col Name of the column containing mean expression
#'   values on a log scale.
#'
#' @return A `ggplot` object showing clustered log fold changes for the selected
#'   genes across retained cell groups.
#'
#' @examples
#' \dontrun{
#' plot_deg_heatmap(
#'     deg_table = degs,
#'     genes = common_degs,
#'     cell_groups = c("Neurons", "Glia"),
#'     sig_p_value_thresh = 0.05,
#'     min_log_mean_expression = -2
#' )
#' }
#' @export
plot_deg_heatmap <- function(deg_table,
                             genes,
                             cell_groups = NULL,
                             sig_p_value_thresh = 0.05,
                             min_log_mean_expression = -2,
                             gene_col = "id",
                             gene_label_col = "gene_short_name",
                             cell_group_col = "cell_group",
                             lfc_col = "perturb_to_ctrl_shrunken_lfc",
                             p_value_col = "perturb_to_ctrl_p_value",
                             present_col = "present_above_thresh",
                             log_mean_expression_col = "log_mean_expression") {
    required_cols <- c(
        gene_col,
        gene_label_col,
        cell_group_col,
        lfc_col,
        p_value_col,
        present_col,
        log_mean_expression_col
    )
    missing_cols <- setdiff(required_cols, names(deg_table))
    if (length(missing_cols) > 0) {
        stop(
            "Missing required columns in `deg_table`: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }
    if (length(genes) == 0) {
        stop("`genes` must contain at least one gene.", call. = FALSE)
    }

    deg_subset <- deg_table %>%
        dplyr::filter(.data[[gene_col]] %in% genes)

    if (nrow(deg_subset) == 0) {
        stop("No rows in `deg_table` matched `genes`.", call. = FALSE)
    }

    if (is.null(cell_groups)) {
        cell_groups <- deg_subset %>%
            dplyr::filter(
                .data[[present_col]],
                .data[[p_value_col]] < sig_p_value_thresh,
                .data[[log_mean_expression_col]] > min_log_mean_expression
            ) %>%
            dplyr::pull(.data[[cell_group_col]]) %>%
            unique()

        if (length(cell_groups) == 0) {
            stop("No cell groups had at least one qualifying significant DEG.", call. = FALSE)
        }
    } else {
        cell_groups <- unique(as.character(cell_groups))
        if (length(cell_groups) == 0) {
            stop("`cell_groups` must contain at least one cell group when provided.", call. = FALSE)
        }

        missing_cell_groups <- setdiff(
            cell_groups,
            unique(as.character(deg_subset[[cell_group_col]]))
        )
        if (length(missing_cell_groups) > 0) {
            stop(
                "Requested `cell_groups` not found in filtered `deg_table`: ",
                paste(missing_cell_groups, collapse = ", "),
                call. = FALSE
            )
        }
    }

    plot_df <- deg_subset %>%
        dplyr::filter(.data[[cell_group_col]] %in% cell_groups) %>%
        dplyr::mutate(
            .gene_label = as.character(.data[[gene_label_col]]),
            .cell_group = as.character(.data[[cell_group_col]]),
            .lfc = .data[[lfc_col]],
            .is_sig = .data[[p_value_col]] < sig_p_value_thresh
        )

    mat_df <- plot_df %>%
        dplyr::select(.gene_label, .cell_group, .lfc) %>%
        tidyr::pivot_wider(
            names_from = .cell_group,
            values_from = .lfc,
            values_fill = 0
        )

    mat <- as.matrix(mat_df[, setdiff(names(mat_df), ".gene_label"), drop = FALSE])
    rownames(mat) <- mat_df$.gene_label
    mat[is.na(mat)] <- 0

    row_order <- if (nrow(mat) > 1) {
        stats::hclust(stats::dist(mat))$order
    } else {
        1L
    }
    col_order <- if (ncol(mat) > 1) {
        stats::hclust(stats::dist(t(mat)))$order
    } else {
        1L
    }

    gene_levels <- rownames(mat)[row_order]
    cell_group_levels <- colnames(mat)[col_order]

    plot_df <- plot_df %>%
        dplyr::mutate(
            .gene_label = factor(.gene_label, levels = gene_levels),
            .cell_group = factor(.cell_group, levels = cell_group_levels)
        )

    plot_df %>%
        ggplot2::ggplot(ggplot2::aes(.cell_group, .gene_label, fill = .lfc)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#006600", mid = "white", high = "#800080") +
        ggplot2::geom_tile(
            data = ~ dplyr::filter(.x, .is_sig),
            fill = NA,
            color = "black",
            linewidth = 0.6
        ) +
        ggplot2::coord_equal() +
        monocle3:::monocle_theme_opts() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::ylab("gene") +
        ggplot2::xlab("cell type") +
        ggplot2::labs(fill = "log(FC)")
}
