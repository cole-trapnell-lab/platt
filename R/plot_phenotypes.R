# Define the base rainbow timepoint colors
rainbow_timepoint_colors <-
    c(
        "18h" = "#DF4828",
        "24h" = "#E78C35",
        "36h" = "#F6C141",
        "48h" = "#4EB265",
        "72h" = "#1965B0"
    )

phenotype_colors <- c(
    "none"           = "#8f8f8f",
    "abundance_gain" = "#e41a1c",
    "abundance_loss" = "#377eb8",
    "identity"       = "#ff8c00",
    "fitness"        = "#ff7f00",
    "apoptosis"      = "#000000",
    "stress"         = "#4daf4a",
    "senescence"     = "#a65628"
)

#' Plot Phenotype Glyphs from an Impact Table
#'
#' Convert a standardized impact table into phenotype glyph fields and render
#' an annotated lineage graph.
#'
#' @param cell_state_graph A cell-state graph object containing node layout in
#'   `@g` and edge/layout metadata in `@layout_info`.
#' @param impact_table A tibble/data frame with one row per cell type and
#'   phenotype annotation columns (see `impact_to_phenos()`).
#' @param filter_by_group Logical. If `TRUE` and `cell_types` is provided, only
#'   keep nodes in the same `group_nodes_by` group as `cell_types`.
#' @param cell_types Optional character vector of cell types to retain.
#' @param show_node_labels Logical. If `TRUE`, label all nodes unless
#'   `label_cell_types` is set through `...`.
#' @param label_font_size Numeric label size used by `ggrepel` (mm units).
#' @param label_font_size_pt Optional numeric label size in points. If set, this
#'   overrides `label_font_size`.
#' @param show_group_labels Logical. If `TRUE`, draw text labels for grouping
#'   boxes.
#' @param group_label_font_size Numeric size for grouping-box labels.
#' @param ... Additional arguments passed to `plot_phenotypes_glyphs()`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' p <- plot_phenotypes_from_impact(
#'   cell_state_graph = cell_state_graph,
#'   impact_table = impact_table,
#'   show_node_labels = TRUE,
#'   label_font_size_pt = 14,
#'   show_group_labels = TRUE
#' )
#' }
#' @export
plot_phenotypes_from_impact <- function(cell_state_graph,
                                        impact_table,
                                        filter_by_group = FALSE,
                                        cell_types = NULL,
                                        show_node_labels = FALSE,
                                        label_font_size = 3,
                                        label_font_size_pt = NULL,
                                        show_group_labels = FALSE,
                                        group_label_font_size = 2,
                                        ... # passthrough to plot_phenotypes_glyphs
) {
    phenos <- impact_to_phenos(impact_table)
    plot_phenotypes_glyphs(
        cell_state_graph,
        phenos_df = phenos,
        filter_by_group = filter_by_group,
        cell_types = cell_types,
        show_node_labels = show_node_labels,
        label_font_size = label_font_size,
        label_font_size_pt = label_font_size_pt,
        show_group_labels = show_group_labels,
        group_label_font_size = group_label_font_size,
        ...
    )
}

#' Convert Impact Table Columns into Phenotype Glyph Fields
#'
#' Translates impact-table phenotype calls (abundance, identity, fitness) into a
#' standardized phenotype table used by `plot_phenotypes_glyphs()`.
#'
#' @param impact_table A data frame with required columns:
#'   `cell_type`, `abundance_code`, `abundance_severity`, `identity_label`,
#'   `fitness_label`. Optional columns include `fitness_evidence`,
#'   `identity_evidence`, and nested `disrupted_pathways` (or legacy
#'   `llm_disrupted_pathways`).
#' @param sev_map Named numeric vector mapping abundance severity to magnitude.
#' @param abundance_code_map Named numeric vector mapping abundance class labels
#'   to signed effect direction/magnitude.
#' @param identity_glyph_map Named character vector mapping identity labels to
#'   glyph symbols.
#'
#' @return A tibble with columns expected by `plot_phenotypes_glyphs()`,
#'   including `cell_group`, `abundance_log2fc`, identity glyph fields, fitness
#'   axes, and optional `dysregulated_genes`.
#' @export
impact_to_phenos <- function(impact_table,
                             sev_map = c(none = 0.25, mild = 0.6, moderate = 1.2, severe = 2.0),
                             abundance_code_map = c(
                                 "A0 No change" = 0,
                                 "A1 Expansion" = +1,
                                 "A2 Depletion" = -1,
                                 "A3 Ablation/Loss" = -1.5,
                                 "A4 Ectopic/extra state" = +1.2
                             ),
                             identity_glyph_map = c(
                                 "I0 Identity intact" = "",
                                 "I1 Maturation delay" = "<<",
                                 "I2 Precocious maturation" = ">>",
                                 "I3 Program failure within identity" = "!!",
                                 "I4 Fate switch / misspecification" = "⇄",
                                 "I5 Identity fragmentation" = ""
                             )) {
    stopifnot(all(c(
        "cell_type", "abundance_code", "abundance_severity",
        "identity_label", "fitness_label"
    ) %in% names(impact_table)))
    tab <- impact_table

    # Normalize text
    norm <- function(x) trimws(as.character(x))
    tab$abundance_code <- norm(tab$abundance_code)
    tab$abundance_severity <- tolower(norm(tab$abundance_severity))
    tab$identity_label <- norm(tab$identity_label)
    tab$fitness_label <- norm(tab$fitness_label)

    # Derive an abundance *proxy* log2FC for coloring (still fixed node size).
    code_sign <- setNames(unname(abundance_code_map), names(abundance_code_map))
    sev_mag <- setNames(unname(sev_map), names(sev_map))
    tab$._abund_sign <- unname(code_sign[tab$abundance_code])
    tab$._sev_mag <- unname(sev_mag[tab$abundance_severity])
    tab$._sev_mag[is.na(tab$._sev_mag)] <- 0
    tab$abundance_log2fc_proxy <- ifelse(is.na(tab$._abund_sign), 0, tab$._abund_sign * tab$._sev_mag)

    # Parse identity to a compact label (I0–I5) if full phrase given
    id_raw <- tab$identity_label
    pick_I <- function(s) {
        if (is.na(s) || s == "") {
            return("I0")
        }
        if (grepl("^I[0-5]", s)) {
            return(sub("^((I[0-5])).*$", "\\1", s))
        }
        if (grepl("delay", s, ignore.case = TRUE)) {
            return("I1")
        }
        if (grepl("precoci", s, ignore.case = TRUE)) {
            return("I2")
        }
        if (grepl("fate switch|misspec", s, ignore.case = TRUE)) {
            return("I4")
        }
        if (grepl("program failure|effector|module", s, ignore.case = TRUE)) {
            return("I3")
        }
        if (grepl("fragment", s, ignore.case = TRUE)) {
            return("I5")
        }
        "I0"
    }
    ident_code <- vapply(id_raw, pick_I, character(1))
    tab$identity_code <- ident_code
    id_map_full <- c(
        setNames(identity_glyph_map, names(identity_glyph_map)),
        c(I0 = "", I1 = "<<", I2 = ">>", I3 = "!!", I4 = "⇄", I5 = "")
    )
    tab$identity_glyph <- unname(id_map_full[ifelse(grepl("^I[0-5]$", ident_code), ident_code, tab$identity_label)])
    tab$identity_glyph[is.na(tab$identity_glyph)] <- ""

    # Parse fitness_label into axes
    fraw <- ifelse(is.na(tab$fitness_label), "", tab$fitness_label)
    tab$F1_dir <- dplyr::case_when(
        grepl("F1", fraw) & grepl("increase|up|\\b\\+\\b", fraw, ignore.case = TRUE) ~ "increase",
        grepl("F1", fraw) & grepl("decrease|down|\\b-\\b", fraw, ignore.case = TRUE) ~ "decrease",
        TRUE ~ NA_character_
    )
    tab$F2_apoptosis <- grepl("F2|apoptos", fraw, ignore.case = TRUE)
    if ("fitness_evidence" %in% names(tab)) {
        ev <- ifelse(is.na(tab$fitness_evidence), "", tab$fitness_evidence)
        get_max_nes <- function(s) {
            # Extract all NES=... values and return the max
            matches <- stringr::str_match_all(s, "NES=([-+]?[0-9]*\\.?[0-9]+)")
            nes_vals <- as.numeric(unlist(lapply(matches, function(m) m[, 2])))
            if (length(nes_vals) && !all(is.na(nes_vals))) max(nes_vals, na.rm = TRUE) else NA_real_
        }
        tab$F3_stress_score <- ifelse(
            grepl("F3|stress|p53|UPR|ROS|hypoxia|interferon", fraw, ignore.case = TRUE),
            vapply(ev, get_max_nes, numeric(1)),
            NA_real_
        )
    } else {
        tab$F3_stress_score <- ifelse(grepl("F3|stress|p53|UPR|ROS|hypoxia|interferon", fraw, ignore.case = TRUE), 1, NA_real_)
    }
    # Parse stress evidence for F3 rows
    parse_stress_evidence <- function(evidence_str) {
        if (is.na(evidence_str) || evidence_str == "") {
            return(NA_character_)
        }
        # Split by semicolon, trim whitespace
        pieces <- strsplit(evidence_str, ";")[[1]]
        pieces <- trimws(pieces)
        # Only keep non-empty
        pieces <- pieces[nzchar(pieces)]
        if (length(pieces) == 0) {
            return(NA_character_)
        }
        # Collapse with <br> for HTML tooltip
        paste(pieces, collapse = "<br>")
    }
    tab$stress_evidence <- ifelse(
        !is.na(tab$F3_stress_score) & tab$F3_stress_score != 0 & "fitness_evidence" %in% names(tab),
        vapply(tab$fitness_evidence, parse_stress_evidence, character(1)),
        NA_character_
    )

    tab$F4_senescence <- grepl("F4|senesc|quiesc", fraw, ignore.case = TRUE)

    # --- Prepare dysregulated genes per cell type ---
    collapse_unique_genes <- function(x) {
        x <- as.character(x)
        x <- stringr::str_squish(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (!length(x)) {
            return(NA_character_)
        }

        # Split compound entries (e.g. "geneA, geneB; geneC") into individual tokens.
        # Then de-duplicate by a normalized key while preserving first-seen display text.
        tokens <- unlist(strsplit(x, "[,;]"), use.names = FALSE)
        tokens <- stringr::str_squish(tokens)
        tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
        if (!length(tokens)) {
            return(NA_character_)
        }

        key <- tokens %>%
            stringr::str_replace_all("[\u2191\u2193]", "") %>%
            stringr::str_replace_all("\\(.*?\\)", "") %>%
            stringr::str_squish() %>%
            tolower()
        tokens_unique <- tokens[!duplicated(key)]
        paste(tokens_unique, collapse = ", ")
    }

    pathway_col <- if ("disrupted_pathways" %in% names(tab)) {
        "disrupted_pathways"
    } else if ("llm_disrupted_pathways" %in% names(tab)) {
        "llm_disrupted_pathways"
    } else {
        NULL
    }

    if (!is.null(pathway_col)) {
        gene_tbl <- tab %>%
            tidyr::unnest(dplyr::all_of(pathway_col)) %>%
            dplyr::select(cell_group = cell_type, dysregulated_genes) %>%
            tidyr::unnest_longer(dysregulated_genes) %>%
            dplyr::mutate(dysregulated_genes = stringr::str_squish(as.character(dysregulated_genes))) %>%
            dplyr::filter(!is.na(dysregulated_genes), nzchar(dysregulated_genes)) %>%
            dplyr::group_by(cell_group) %>%
            dplyr::summarise(
                dysregulated_genes = collapse_unique_genes(dysregulated_genes),
                .groups = "drop"
            )
    } else {
        gene_tbl <- tibble::tibble(cell_group = character(), dysregulated_genes = character())
    }

    phenos <- dplyr::tibble(
        cell_group         = tab$cell_type,
        abundance_code     = tab$abundance_code,
        abundance_severity = tab$abundance_severity,
        abundance_log2fc   = tab$abundance_log2fc_proxy,
        abundance_q        = NA_real_,
        identity_label     = tab$identity_code,
        identity_glyph     = tab$identity_glyph,
        F1_dir             = tab$F1_dir,
        F2_apoptosis       = tab$F2_apoptosis,
        F3_stress_score    = tab$F3_stress_score,
        F4_senescence      = tab$F4_senescence,
        stress_evidence    = tab$stress_evidence,
        identity_evidence  = if ("identity_evidence" %in% names(tab)) as.character(tab$identity_evidence) else NA_character_,
        expectation        = if ("expectation" %in% names(tab)) as.character(tab$expectation) else NA_character_,
        rationale          = if ("rationale" %in% names(tab)) as.character(tab$rationale) else NA_character_
    )

    # Join dysregulated genes if available
    phenos <- phenos %>%
        dplyr::left_join(gene_tbl, by = "cell_group")

    return(phenos)
}



#' Plot Lineage Graph with Phenotype Overlays
#'
#' Draws a lineage graph with phenotype-driven node fills, optional identity
#' glyphs/fitness badges, and optional node/group labels.
#'
#' @param cell_state_graph A cell-state graph object with node coordinates
#'   (`@g`) and layout metadata (`@layout_info`).
#' @param phenos_df Phenotype table, typically from `impact_to_phenos()`.
#' @param map Named list mapping expected roles to `phenos_df` columns:
#'   `id`, `lfc`, `q`, `ident`, `glyph`, `f1`, `f2`, `f3`, `f4`.
#' @param lfc_cap Numeric cap applied to abundance log2FC proxy.
#' @param arrow_unit Numeric arrowhead size for edges (points).
#' @param node_size Base node size scalar.
#' @param con_colour Edge/node outline color.
#' @param legend_position Legend position passed to `theme(legend.position=...)`.
#'   Default is `"none"` to preserve historical no-legend behavior.
#'   (e.g. `"right"`, `"bottom"`, `"top"`, `"left"`, or `"none"` to hide).
#' @param label_cell_types Optional character vector of node names to label, or
#'   `"all"`.
#' @param show_node_labels Logical. If `TRUE` and `label_cell_types` is `NULL`,
#'   labels all nodes.
#' @param label_font_size Numeric node label size in mm units (`ggrepel` scale).
#' @param label_font_size_pt Optional node label size in points; overrides
#'   `label_font_size` when provided.
#' @param stress_cap Numeric cap for stress-score alpha scaling.
#' @param filter_by_group Logical. If `TRUE` and `cell_types` is provided,
#'   restrict to the shared grouping-box group.
#' @param cell_types Optional character vector of node names to retain.
#' @param draw_group_boxes Logical. Draw grouping boxes from layout metadata.
#' @param show_group_labels Logical. Draw text labels above grouping boxes.
#' @param group_label_font_size Numeric grouping-box label size.
#' @param node_overlay One of `"none"`, `"glyphs"`, `"badges"`, `"both"`.
#' @param glyph_color Optional fixed color for identity glyph text.
#' @param badge_color Fill/stroke color for fitness badges.
#' @param badge_outline_color Outline color for filled badge shapes.
#' @param render_mode One of `"tissue"` or `"global"`. Global mode hides
#'   glyph/badge overlays and can draw a center marker for identity phenotypes.
#' @param global_identity_marker Logical; when `TRUE`, draw a small black center
#'   marker for nodes with non-intact identity (used for global view).
#'
#' @return A `ggplot` object.
#' @export
plot_phenotypes_glyphs <- function(cell_state_graph,
                                   phenos_df,
                                   map = list(
                                       id = "cell_group",
                                       lfc = "abundance_log2fc",
                                       q = "abundance_q",
                                       ident = "identity_label",
                                       glyph = "identity_glyph",
                                       f1 = "F1_dir",
                                       f2 = "F2_apoptosis",
                                       f3 = "F3_stress_score",
                                       f4 = "F4_senescence"
                                   ),
                                   lfc_cap = 2,
                                   arrow_unit = 3,
                                   node_size = 2.2,
                                   con_colour = "darkgrey",
                                   legend_position = "none",
                                   label_cell_types = NULL,
                                   show_node_labels = FALSE,
                                   label_font_size = 3,
                                   label_font_size_pt = NULL,
                                   stress_cap = 2.5,
                                   filter_by_group = FALSE,
                                   cell_types = NULL,
                                   draw_group_boxes = TRUE,
                                   show_group_labels = FALSE,
                                   group_label_font_size = 2,
                                   node_overlay = c("none", "glyphs", "badges", "both"),
                                   glyph_color = NULL,
                                   badge_color = "black",
                                   badge_outline_color = "white",
                                   render_mode = c("tissue", "global"),
                                   global_identity_marker = NULL) {
    node_overlay <- match.arg(node_overlay)
    render_mode <- match.arg(render_mode)
    show_node_glyphs <- node_overlay %in% c("glyphs", "both")
    show_node_badges <- node_overlay %in% c("badges", "both")
    if (is.null(global_identity_marker)) {
        global_identity_marker <- identical(render_mode, "global")
    }
    if (identical(render_mode, "global")) {
        show_node_glyphs <- FALSE
        show_node_badges <- FALSE
    }

    g <- cell_state_graph@g %>% dplyr::mutate(name = stringr::str_trim(name))
    bezier_df <- cell_state_graph@layout_info$bezier_df
    grouping_df <- cell_state_graph@layout_info$grouping_df

    # --- Filter by group if requested ---
    if (filter_by_group && !is.null(cell_types)) {
        if (is.null(grouping_df)) stop("grouping_df is not available in the cell_state_graph.")
        group_val <- grouping_df %>%
            dplyr::filter(id %in% cell_types) %>%
            dplyr::pull(group_nodes_by) %>%
            unique()
        if (length(group_val) != 1) stop("cell_types must all belong to a single group.")
        g <- g %>% dplyr::filter(group_nodes_by == group_val)
        phenos_df <- phenos_df %>% dplyr::filter(.data[[map$id]] %in% g$name)
        bezier_df <- bezier_df %>% dplyr::filter(from %in% g$name & to %in% g$name)
    } else if (!is.null(cell_types)) {
        g <- g %>% dplyr::filter(name %in% cell_types)
        phenos_df <- phenos_df %>% dplyr::filter(.data[[map$id]] %in% cell_types)
        bezier_df <- bezier_df %>% dplyr::filter(from %in% g$name & to %in% g$name)
    }

    pdf <- phenos_df %>%
        dplyr::transmute(
            name = .data[[map$id]],
            lfc = as.numeric(.data[[map$lfc]]),
            q = if (map$q %in% names(phenos_df)) as.numeric(.data[[map$q]]) else NA_real_,
            abundance_code = if ("abundance_code" %in% names(phenos_df)) as.character(.data[["abundance_code"]]) else NA_character_,
            ident = if (map$ident %in% names(phenos_df)) as.character(.data[[map$ident]]) else "I0",
            glyph = if (map$glyph %in% names(phenos_df)) as.character(.data[[map$glyph]]) else "",
            f1_dir = if (map$f1 %in% names(phenos_df)) as.character(.data[[map$f1]]) else NA_character_,
            fitness_label = if ("fitness_label" %in% names(phenos_df)) as.character(.data[["fitness_label"]]) else NA_character_,
            f2 = if (map$f2 %in% names(phenos_df)) as.logical(.data[[map$f2]]) else NA,
            f3 = if (map$f3 %in% names(phenos_df)) as.numeric(.data[[map$f3]]) else NA_real_,
            f4 = if (map$f4 %in% names(phenos_df)) as.logical(.data[[map$f4]]) else NA,
            dysregulated_genes = if ("dysregulated_genes" %in% names(phenos_df)) as.character(.data[["dysregulated_genes"]]) else NA_character_,
            identity_evidence = if ("identity_evidence" %in% names(phenos_df)) as.character(.data[["identity_evidence"]]) else NA_character_,
            stress_evidence = if ("stress_evidence" %in% names(phenos_df)) as.character(.data[["stress_evidence"]]) else NA_character_,
            expectation = if ("expectation" %in% names(phenos_df)) as.character(.data[["expectation"]]) else NA_character_,
            rationale = if ("rationale" %in% names(phenos_df)) as.character(.data[["rationale"]]) else NA_character_
        )

    g <- g %>%
        dplyr::left_join(pdf, by = "name")

    # Helper to show phenotype only if not NA or "no phenotype" code
    show_pheno <- function(val, none_codes) {
        ifelse(is.na(val) | val %in% none_codes, "", val)
    }

    g <- g %>%
        dplyr::mutate(
            identity_display = dplyr::case_when(
                is.na(ident) ~ "",
                ident == "I0" ~ "Identity intact",
                ident == "I1" ~ "Maturation delay",
                ident == "I2" ~ "Precocious maturation",
                ident == "I3" ~ "Program failure within identity",
                ident == "I4" ~ "Fate switch / misspecification",
                ident == "I5" ~ "Identity fragmentation",
                TRUE ~ ident
            ),
            abundance_str = show_pheno(abundance_code, c("A0 No change", NA)),
            identity_str = show_pheno(identity_display, c("I0", "I0 Identity intact", "Identity intact", NA)),
            fitness_str = show_pheno(f1_dir, c("F0 No significant phenotype", "F0 Normal", NA)),
            badge_caption = paste0(
                ifelse(!is.na(f1_dir) & f1_dir == "increase", "F1_up, ", ""),
                ifelse(!is.na(f1_dir) & f1_dir == "decrease", "F1_down, ", ""),
                ifelse(!is.na(f2) & f2, "F2, ", ""),
                ifelse(!is.na(f3) & f3 != 0, "F3, ", ""),
                ifelse(!is.na(f4) & f4, "F4, ", "")
            ),
            badge_caption = stringr::str_replace(badge_caption, ",\\s*$", ""),
            tooltip = paste0(
                "<b>", name, "</b><br>",
                ifelse(identity_str != "", paste0(ifelse(glyph != "", paste0(glyph, " "), ""), "Major transcriptional phenotype: ", identity_str, "<br>"), ""),
                ifelse(abundance_str != "", paste0("Abundance: ", abundance_str, "<br>"), ""),
                ifelse(!identical(render_mode, "tissue") & !is.na(lfc), paste0("Abundance logFC: ", formatC(lfc, digits = 2, format = "f"), "<br>"), ""),
                ifelse(!is.na(q), paste0("Abundance logFC q-value: ", formatC(q, digits = 2, format = "e"), "<br>"), ""),
                ifelse(!is.na(expectation) & expectation != "", paste0("Expected: ", expectation, "<br>"), ""),
                ifelse(!is.na(rationale) & rationale != "", paste0("Expectation rationale: ", stringr::str_replace_all(stringr::str_wrap(stringr::str_trunc(rationale, 300), width = 50), "\n", "<br>"), "<br>"), ""),
                ifelse(
                    identity_str != "" & !is.na(identity_evidence) & identity_evidence != "",
                    paste0(
                        "Identity evidence: ",
                        stringr::str_replace_all(
                            stringr::str_wrap(stringr::str_trunc(identity_evidence, 300), width = 50),
                            "\n", "<br>"
                        ),
                        "<br>"
                    ),
                    ""
                ),
                ifelse(!is.na(f1_dir) & f1_dir == "increase", "▲ Proliferation: increase<br>", ""),
                ifelse(!is.na(f1_dir) & f1_dir == "decrease", "▼ Proliferation: decrease<br>", ""),
                ifelse(!is.na(f2) & f2, "─ Apoptosis: yes<br>", ""),
                ifelse(!is.na(f3) & f3 != 0, "* Stress response<br>", ""),
                ifelse(!is.na(stress_evidence) & stress_evidence != "", paste0("Stress evidence:<br>", stress_evidence, "<br>"), ""),
                ifelse(!is.na(f4) & f4, "□ Senescence: yes<br>", ""),
                ifelse(
                    !is.na(dysregulated_genes) & dysregulated_genes != "",
                    paste0(
                        "Dysregulated genes: ",
                        # First truncate and wrap, then replace arrows with HTML
                        stringr::str_trunc(dysregulated_genes, 300) %>%
                            stringr::str_wrap(width = 50) %>%
                            stringr::str_replace_all("\n", "<br>") %>%
                            stringr::str_replace_all("↑", "<span style='color:red;'>↑</span>") %>%
                            stringr::str_replace_all("↓", "<span style='color:blue;'>↓</span>"),
                        "<br>"
                    ),
                    ""
                )
            ),
            lfc_capped = pmin(pmax(lfc, -lfc_cap), lfc_cap),
            glyph = dplyr::coalesce(glyph, ""),
            glyph_col = ifelse(is.na(lfc), "black", ifelse(abs(lfc) >= 0.8, "white", "black")),
            f3_alpha = dplyr::case_when(
                is.na(f3) ~ 0,
                TRUE ~ scales::rescale(pmin(abs(f3), stress_cap), to = c(0.25, 1))
            ),
            node_size_plot = node_size * 3.2
        )

    color_priority <- c("abundance_loss", "identity", "abundance_gain", "none")
    g <- g %>%
        dplyr::mutate(
            primary_phenotype = dplyr::case_when(
                !is.na(abundance_code) & grepl("^A2|^A3", abundance_code) ~ "abundance_loss",
                !is.na(ident) & !(ident %in% c("I0", "I0 Identity intact", "Identity intact", "", NA)) ~ "identity",
                !is.na(abundance_code) & grepl("^A1|^A4", abundance_code) ~ "abundance_gain",
                TRUE ~ "none"
            ),
            has_identity_change = !is.na(ident) & !(ident %in% c("I0", "I0 Identity intact", "Identity intact", "", NA)),
            is_expected = !is.na(expectation) & tolower(expectation) %in% c("expected", "expected change")
        )

    p <- ggplot2::ggplot(ggplot2::aes(x, y), data = g) +
        ggplot2::geom_path(
            ggplot2::aes(x, y, group = edge_name),
            linewidth = 0.25, colour = con_colour, data = dplyr::distinct(bezier_df),
            arrow = ggplot2::arrow(angle = 30, length = grid::unit(arrow_unit, "pt"), type = "closed"),
            linejoin = "mitre"
        ) +
        ggnetwork::theme_blank() +
        hooke_theme_opts()

    # Only draw group boxes if requested
    if (draw_group_boxes && !is.null(grouping_df) && !identical(grouping_df$group_nodes_by, grouping_df$id)) {
        p <- p + ggforce::geom_mark_rect(
            ggplot2::aes(x, y, group = group_nodes_by, color = I("lightgrey")),
            size = 0.25, radius = grid::unit(0.5, "mm"),
            expand = grid::unit(1, "mm"),
            con.type = "straight", con.colour = "lightgrey",
            con.size = 0.25, con.border = "one", na.rm = TRUE, data = g
        )
    }
    if (show_group_labels && !is.null(grouping_df) && !identical(grouping_df$group_nodes_by, grouping_df$id)) {
        yrange <- max(g$y, na.rm = TRUE) - min(g$y, na.rm = TRUE)
        group_labels <- g %>%
            dplyr::select(x, y, group_nodes_by) %>%
            dplyr::distinct() %>%
            dplyr::group_by(group_nodes_by) %>%
            dplyr::summarise(
                x = mean(x, na.rm = TRUE),
                y = max(y, na.rm = TRUE) + 0.02 * yrange,
                .groups = "drop"
            )
        p <- p + ggplot2::geom_text(
            data = group_labels,
            ggplot2::aes(x = x, y = y, label = group_nodes_by),
            size = group_label_font_size,
            color = "black"
        )
    }

    # Single node fill with explicit priority: loss > identity > gain > none.
    p <- p +
        ggiraph::geom_point_interactive(
            ggplot2::aes(x = x, y = y, tooltip = tooltip, fill = primary_phenotype, size = node_size_plot),
            data = g,
            shape = 21,
            color = if (identical(render_mode, "global")) NA else con_colour,
            linewidth = if (identical(render_mode, "global")) 0 else 0.25
        ) +
        ggplot2::scale_size_identity() +
        ggplot2::scale_fill_manual(
            values = phenotype_colors,
            breaks = color_priority,
            labels = c(
                abundance_loss = "Abundance decrease",
                identity = "Major transcriptional phenotype",
                abundance_gain = "Abundance increase",
                none = "No phenotype"
            ),
            name = if (identical(render_mode, "global")) NULL else "Node color",
            guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4, alpha = 1))
        )

    expected_nodes <- g %>% dplyr::filter(is_expected)
    if (nrow(expected_nodes) > 0) {
        expected_nodes$expected_outline <- "Expected (black outline)"
        p <- p +
            ggplot2::geom_point(
                data = expected_nodes,
                ggplot2::aes(x = x, y = y, color = expected_outline, size = node_size_plot),
                shape = 21,
                fill = NA,
                stroke = if (identical(render_mode, "global")) 0.035 else 0.35,
                show.legend = TRUE
            ) +
            ggplot2::scale_color_manual(
                name = "Markers",
                values = c("Expected (black outline)" = "black")
            )
    }

    p <- p +
        ggplot2::coord_equal(clip = "off") +
        ggplot2::theme(
            legend.position = legend_position,
            plot.margin = ggplot2::margin(10, 10, 10, 10)
        )

    if (isTRUE(global_identity_marker)) {
        p <- p +
            ggiraph::geom_point_interactive(
                data = g %>% dplyr::filter(has_identity_change),
                ggplot2::aes(x = x, y = y, tooltip = tooltip),
                shape = 16,
                size = node_size * 0.275,
                color = "black"
            )
    }

    if (show_node_glyphs) {
        glyph_draw_color <- if (is.null(glyph_color)) g$glyph_col else glyph_color
        p <- p +
            ggiraph::geom_text_interactive(
                ggplot2::aes(x = x, y = y, label = glyph, tooltip = tooltip),
                data = g, size = node_size * 1.4, color = glyph_draw_color, fontface = "bold", vjust = 0.35
            )
    }

    if (show_node_badges) {
        mk_badges <- function(df) {
            offs <- tibble::tibble(
                badge = c("F1_up", "F1_down", "F2", "F3", "F4"),
                dx = c(0.00, 0.00, -0.23, 0.23, 0.00),
                dy = c(0.27, -0.27, 0.00, 0.00, 0.36)
            )
            base <- df %>% dplyr::select(dplyr::any_of(c("name", "x", "y", "node_size_plot", "node_size", "f1_dir", "f2", "f3", "f4", "f3_alpha")))
            if (!("node_size_plot" %in% names(base))) {
                if ("node_size" %in% names(base)) {
                    base <- base %>% dplyr::mutate(node_size_plot = node_size * 3.2)
                } else {
                    base <- base %>% dplyr::mutate(node_size_plot = 3.2)
                }
            }
            b1u <- base %>%
                dplyr::filter(f1_dir == "increase") %>%
                dplyr::mutate(badge = "F1_up")
            b1d <- base %>%
                dplyr::filter(f1_dir == "decrease") %>%
                dplyr::mutate(badge = "F1_down")
            b2 <- base %>%
                dplyr::filter(isTRUE(f2)) %>%
                dplyr::mutate(badge = "F2")
            b3 <- base %>%
                dplyr::filter(!is.na(f3) & f3_alpha > 0) %>%
                dplyr::mutate(badge = "F3")
            b4 <- base %>%
                dplyr::filter(isTRUE(f4)) %>%
                dplyr::mutate(badge = "F4")
            dplyr::bind_rows(b1u, b1d, b2, b3, b4) %>%
                dplyr::left_join(offs, by = "badge") %>%
                dplyr::mutate(
                    # Keep badge placement proportional to the actual rendered node size.
                    bx = x + dx * node_size_plot,
                    by = y + dy * node_size_plot
                )
        }

        badges <- mk_badges(g)
        if (nrow(badges)) {
            p <- p +
                ggnewscale::new_scale_color() +
                ggplot2::geom_point(
                    data = badges %>% dplyr::filter(badge == "F1_up"),
                    ggplot2::aes(bx, by), shape = 24, size = 1.9, fill = badge_color, color = badge_outline_color, stroke = 0.25
                ) +
                ggplot2::geom_point(
                    data = badges %>% dplyr::filter(badge == "F1_down"),
                    ggplot2::aes(bx, by), shape = 25, size = 1.9, fill = badge_color, color = badge_outline_color, stroke = 0.25
                ) +
                ggplot2::geom_point(
                    data = badges %>% dplyr::filter(badge == "F2"),
                    ggplot2::aes(bx, by), shape = 95, size = 1.9, color = badge_color, stroke = 0.6
                ) +
                ggplot2::geom_text(
                    data = badges %>% dplyr::filter(badge == "F3"),
                    ggplot2::aes(x = bx, y = by, alpha = pmin(1, f3_alpha)),
                    label = "*", size = 4.6, color = badge_color, fontface = "bold"
                ) +
                ggplot2::geom_point(
                    data = badges %>% dplyr::filter(badge == "F4"),
                    ggplot2::aes(bx, by),
                    shape = 22, size = 1.9, fill = badge_color, color = badge_outline_color, stroke = 0.25
                ) +
                ggplot2::scale_alpha_continuous(range = c(0.25, 1), guide = "none")
        }
    }

    label_target <- label_cell_types
    if (isTRUE(show_node_labels) && is.null(label_target)) {
        label_target <- "all"
    }
    if (!is.null(label_target)) {
        lbl <- if (identical(label_target, "all")) g else g %>% dplyr::filter(name %in% label_target)
        lbl <- lbl %>% dplyr::distinct(name, x, y)
        label_size_mm <- if (is.null(label_font_size_pt)) label_font_size else label_font_size_pt / ggplot2::.pt
        p <- p + ggrepel::geom_text_repel(
            data = lbl, ggplot2::aes(x, y, label = name),
            size = label_size_mm, color = "black",
            box.padding = 0.3, point.padding = 0.3, segment.color = "grey50"
        )
    }

    p
}
