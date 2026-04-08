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
    "none"           = "#ffffff",
    "abundance_gain" = "#f6c141",
    "abundance_loss" = "#e41a1c",
    "identity"       = "#ff8c00",
    "fitness"        = "#f6c141",
    "apoptosis"      = "#000000",
    "stress"         = "#4daf4a",
    "senescence"     = "#a65628",
    "not present"    = "#cccccc"
)


format_group_display_name <- function(x) {
    x <- trimws(as.character(x))
    x <- dplyr::na_if(x, "")
    out <- stringr::str_to_title(stringr::str_replace_all(x, "_", " "))
    out[out == "Cns Other"] <- "CNS Other"
    out
}

format_abundance_change_percent <- function(lfc) {
    pct <- (exp(lfc) - 1) * 100
    dplyr::if_else(
        is.na(lfc),
        NA_character_,
        paste0(formatC(pct, digits = 0, format = "f", flag = "+"), "%")
    )
}

abundance_power_status_label <- function(power, powered_thresh = 0.8) {
    dplyr::case_when(
        is.na(power) ~ NA_character_,
        power >= powered_thresh ~ "powered",
        TRUE ~ "underpowered"
    )
}

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
#'     cell_state_graph = cell_state_graph,
#'     impact_table = impact_table,
#'     show_node_labels = TRUE,
#'     label_font_size_pt = 14,
#'     show_group_labels = TRUE
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
                                 "I4 Fate switch / misspecification" = "⬌",
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
        c(I0 = "", I1 = "<<", I2 = ">>", I3 = "!!", I4 = "⬌", I5 = "")
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
        pathway_long <- tab %>%
            tidyr::unnest(dplyr::all_of(pathway_col))
        gene_tbl <- pathway_long %>%
            dplyr::select(cell_group = cell_type, dysregulated_genes) %>%
            tidyr::unnest_longer(dysregulated_genes) %>%
            dplyr::mutate(dysregulated_genes = stringr::str_squish(as.character(dysregulated_genes))) %>%
            dplyr::filter(!is.na(dysregulated_genes), nzchar(dysregulated_genes)) %>%
            dplyr::group_by(cell_group) %>%
            dplyr::summarise(
                dysregulated_genes = collapse_unique_genes(dysregulated_genes),
                .groups = "drop"
            )
        pathway_tbl <- pathway_long %>%
            dplyr::transmute(
                cell_group = cell_type,
                pathway_name = as.character(name)
            ) %>%
            dplyr::filter(!is.na(pathway_name), nzchar(stringr::str_squish(pathway_name))) %>%
            dplyr::mutate(pathway_name = stringr::str_squish(pathway_name)) %>%
            dplyr::group_by(cell_group) %>%
            dplyr::summarise(
                pathway_names = paste(unique(pathway_name), collapse = "; "),
                .groups = "drop"
            )
        pathway_gene_tbl <- pathway_long %>%
            dplyr::transmute(
                cell_group = as.character(cell_type),
                pathway_name = as.character(name),
                dysregulated_genes = if ("dysregulated_genes" %in% names(pathway_long)) dysregulated_genes else NA
            ) %>%
            dplyr::mutate(
                pathway_name = stringr::str_squish(pathway_name),
                pathway_name = ifelse(is.na(pathway_name) | !nzchar(pathway_name), "Unknown pathway", pathway_name)
            )
        if (is.list(pathway_gene_tbl$dysregulated_genes)) {
            pathway_gene_tbl <- pathway_gene_tbl %>%
                dplyr::mutate(genes_text = purrr::map_chr(dysregulated_genes, ~ collapse_unique_genes(as.character(.x))))
        } else {
            pathway_gene_tbl <- pathway_gene_tbl %>%
                dplyr::mutate(genes_text = purrr::map_chr(dysregulated_genes, ~ collapse_unique_genes(as.character(.x))))
        }
        pathway_gene_tbl <- pathway_gene_tbl %>%
            dplyr::group_by(cell_group, pathway_name) %>%
            dplyr::summarise(
                genes_text = collapse_unique_genes(genes_text),
                .groups = "drop"
            ) %>%
            dplyr::group_by(cell_group) %>%
            dplyr::arrange(pathway_name, .by_group = TRUE) %>%
            dplyr::mutate(
                pathway_line = paste0(
                    "Pathway: ", pathway_name, "<br>&nbsp;&nbsp;Genes dysregulated: ",
                    ifelse(is.na(genes_text) | !nzchar(genes_text), "NA", genes_text)
                )
            ) %>%
            dplyr::summarise(pathway_gene_details = paste(pathway_line, collapse = "<br>"), .groups = "drop")
    } else {
        gene_tbl <- tibble::tibble(cell_group = character(), dysregulated_genes = character())
        pathway_tbl <- tibble::tibble(cell_group = character(), pathway_names = character())
        pathway_gene_tbl <- tibble::tibble(cell_group = character(), pathway_gene_details = character())
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
        effect_type        = if ("effect_type" %in% names(tab)) as.character(tab$effect_type) else NA_character_,
        sulston_display    = if ("sulston_display" %in% names(tab)) as.character(tab$sulston_display) else NA_character_,
        expectation        = if ("expectation" %in% names(tab)) as.character(tab$expectation) else NA_character_,
        rationale          = if ("rationale" %in% names(tab)) as.character(tab$rationale) else NA_character_
    )

    # Join dysregulated genes if available
    phenos <- phenos %>%
        dplyr::left_join(gene_tbl, by = "cell_group") %>%
        dplyr::left_join(pathway_tbl, by = "cell_group") %>%
        dplyr::left_join(pathway_gene_tbl, by = "cell_group")

    return(phenos)
}


phenotype_tooltip_builder <- function(g, render_mode = c("tissue", "global")) {
    render_mode <- match.arg(render_mode)
    g <- g %>%
        dplyr::mutate(
            group_str = dplyr::coalesce(
                if ("sulston_display" %in% names(g)) format_group_display_name(sulston_display) else NA_character_,
                format_group_display_name(group_nodes_by),
                ""
            ),
            abundance_code = if ("abundance_code" %in% names(g)) as.character(abundance_code) else NA_character_,
            lfc = if ("lfc" %in% names(g)) as.numeric(lfc) else NA_real_,
            q = if ("q" %in% names(g)) as.numeric(q) else NA_real_,
            power = if ("power" %in% names(g)) as.numeric(power) else NA_real_,
            powered_thresh = if ("powered_thresh" %in% names(g)) as.numeric(powered_thresh) else 0.8,
            present_above_thresh = if ("present_above_thresh" %in% names(g)) as.logical(present_above_thresh) else NA,
            present_above_thresh_flag = if ("present_above_thresh_flag" %in% names(g)) as.logical(present_above_thresh_flag) else dplyr::coalesce(present_above_thresh, FALSE),
            abundance_text = dplyr::case_when(
                (!is.na(abundance_code) & abundance_code == "A0 No change" & !is.na(power)) ~ "0%",
                ((is.na(abundance_code) | abundance_code == "") & !is.na(power)) ~ "0%",
                is.na(lfc) ~ "",
                TRUE ~ format_abundance_change_percent(lfc)
            ),
            abundance_power_status = abundance_power_status_label(power, powered_thresh),
            abundance_display = dplyr::case_when(
                abundance_text == "" ~ "",
                is.na(abundance_power_status) ~ abundance_text,
                TRUE ~ paste0(abundance_text, " (", abundance_power_status, ")")
            ),
            abundance_display_color = ifelse(
                abundance_text == "",
                NA_character_,
                ifelse(
                    !is.na(q) & q < 0.05 & abundance_text != "0%" & substr(abundance_text, 1, 1) == "-",
                    "#377eb8",
                    ifelse(!is.na(q) & q < 0.05 & abundance_text != "0%" & substr(abundance_text, 1, 1) == "+", "#e41a1c", NA_character_)
                )
            ),
            identity_tooltip_label = if ("identity_tooltip_label" %in% names(g)) as.character(identity_tooltip_label) else "",
            identity_display = if ("identity_display" %in% names(g)) as.character(identity_display) else as.character(identity_tooltip_label),
            glyph = dplyr::coalesce(as.character(glyph), ""),
            effect_type = if ("effect_type" %in% names(g)) as.character(effect_type) else "",
            expectation = if ("expectation" %in% names(g)) as.character(expectation) else "",
            expectation_text = if ("expectation_text" %in% names(g)) as.character(expectation_text) else "",
            expectation_norm = tolower(trimws(dplyr::coalesce(expectation, ""))),
            expectation_display = dplyr::case_when(
                identical(render_mode, "global") ~ expectation_text,
                expectation_text == "" ~ "",
                TRUE ~ stringr::str_replace(
                    expectation_text,
                    regex("see below for more information\\.?$", ignore_case = TRUE),
                    "See table for more information about expectation."
                )
            ),
            f1_dir = if ("f1_dir" %in% names(g)) as.character(f1_dir) else NA_character_,
            f2 = if ("f2" %in% names(g)) as.logical(f2) else FALSE,
            f3 = if ("f3" %in% names(g)) as.numeric(f3) else NA_real_,
            f4 = if ("f4" %in% names(g)) as.logical(f4) else FALSE,
            pathway_gene_details = if ("pathway_gene_details" %in% names(g)) as.character(pathway_gene_details) else NA_character_,
            pathway_names = if ("pathway_names" %in% names(g)) as.character(pathway_names) else NA_character_
        )

    paste0(
        "<b>Cell type: ", g$name, "</b>", ifelse(g$group_str != "", paste0(" | Tissue: ", g$group_str), ""), "<br>",
        ifelse(
            !is.na(g$present_above_thresh_flag) & !g$present_above_thresh_flag,
            "This cell type was not present at the timepoints sampled in this experimental design.<br>",
            ""
        ),
        ifelse(
            g$abundance_display != "",
            paste0(
                "<span",
                ifelse(!is.na(g$abundance_display_color), paste0(" style='color:", g$abundance_display_color, ";'"), ""),
                ">Abundance phenotype: ", g$abundance_display, "</span><br>"
            ),
            ""
        ),
        ifelse(
            !is.na(g$identity_display) & g$identity_display != "" & !(g$identity_display %in% c("Identity intact")),
            paste0(
                ifelse(g$glyph != "", paste0(g$glyph, " "), ""),
                stringr::str_to_lower(g$identity_display),
                ifelse(!is.na(g$effect_type) & g$effect_type != "", paste0(" (", g$effect_type, ")"), ""),
                "<br>"
            ),
            ""
        ),
        ifelse(
            g$expectation_display != "",
            paste0(
                ifelse(g$expectation_norm %in% c("expected", "expected change"), "○ ", ""),
                g$expectation_display,
                "<br>"
            ),
            ""
        ),
        ifelse(!is.na(g$f1_dir) & g$f1_dir == "increase", "▲ proliferation increase<br>", ""),
        ifelse(!is.na(g$f1_dir) & g$f1_dir == "decrease", "▼ proliferation decrease<br>", ""),
        ifelse(!is.na(g$f2) & g$f2, "─ apoptosis<br>", ""),
        ifelse(!is.na(g$f3) & g$f3 != 0, "* stress<br>", ""),
        ifelse(!is.na(g$f4) & g$f4, "□ senescence<br>", ""),
        ifelse(!identical(render_mode, "global") & !is.na(g$pathway_gene_details) & g$pathway_gene_details != "", g$pathway_gene_details, ""),
        ifelse(!identical(render_mode, "global") & (is.na(g$pathway_gene_details) | g$pathway_gene_details == "") & !is.na(g$pathway_names) & g$pathway_names != "", paste0("Pathway disrupted: ", stringr::str_trunc(g$pathway_names, 240)), "")
    )
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
                                   legend_position = "right",
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
                                   node_overlay = c("glyphs", "none"),
                                   glyph_color = "white",
                                   badge_color = "black",
                                   badge_outline_color = "white",
                                   render_mode = c("tissue", "global"),
                                   tooltip_builder = NULL,
                                   interactive = FALSE,
                                   width = NULL,
                                   height = NULL) {
    node_overlay <- match.arg(node_overlay)
    render_mode <- match.arg(render_mode)
    show_node_glyphs <- node_overlay %in% c("glyphs", "both")
    global_identity_marker <- identical(render_mode, "global")

    if (identical(render_mode, "global")) {
        show_node_glyphs <- FALSE
    }

    g <- cell_state_graph@g %>% dplyr::mutate(name = stringr::str_trim(name))
    g <- g %>%
        select(x, y, name, color_nodes_by, label_nodes_by, group_nodes_by) %>%
        distinct()
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
            power = if ("power" %in% names(phenos_df)) as.numeric(.data[["power"]]) else NA_real_,
            powered_thresh = if ("powered_thresh" %in% names(phenos_df)) as.numeric(.data[["powered_thresh"]]) else 0.8,
            present_above_thresh = if ("present_above_thresh" %in% names(phenos_df)) as.logical(.data[["present_above_thresh"]]) else NA,
            abundance_code = if ("abundance_code" %in% names(phenos_df)) as.character(.data[["abundance_code"]]) else NA_character_,
            ident = if (map$ident %in% names(phenos_df)) as.character(.data[[map$ident]]) else "I0",
            glyph = if (map$glyph %in% names(phenos_df)) as.character(.data[[map$glyph]]) else "",
            f1_dir = if (map$f1 %in% names(phenos_df)) as.character(.data[[map$f1]]) else NA_character_,
            fitness_label = if ("fitness_label" %in% names(phenos_df)) as.character(.data[["fitness_label"]]) else NA_character_,
            f2 = if (map$f2 %in% names(phenos_df)) as.logical(.data[[map$f2]]) else NA,
            f3 = if (map$f3 %in% names(phenos_df)) as.numeric(.data[[map$f3]]) else NA_real_,
            f4 = if (map$f4 %in% names(phenos_df)) as.logical(.data[[map$f4]]) else NA,
            dysregulated_genes = if ("dysregulated_genes" %in% names(phenos_df)) as.character(.data[["dysregulated_genes"]]) else NA_character_,
            pathway_names = if ("pathway_names" %in% names(phenos_df)) as.character(.data[["pathway_names"]]) else NA_character_,
            pathway_gene_details = if ("pathway_gene_details" %in% names(phenos_df)) as.character(.data[["pathway_gene_details"]]) else NA_character_,
            identity_evidence = if ("identity_evidence" %in% names(phenos_df)) as.character(.data[["identity_evidence"]]) else NA_character_,
            stress_evidence = if ("stress_evidence" %in% names(phenos_df)) as.character(.data[["stress_evidence"]]) else NA_character_,
            effect_type = if ("effect_type" %in% names(phenos_df)) as.character(.data[["effect_type"]]) else NA_character_,
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
            identity_tooltip_label = dplyr::case_when(
                is.na(ident) ~ "",
                ident == "I0" ~ "",
                ident == "I1" ~ "maturation delay",
                ident == "I2" ~ "precocious maturation",
                ident == "I3" ~ "identity",
                ident == "I4" ~ "fate switch / misspecification",
                ident == "I5" ~ "identity fragmentation",
                TRUE ~ tolower(identity_display)
            ),
            abundance_str = show_pheno(abundance_code, c("A0 No change", NA)),
            identity_str = show_pheno(identity_display, c("I0", "I0 Identity intact", "Identity intact", NA)),
            group_str = dplyr::coalesce(as.character(group_nodes_by), ""),
            abundance_dir = dplyr::case_when(
                !is.na(lfc) & lfc < 0 ~ "decrease",
                !is.na(lfc) & lfc > 0 ~ "increase",
                !is.na(abundance_code) & grepl("^A2|^A3", abundance_code) ~ "decrease",
                !is.na(abundance_code) & grepl("^A1|^A4", abundance_code) ~ "increase",
                TRUE ~ ""
            ),
            abundance_fc = ifelse(!is.na(lfc), exp(abs(lfc)), NA_real_),
            abundance_text = dplyr::case_when(
                abundance_str == "" ~ "",
                !is.na(abundance_fc) & abundance_dir != "" ~ paste0(formatC(abundance_fc, digits = 1, format = "f"), "x ", abundance_dir),
                abundance_dir != "" ~ abundance_dir,
                TRUE ~ ""
            ),
            expectation_norm = tolower(trimws(dplyr::coalesce(expectation, ""))),
            expectation_text = dplyr::case_when(
                expectation_norm %in% c("expected", "expected change") ~ "Expectation: a phenotype is expected, see tissue specific plot for more information",
                expectation_norm != "" ~ "Expectation: no phenotype expected, see tissue specific plot for more information",
                TRUE ~ ""
            ),
            fitness_str = show_pheno(f1_dir, c("F0 No significant phenotype", "F0 Normal", NA)),
            badge_caption = paste0(
                ifelse(!is.na(f1_dir) & f1_dir == "increase", "F1_up, ", ""),
                ifelse(!is.na(f1_dir) & f1_dir == "decrease", "F1_down, ", ""),
                ifelse(!is.na(f2) & f2, "F2, ", ""),
                ifelse(!is.na(f3) & f3 != 0, "F3, ", ""),
                ifelse(!is.na(f4) & f4, "F4, ", "")
            ),
            badge_caption = stringr::str_replace(badge_caption, ",\\s*$", ""),
            lfc_capped = pmin(pmax(lfc, -lfc_cap), lfc_cap),
            glyph = dplyr::coalesce(glyph, ""),
            glyph_col = ifelse(is.na(lfc), "black", ifelse(abs(lfc) >= 0.8, "white", "black")),
            f3_alpha = dplyr::case_when(
                is.na(f3) ~ 0,
                TRUE ~ scales::rescale(pmin(abs(f3), stress_cap), to = c(0.25, 1))
            ),
            power_status = dplyr::case_when(
                !is.na(power) & power >= powered_thresh ~ "Powered",
                TRUE ~ "Underpowered"
            ),
            node_size_plot = dplyr::case_when(
                identical(render_mode, "global") & power_status == "Powered" ~ node_size * 3.8,
                identical(render_mode, "global") ~ node_size * 1.9,
                power_status == "Powered" ~ node_size * 3.2,
                TRUE ~ node_size * 2.0
            )
        )

    g$.tooltip <- NA_character_
    if (isTRUE(interactive) && !is.null(tooltip_builder)) {
        tooltip_vec <- tooltip_builder(g, render_mode)
        if (length(tooltip_vec) != nrow(g)) {
            stop("tooltip_builder must return one tooltip per plotted node.")
        }
        g$.tooltip <- as.character(tooltip_vec)
    } else if (isTRUE(interactive)) {
        g$.tooltip <- phenotype_tooltip_builder(g)
    }

    color_priority <- c("severe", "medium", "mild", "none", "not_present")
    g <- g %>%
        dplyr::mutate(
            has_fitness_change = (!is.na(f1_dir) & f1_dir %in% c("increase", "decrease")) |
                dplyr::coalesce(f2, FALSE) |
                (!is.na(f3) & f3 != 0) |
                dplyr::coalesce(f4, FALSE),
            primary_phenotype = dplyr::case_when(
                !is.na(abundance_code) & grepl("^A2|^A3", abundance_code) ~ "abundance_loss",
                !is.na(ident) & !(ident %in% c("I0", "I0 Identity intact", "Identity intact", "", NA)) ~ "identity",
                !is.na(abundance_code) & grepl("^A1|^A4", abundance_code) ~ "abundance_gain",
                has_fitness_change ~ "fitness",
                TRUE ~ "none"
            ),
            severity_fill = dplyr::case_when(
                primary_phenotype == "abundance_loss" ~ "severe",
                primary_phenotype == "identity" ~ "medium",
                primary_phenotype %in% c("abundance_gain", "fitness") ~ "mild",
                TRUE ~ "none"
            ),
            has_identity_change = !is.na(ident) & !(ident %in% c("I0", "I0 Identity intact", "Identity intact", "", NA)),
            expectation_norm = tolower(trimws(dplyr::coalesce(expectation, ""))),
            is_expected = expectation_norm %in% c("expected", "expected change"),
            expected_shape = dplyr::if_else(is_expected, "Expected phenotype", "Observed phenotype"),
            present_above_thresh_flag = dplyr::coalesce(present_above_thresh, TRUE),
            effect_type_norm = tolower(trimws(dplyr::coalesce(effect_type, ""))),
            autonomy_display = dplyr::case_when(
                effect_type_norm %in% c("cell-autonomous", "cell autonomous") ~ "Cell-autonomous",
                stringr::str_detect(effect_type_norm, "non-autonomous|non autonomous") ~ "Non-autonomous",
                TRUE ~ "Cell-autonomous"
            )
        )

    g <- g %>% mutate(severity_fill = ifelse(present_above_thresh_flag == FALSE, "not_present", severity_fill))

    if (identical(render_mode, "global")) {
        # Mild horizontal spread keeps the global plot readable without shrinking nodes.
        g <- g %>% dplyr::mutate(x = x * 1.35)
        if ("x" %in% names(bezier_df)) {
            bezier_df <- bezier_df %>% dplyr::mutate(x = x * 1.35)
        }
    }

    g_draw <- g %>%
        dplyr::mutate(.draw_order = ifelse(primary_phenotype == "none", 0L, 1L)) %>%
        dplyr::arrange(.draw_order)

    edge_arrow_unit <- if (identical(render_mode, "tissue")) max(arrow_unit, 4) else max(arrow_unit, 5)
    edge_linewidth <- if (identical(render_mode, "tissue")) 0.35 else 0.45
    edge_colour <- if (identical(render_mode, "global")) "#6f6f6f" else con_colour

    p <- ggplot2::ggplot(ggplot2::aes(x, y), data = g) +
        ggplot2::geom_path(
            ggplot2::aes(x, y, group = edge_name),
            linewidth = edge_linewidth, colour = edge_colour, data = dplyr::distinct(bezier_df),
            arrow = ggplot2::arrow(angle = 30, length = grid::unit(edge_arrow_unit, "pt"), type = "closed"),
            linejoin = "mitre", lineend = "round"
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
        p <- p + ggrepel::geom_text_repel(
            data = group_labels,
            ggplot2::aes(x = x, y = y, label = group_nodes_by),
            size = group_label_font_size,
            color = "black",
            box.padding = 0.3,
            point.padding = 0.1,
            segment.color = "grey50"
        )
    }


    # Single node fill with explicit priority: loss > identity/gain > fitness > none.
    if (isTRUE(interactive)) {
        p <- p +
            ggiraph::geom_point_interactive(
                ggplot2::aes(x = x, y = y, tooltip = .tooltip, fill = severity_fill, shape = expected_shape, size = power_status, alpha = autonomy_display),
                data = g_draw,
                color = "black",
                stroke = if (identical(render_mode, "global")) 0.7 else 0.5
            )
    } else {
        p <- p +
            ggplot2::geom_point(
                ggplot2::aes(x = x, y = y, fill = severity_fill, shape = expected_shape, size = power_status, alpha = autonomy_display),
                data = g_draw,
                color = "black",
                stroke = if (identical(render_mode, "global")) 0.7 else 0.5
            )
    }

    p <- p +
        ggplot2::scale_alpha_manual(
            name = "Autonomy",
            values = c(
                "Cell-autonomous" = 1,
                "Non-autonomous" = 0.45
            ),
            breaks = c("Cell-autonomous", "Non-autonomous"),
            drop = FALSE
        ) +
        ggplot2::scale_size_identity(guide = "none") +
        ggplot2::scale_fill_manual(
            values = c(
                severe = unname(phenotype_colors["abundance_loss"]),
                medium = unname(phenotype_colors["identity"]),
                mild = unname(phenotype_colors["abundance_gain"]),
                none = unname(phenotype_colors["none"]),
                not_present = unname(phenotype_colors["not_present"])
            ),
            breaks = color_priority,
            labels = c(
                severe = "Severe phenotype",
                medium = "Medium phenotype",
                mild = "Mild phenotype",
                none = "No phenotype",
                not_present = "Not present"
            ),
            name = if (identical(render_mode, "global")) NULL else "Node color",
            guide = "none"
        )

    fill_legend_df <- data.frame(
        x = NA_real_,
        y = NA_real_,
        severity_fill = factor(
            c("severe", "medium", "mild", "none", "not_present"),
            levels = c("severe", "medium", "mild", "none", "not_present")
        )
    )

    p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(
            data = fill_legend_df,
            ggplot2::aes(x = x, y = y, fill = severity_fill),
            inherit.aes = FALSE,
            shape = 21,
            size = 4,
            colour = "black",
            stroke = 0.5,
            show.legend = c(color = TRUE, size = FALSE)
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                severe = unname(phenotype_colors["abundance_loss"]),
                medium = unname(phenotype_colors["identity"]),
                mild = unname(phenotype_colors["abundance_gain"]),
                none = unname(phenotype_colors["none"]),
                not_present = unname(phenotype_colors["not_present"])
            ),
            breaks = c("severe", "medium", "mild", "none", "not_present"),
            labels = c(
                severe = "Severe phenotype",
                medium = "Medium phenotype",
                mild = "Mild phenotype",
                none = "No phenotype",
                not_present = "Not present"
            ),
            drop = FALSE,
            name = "Node color",
            guide = ggplot2::guide_legend(
                order = 1,
                override.aes = list(
                    shape = 21,
                    size = 4,
                    colour = "black",
                    stroke = 0.5,
                    alpha = 1
                )
            )
        )

    size_legend_df <- data.frame(
        x = NA_real_,
        y = NA_real_,
        size_value = c("Powered", "Underpowered")
    )

    p <- p +
        ggplot2::geom_point(
            data = size_legend_df,
            ggplot2::aes(x = x, y = y, size = size_value),
            inherit.aes = FALSE,
            shape = 21,
            fill = "black",
            stroke = 0.5,
            show.legend = c(color = FALSE, size = TRUE)
        ) +
        ggplot2::scale_size_manual(
            values = c("Powered" = node_size * 1.6, "Underpowered" = node_size * 0.8),
            guide = guide_legend(
                override.aes = list(size = c(1.6 * 2, 0.8 * 2))
            ),
            name = "Size"
        )

    shape_legend_df <- data.frame(
        x = NA_real_,
        y = NA_real_,
        shape_value = c("Observed phenotype", "Expected phenotype")
    )

    if ("Expected phenotype" %in% unique(g_draw$expected_shape)) {
        p <- p +
            ggplot2::geom_point(
                data = shape_legend_df,
                ggplot2::aes(x = x, y = y, shape = shape_value),
                inherit.aes = FALSE,
                color = "black",
                fill = "white",
                stroke = 0.5,
                size = 3,
                show.legend = c(color = FALSE, size = TRUE)
            ) +
            ggplot2::scale_shape_manual(
                values = c("Observed phenotype" = 21, "Expected phenotype" = 22),
                name = "Shape"
            )
    } else {
        p <- p +
            ggplot2::scale_shape_manual(
                values = c("Observed phenotype" = 21, "Expected phenotype" = 22),
                name = "Shape"
            ) + guides(shape = "none")
    }


    if (isTRUE(global_identity_marker)) {
        identity_marker_legend_df <- data.frame(
            x = NA_real_,
            y = NA_real_,
            identity_marker_label = "Phenotype detected",
            identity_marker_glyph = "*",
            stringsAsFactors = FALSE
        )
        p <- p +
            geom_text(
                data = g,
                aes(x = x, y = y, label = ifelse(has_identity_change, "*", "")),
                show.legend = FALSE,
                color = glyph_color,
                size = node_size,
                vjust = 0.8,
                hjust = 0.5
            ) +
            ggnewscale::new_scale_color() +
            geom_text(
                data = identity_marker_legend_df,
                aes(x = x, y = y, label = identity_marker_glyph, color = identity_marker_label),
                inherit.aes = FALSE,
                show.legend = TRUE,
                size = node_size,
                vjust = 0.8,
                hjust = 0.5
            ) +
            scale_color_manual(
                name = "Transcriptional Identity",
                values = c("Phenotype detected" = "black"),
                guide = guide_legend(
                    override.aes = list(label = "*", size = 6),
                    order = 1
                )
            )
    }

    if (show_node_glyphs) {
        glyph_draw_color <- if (is.null(glyph_color)) g$glyph_col else glyph_color
        g_draw$glyph_size <- case_when(
            g_draw$power_status == "Powered" ~ node_size * 1.6 * 0.5,
            TRUE ~ node_size * 0.8 * 0.6
        )
        if (isTRUE(interactive)) {
            p <- p +
                ggnewscale::new_scale("size") +
                ggiraph::geom_text_interactive(
                    ggplot2::aes(x = x, y = y, label = glyph, tooltip = .tooltip, size = glyph_size, alpha = autonomy_display),
                    data = g_draw,
                    color = glyph_draw_color,
                    fontface = "bold", vjust = 0.43, hjust = 0.5,
                    show.legend = FALSE
                ) + ggplot2::scale_size_identity(guide = "none")
        } else {
            p <- p +
                ggnewscale::new_scale("size") +
                ggplot2::geom_text(
                    ggplot2::aes(x = x, y = y, label = glyph, size = glyph_size, alpha = autonomy_display),
                    data = g_draw,
                    color = glyph_draw_color,
                    fontface = "bold", vjust = 0.43, hjust = 0.5, show.legend = FALSE
                ) + ggplot2::scale_size_identity(guide = "none")

            glyph_legend_df <- data.frame(
                x = NA_real_,
                y = NA_real_,
                glyph = c("", "<<", ">>", "!!", "<>", "##"),
                glyph_type = c(
                    "Identity intact",
                    "Maturation delay",
                    "Precocious maturation",
                    "Program failure",
                    "Fate switch / misspecification",
                    "Identity fragmentation"
                ),
                stringsAsFactors = FALSE
            )
            glyph_order <- c(
                "Identity intact",
                "Maturation delay",
                "Precocious maturation",
                "Program failure",
                "Fate switch / misspecification",
                "Identity fragmentation"
            )

            glyph_legend_df$glyph_type <- factor(
                glyph_legend_df$glyph_type,
                levels = glyph_order
            )

            p <- p +
                ggnewscale::new_scale_color() +
                geom_text(
                    data = glyph_legend_df,
                    aes(x = x, y = y, label = glyph, color = glyph_type),
                    inherit.aes = FALSE,
                    fontface = "bold",
                    show.legend = TRUE
                ) +
                scale_color_manual(
                    name = "Glyphs",
                    breaks = glyph_order,
                    values = stats::setNames(rep("black", length(glyph_order)), glyph_order),
                    labels = c(
                        "Identity intact",
                        "Maturation delay",
                        "Precocious maturation",
                        "Program failure",
                        "Fate switch / misspecification",
                        "Identity fragmentation"
                    ),
                    guide = guide_legend(
                        override.aes = list(
                            label = c("", "<<", ">>", "!!", "<>", "##"),
                            size = 3,
                            colour = "black"
                        ),
                        order = 1
                    )
                )
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

    if (isTRUE(interactive)) {
        if (is.null(width) || is.null(height)) {
            p <- p + ggplot2::coord_equal() + theme(legend.position = legend_position)
            ggiraph::girafe(
                ggobj = p
            )
        } else {
            p <- p + theme(legend.position = legend_position)
            p <- ggiraph::girafe(
                ggobj = p,
                width_svg = width, height_svg = height
            )
        }
    } else {
        p + theme(legend.position = legend_position)
    }
}
