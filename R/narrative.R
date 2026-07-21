load_ai_precompute_for_cell_types <- function(cell_types, ai_notes_path = "../../../ai_notes/cell_types/") {
  breakFun <- function(x) {
    if (nchar(x) == 0) {
      return("\n\n")
    } else {
      return(x)
    }
  }

  read_cell_type_background <- function(ct, ai_notes_path, bg_type = "genetic_req.md") {
    ai_bg_path <- paste(ai_notes_path, blogdown:::dash_filename(ct), bg_type, sep = "/")
    if (fs::file_exists(ai_bg_path)) {
      storeLines <- readLines(ai_bg_path, warn = FALSE)
      bg_content <- paste0(storeLines, collapse = "\n")
      return(bg_content)
    } else {
      message(paste("No", bg_type, "file found for cell type", ct))
      return(NA_character_)
    }
  }

  # Load sig_pathways.rds for each cell type
  read_sig_pathways <- function(ct, ai_notes_path) {
    rds_path <- file.path(ai_notes_path, blogdown:::dash_filename(ct), "sig_pathways.rds")
    if (fs::file_exists(rds_path)) {
      tryCatch(
        readRDS(rds_path),
        error = function(e) {
          message(paste("Error loading sig_pathways.rds for", ct, ":", e$message))
          return(NULL)
        }
      )
    } else {
      message(paste("No sig_pathways.rds file found for cell type", ct))
      return(NULL)
    }
  }

  ai_gen_req_bg <- lapply(cell_types, read_cell_type_background, ai_notes_path = ai_notes_path, bg_type = "genetic_req.md")
  names(ai_gen_req_bg) <- cell_types

  ai_kinetic_bg <- lapply(cell_types, read_cell_type_background, ai_notes_path = ai_notes_path, bg_type = "kinetic_summary.md")
  names(ai_kinetic_bg) <- cell_types

  ai_lineage_bg <- lapply(cell_types, read_cell_type_background, ai_notes_path = ai_notes_path, bg_type = "lineage_summary.md")
  names(ai_lineage_bg) <- cell_types

  ai_pathway_bg <- lapply(cell_types, read_cell_type_background, ai_notes_path = ai_notes_path, bg_type = "pathway_summary.md")
  names(ai_pathway_bg) <- cell_types

  ai_regulator_pathway_bg <- lapply(cell_types, read_cell_type_background, ai_notes_path = ai_notes_path, bg_type = "regulator_pathway_summary.md")
  names(ai_regulator_pathway_bg) <- cell_types

  ai_summary_bg <- lapply(cell_types, read_cell_type_background, ai_notes_path = ai_notes_path, bg_type = "summary.md")
  names(ai_summary_bg) <- cell_types

  # Load sig_pathways.rds for each cell type
  sig_pathways <- lapply(cell_types, read_sig_pathways, ai_notes_path = ai_notes_path)
  names(sig_pathways) <- cell_types

  bg_tibble <- tibble(
    cell_type = names(ai_gen_req_bg),
    genetic_req_background = unlist(ai_gen_req_bg),
    kinetic_background = unlist(ai_kinetic_bg),
    lineage_background = unlist(ai_lineage_bg),
    pathway_background = unlist(ai_pathway_bg),
    regulator_pathway_background = unlist(ai_regulator_pathway_bg),
    summary_background = unlist(ai_summary_bg),
    sig_pathways = sig_pathways
  )

  return(bg_tibble)
}

# Repeated-glyph arrow encoding the MAGNITUDE of a fold change, so larger
# effects get more arrows and a gene list sorted by |LFC| reads with a
# monotonically shrinking arrow stack. Direction comes from dysreg_type (or the
# sign of the shrunken log2 fold change); the count is bucketed by |log2FC|
# against `breaks` (default: <1 -> 1 arrow, 1-2 -> 2, >=2 -> 3). No detectable
# direction -> "" (no arrow). Vectorized. Kept as one helper so format_goi and
# format_pathway_regulatory_genes speak the same arrow language.
magnitude_arrow <- function(lfc, dysreg_type = NA_character_, breaks = c(1, 2)) {
  up <- "\u2191"
  down <- "\u2193"
  dir <- dplyr::case_when(
    !is.na(dysreg_type) & dysreg_type == "Overexpressed" ~ 1,
    !is.na(dysreg_type) & dysreg_type == "Underexpressed" ~ -1,
    !is.na(lfc) & lfc > 0 ~ 1,
    !is.na(lfc) & lfc < 0 ~ -1,
    TRUE ~ 0
  )
  n <- ifelse(
    is.na(lfc), 1L,
    1L + as.integer(abs(lfc) >= breaks[1]) + as.integer(abs(lfc) >= breaks[2])
  )
  glyph <- dplyr::case_when(dir > 0 ~ up, dir < 0 ~ down, TRUE ~ "")
  ifelse(glyph == "", "", strrep(glyph, n))
}

# Build the per-cell "genes of interest" line. Two compartments -- regulators
# (transcription factors, i.e. genes in `genes_of_interest`) first, then all
# other dysregulated genes -- so the prose the LLM writes (which may name any
# significant DEG, not just regulators) is anchored to something the reader can
# see. Within each compartment genes are sorted by |log2FC| descending (largest
# effect first; NA/zero LFC last, alphabetical), carry magnitude arrows, and are
# capped at `max_per_compartment` with a trailing "[+K more]".
format_goi <- function(cell_type_degs, genes_of_interest = NULL,
                       max_per_compartment = 15) {
  empty <- list(goi_vector_sorted = character(0), goi_line = "none")
  if (
    is.null(cell_type_degs) ||
      !is.data.frame(cell_type_degs) ||
      nrow(cell_type_degs) == 0 ||
      !"gene_short_name" %in% colnames(cell_type_degs)
  ) {
    return(empty)
  }
  df <- cell_type_degs
  if (!"dysreg_type" %in% colnames(df)) df$dysreg_type <- NA_character_
  if (!"perturb_to_ctrl_shrunken_lfc" %in% colnames(df)) {
    df$perturb_to_ctrl_shrunken_lfc <- NA_real_
  }
  df <- df %>%
    dplyr::filter(!is.na(gene_short_name), nzchar(gene_short_name)) %>%
    dplyr::distinct(gene_short_name, .keep_all = TRUE) %>%
    dplyr::mutate(
      lfc = perturb_to_ctrl_shrunken_lfc,
      abs_lfc = abs(lfc),
      arrow = magnitude_arrow(lfc, dysreg_type),
      gene_arrow = paste0(gene_short_name, arrow),
      # Case-insensitive: some regulators are listed lowercase in the panel but
      # appear uppercase in the DEG table (e.g. RNF14 vs rnf14); a case-sensitive
      # match wrongly dropped them into the non-regulator compartment.
      is_regulator = if (!is.null(genes_of_interest)) {
        tolower(gene_short_name) %in% tolower(genes_of_interest)
      } else {
        FALSE
      }
    )
  if (nrow(df) == 0) return(empty)

  render_compartment <- function(sub) {
    if (nrow(sub) == 0) return(character(0))
    # arrange() sends NA abs_lfc to the end even under desc(), so NA/zero-effect
    # genes fall to the bottom of the compartment with an alphabetical tiebreak.
    genes <- sub %>%
      dplyr::arrange(dplyr::desc(abs_lfc), gene_short_name) %>%
      dplyr::pull(gene_arrow)
    if (length(genes) > max_per_compartment) {
      extra <- length(genes) - max_per_compartment
      genes <- c(genes[seq_len(max_per_compartment)], paste0("[+", extra, " more]"))
    }
    genes
  }

  reg <- render_compartment(df %>% dplyr::filter(is_regulator))
  oth <- render_compartment(df %>% dplyr::filter(!is_regulator))
  goi_vector_sorted <- c(reg, oth)
  if (length(goi_vector_sorted) == 0) return(empty)
  list(
    goi_vector_sorted = goi_vector_sorted,
    goi_line = paste(goi_vector_sorted, collapse = ", ")
  )
}

format_pathway_regulatory_genes <- function(goi_in_pathway, degs_of_interest) {
  if (is.null(degs_of_interest) ||
    !is.data.frame(degs_of_interest) ||
    !"gene_short_name" %in% colnames(degs_of_interest) ||
    is.null(goi_in_pathway) || length(goi_in_pathway) == 0) {
    return("")
  }
  # Ensure columns exist
  if (!"dysreg_type" %in% colnames(degs_of_interest)) {
    degs_of_interest$dysreg_type <- NA_character_
  }
  if (!"perturb_to_ctrl_shrunken_lfc" %in% colnames(degs_of_interest)) {
    degs_of_interest$perturb_to_ctrl_shrunken_lfc <- NA_real_
  }
  goi_df <- degs_of_interest %>%
    filter(gene_short_name %in% goi_in_pathway) %>%
    mutate(
      lfc = perturb_to_ctrl_shrunken_lfc,
      arrow = magnitude_arrow(lfc, dysreg_type),
      gene_arrow = paste0(gene_short_name, arrow)
    ) %>%
    # |LFC| descending so the arrow stack shrinks down the list (NA last).
    arrange(dplyr::desc(abs(lfc)), gene_short_name)
  paste(goi_df$gene_arrow, collapse = ", ")
}

# Strip the leading A#/F#/I# code token from a phenotype label so an ancestor's
# phenotype, when quoted in a no-signal cell's summary, reads in plain words --
# the reader never sees the coding scheme. "F1 Proliferation change" ->
# "proliferation change"; a plain word like "depleted" is returned unchanged.
humanize_pheno_label <- function(x) {
  x <- as.character(x)
  tolower(sub("^\\s*[AFI][0-9]\\s+", "", x))
}

# Build a gene-allowlist map for the LLM impact call.
# Returns a named list keyed by lowercase gene symbol, with values
# "up" / "down" / "any" based on the sign of perturb_to_ctrl_shrunken_lfc.
# Empty list when the DEG tibble is absent or has no gene column.
build_allowed_degs_map <- function(deg_tbl) {
  if (is.null(deg_tbl) || !is.data.frame(deg_tbl) || nrow(deg_tbl) == 0) {
    return(list())
  }
  if (!"gene_short_name" %in% colnames(deg_tbl)) {
    return(list())
  }
  lfc_col <- "perturb_to_ctrl_shrunken_lfc"
  tbl <- deg_tbl %>%
    dplyr::filter(!is.na(gene_short_name), nzchar(gene_short_name)) %>%
    dplyr::distinct(gene_short_name, .keep_all = TRUE)
  if (nrow(tbl) == 0) {
    return(list())
  }
  direction <- if (lfc_col %in% colnames(tbl)) {
    lfc <- tbl[[lfc_col]]
    dplyr::case_when(
      is.na(lfc) ~ "any",
      lfc > 0 ~ "up",
      lfc < 0 ~ "down",
      TRUE ~ "any"
    )
  } else {
    rep("any", nrow(tbl))
  }
  setNames(as.list(direction), tolower(tbl$gene_short_name))
}

# Helper: Build pathway_tbl for an ancestor
build_pathway_tbl <- function(info, sig_p_val_thresh, top_n_pathways) {
  pathway_tbl <- tibble()
  if (!is.null(info$deg_enrichment_in_specific_pathways) && nrow(info$deg_enrichment_in_specific_pathways) > 0) {
    sig_pathways <- info$deg_enrichment_in_specific_pathways %>%
      filter(pval < sig_p_val_thresh) %>%
      group_by(ontology) %>%
      slice_min(order_by = pval, n = top_n_pathways, with_ties = FALSE) %>%
      arrange(ontology, pval)
    if (nrow(sig_pathways) > 0) {
      pathway_tbl <- sig_pathways %>%
        mutate(
          regulatory_genes = purrr::map_chr(
            genes_of_interest_in_pathway,
            ~ format_pathway_regulatory_genes(.x, info$degs_of_interest)
          )
        )
    }
  }
  pathway_tbl
}

get_descendants <- function(cell_type, combined_psg) {
  g <- coerce_state_graph(combined_psg)
  if (!cell_type %in% igraph::V(g)$name) {
    return(character(0))
  }
  igraph::subcomponent(g, cell_type, mode = "out") %>%
    names() %>%
    setdiff(cell_type)
}

get_roots <- function(ct, combined_psg) {
  g <- coerce_state_graph(combined_psg)
  # Find all vertices in the connected component containing ct
  if (!ct %in% igraph::V(g)$name) {
    return(character(0))
  }
  component_nodes <- igraph::subcomponent(g, ct, mode = "all") %>% names()
  # Roots are nodes in this component with no incoming edges
  roots <- component_nodes[
    sapply(component_nodes, function(v) length(igraph::neighbors(g, v, mode = "in")) == 0)
  ]
  roots
}

humanize_pathway_name <- function(pathway) {
  # Remove ontology prefix (e.g., "GOBP_", "GOMF_", "GOCC_")
  pathway <- sub("^(GOBP_|GOMF_|GOCC_)", "", pathway)
  # Replace underscores with spaces
  pathway <- gsub("_", " ", pathway)
  # Make lower case, then capitalize first letter
  pathway <- stringr::str_to_sentence(pathway)
  pathway
}

summarize_cell_type_impact <- function(
  ct,
  perturbation_description,
  target_gene_expression,
  dact_results,
  degs,
  ref_expression,
  ontologies,
  ai_notes_path,
  sig_p_val_thresh = 0.05,
  genes_of_interest = NULL,
  empirical_p_thresh = 1.0,
  power_thresh = 0.8,
  top_n_pathways = 5,
  abundance_phenotypes = NULL,
  fitness_phenotypes = NULL,
  identity_phenotypes = NULL # <-- NEW ARGUMENT
) {
  # 1. Abundance change (from abundance_phenotypes if available)
  abundance_row <- if (!is.null(abundance_phenotypes)) {
    abundance_phenotypes %>% filter(cell_group == ct)
  } else {
    tibble()
  }
  abundance_code <- if (nrow(abundance_row) > 0) abundance_row$abundance_code[1] else NA_character_
  abundance_severity <- if (nrow(abundance_row) > 0) abundance_row$abundance_severity[1] else NA_character_

  # Fallback to dact_results if abundance_phenotypes is not available
  dact_row <- dact_results %>% filter(cell_group == ct)
  abundance_summary <- if (!is.na(abundance_code)) {
    abundance_code
  } else if (nrow(dact_row) > 0 && all(c("delta_q_value", "delta_log_abund", "power") %in% colnames(dact_row))) {
    if (dact_row$delta_q_value[1] < sig_p_val_thresh) {
      if (dact_row$delta_log_abund[1] < 0) {
        "depleted"
      } else if (dact_row$delta_log_abund[1] > 0) {
        "enriched"
      } else {
        "significant_no_direction"
      }
    } else {
      if (dact_row$power[1] >= power_thresh) {
        "no_change"
      } else {
        "unknown (underpowered)"
      }
    }
  } else {
    NA_character_
  }

  # 2. Fitness phenotype (from fitness_phenotypes if available)
  fitness_rows <- if (!is.null(fitness_phenotypes)) {
    fitness_phenotypes %>% filter(cell_group == ct)
  } else {
    tibble()
  }
  fitness_label <- if (nrow(fitness_rows) > 0) paste(unique(fitness_rows$fitness_label), collapse = "; ") else NA_character_
  fitness_severity <- if (nrow(fitness_rows) > 0) paste(unique(fitness_rows$severity), collapse = "; ") else NA_character_
  fitness_evidence <- if (nrow(fitness_rows) > 0) paste(unique(fitness_rows$evidence), collapse = "; ") else NA_character_

  # 3. Identity/maturation phenotype (from identity_phenotypes if available)
  identity_rows <- if (!is.null(identity_phenotypes)) {
    identity_phenotypes %>% filter(cell_group == ct)
  } else {
    tibble()
  }
  identity_label <- if (nrow(identity_rows) > 0) paste(unique(identity_rows$identity_label), collapse = "; ") else NA_character_
  identity_evidence <- if (nrow(identity_rows) > 0) paste(unique(identity_rows$evidence), collapse = "; ") else NA_character_

  # 4. DEGs and key regulatory genes.
  # Foreground for BOTH the FORA/GO enrichment and the named genes-of-interest is
  # the artifact-aware set: genes whose standardized effect stands out beyond the
  # control-sampling-artifact null (empirical_p <= empirical_p_thresh). Unlike a
  # nominal-p foreground, this drops the low-expression artifact calls that
  # otherwise drown real on-lineage programs; unlike the BH empirical_fdr, it
  # resolves per gene and does not hit the per-cell FDR floor. Falls back to
  # nominal p when the empirical_p column is absent (undecorated tables).
  cell_type_degs <- degs %>% filter(cell_group == ct)
  # Whether the abundance analysis judged this cell type reliably present in the
  # experiment (longest-contiguous window above the abundance threshold; see
  # assembly_utils.R). Captured from the full cell DEGs BEFORE the empirical_p
  # filter so it survives even when the significant set is empty. Cells that are
  # not present-above-threshold cannot be captured by this perturbation and are
  # dropped from the impact table downstream.
  cell_present_above_thresh <- "present_above_thresh" %in% names(cell_type_degs) &&
    isTRUE(any(cell_type_degs$present_above_thresh, na.rm = TRUE))
  if (!is.null(empirical_p_thresh) && empirical_p_thresh < 1 &&
      "empirical_p" %in% colnames(cell_type_degs)) {
    cell_type_degs <- cell_type_degs %>%
      filter(!is.na(empirical_p), empirical_p <= empirical_p_thresh)
  } else {
    cell_type_degs <- cell_type_degs %>%
      filter(perturb_to_ctrl_p_value < sig_p_val_thresh)
  }
  degs_of_interest <- if (!is.null(genes_of_interest)) {
    cell_type_degs %>% filter(tolower(gene_short_name) %in% tolower(genes_of_interest))
  } else {
    cell_type_degs
  }
  # goi_line now spans ALL significant DEGs (regulators first, then the rest) so
  # the LLM prose -- which may name any significant DEG -- stays anchored to
  # something the reader sees. `degs_of_interest` (regulators only) is still used
  # below for pathway regulatory-gene chips and for the has_regulator_goi gate.
  # max_per_compartment is tunable (see plan open items).
  goi <- format_goi(cell_type_degs, genes_of_interest, max_per_compartment = 15)

  # 5. Pathway enrichment (using precomputed or run on the fly)
  precompute_data <- load_ai_precompute_for_cell_types(cell_types = ct, ai_notes_path = ai_notes_path)
  pathway_background <- precompute_data$pathway_background
  sig_pathways <- precompute_data$sig_pathways[[1]]
  if (is.null(sig_pathways) || !is.data.frame(sig_pathways)) {
    sig_pathways <- tibble::tibble(pathway = character())
  }
  ancestor_specific_pathways <- sig_pathways %>%
    distinct(pathway) %>%
    rename(gs_name = pathway)
  deg_genes <- cell_type_degs$gene_short_name
  fora_res_list <- list()
  for (ont in names(ontologies)) {
    gene_set <- ontologies[[ont]]
    pathway_gene_sets <- split(gene_set$gene_short_name, gene_set$gs_name)
    if (!is.null(ancestor_specific_pathways) && nrow(ancestor_specific_pathways) > 0) {
      universe_genes <- unique(unlist(pathway_gene_sets[names(pathway_gene_sets) %in% ancestor_specific_pathways$gs_name]))
      if (length(universe_genes) > 0) {
        fora_res <- suppressWarnings(
          fgsea::fora(
            genes = deg_genes,
            pathways = pathway_gene_sets[names(pathway_gene_sets) %in% ancestor_specific_pathways$gs_name],
            universe = universe_genes
          )
        )
        fora_res_list[[ont]] <- fora_res %>%
          arrange(pval) %>%
          mutate(ontology = ont)
      }
    }
  }
  fora_res <- bind_rows(fora_res_list)

  # 6. Target genes expressed in this cell type
  expressed_targets <- character(0)
  if (!is.null(target_gene_expression)) {
    ancestor_expr <- target_gene_expression %>%
      filter(cell_group == ct)
    if (nrow(ancestor_expr) > 0) {
      expressed_targets <- unique(ancestor_expr$gene_short_name)
    }
  }
  expressed_targets_line <- if (length(expressed_targets) > 0) {
    paste(expressed_targets, collapse = ", ")
  } else {
    "none"
  }

  # 7. Pathway summary lines (top N by adjusted p-value, FDR-significant only).
  # Was raw `pval`; switched to `padj` so pathways shown to the LLM are
  # multiple-testing-corrected. Empirically this drops ~85-95% of the prior
  # candidate set (most fora hits in this run had padj ~ 0.99 despite raw
  # pval < 0.05), cutting noise the LLM was previously grasping at.
  pathway_lines <- ""
  pathway_genes_in_pathways <- character(0)
  if (!is.null(fora_res) && nrow(fora_res) > 0) {
    sig_pathways <- fora_res %>%
      filter(padj < sig_p_val_thresh) %>%
      arrange(padj) %>%
      group_by(ontology) %>%
      slice_head(n = top_n_pathways) %>%
      ungroup()
    if (nrow(sig_pathways) > 0) {
      pathway_lines <- paste(
        apply(sig_pathways, 1, function(row) {
          pathway_name <- humanize_pathway_name(row[["pathway"]])
          ontology <- row[["ontology"]]
          padj <- signif(as.numeric(row[["padj"]]), 2)
          # Genes of interest in this pathway
          pathway_genes <- if ("overlapGenes" %in% names(row)) row[["overlapGenes"]] else character(0)
          goi_in_pathway <- intersect(pathway_genes, degs_of_interest$gene_short_name)
          # Get arrows for these genes
          goi_arrows <- ""
          if (length(goi_in_pathway) > 0) {
            goi_arrows <- format_pathway_regulatory_genes(goi_in_pathway, degs_of_interest)
          }
          # Track all genes mentioned in pathways
          pathway_genes_in_pathways <<- union(pathway_genes_in_pathways, goi_in_pathway)
          paste0(
            "- [", ontology, "] ", pathway_name,
            " (padj=", padj,
            if (goi_arrows != "") paste0(", regulatory genes: ", goi_arrows) else "",
            ")"
          )
        }),
        collapse = "\n"
      )
      pathway_lines <- paste0("Top perturbed pathways:\n", pathway_lines)
    } else {
      pathway_lines <- "Top perturbed pathways: none"
    }
  } else {
    pathway_lines <- "Top perturbed pathways: none"
  }

  # 8. Other regulatory genes (not in any enriched pathway)
  other_regulatory_genes <- setdiff(degs_of_interest$gene_short_name, pathway_genes_in_pathways)
  other_regulatory_genes_line <- ""
  if (length(other_regulatory_genes) > 0) {
    other_regulatory_genes_line <- paste0(
      "Other regulatory genes: ",
      format_pathway_regulatory_genes(other_regulatory_genes, degs_of_interest)
    )
  }

  # 9. Build summary for this cell type
  perturbed_genes <- if (!is.null(target_gene_expression) && "gene_short_name" %in% colnames(target_gene_expression)) {
    unique(target_gene_expression$gene_short_name)
  } else {
    character(0)
  }
  perturbed_genes_line <- if (length(perturbed_genes) > 0) {
    paste0("Perturbed genes in this experiment: ", paste(perturbed_genes, collapse = ", "), ".")
  } else {
    "Perturbed genes in this experiment: [not specified]."
  }

  summary_text <- paste(
    perturbed_genes_line,
    paste0("Perturbation: ", perturbation_description),
    paste0("Cell type: ", ct),
    paste0("Abundance change: ", abundance_summary),
    paste0("Abundance severity: ", abundance_severity),
    paste0("Fitness label: ", fitness_label),
    paste0("Fitness severity: ", fitness_severity),
    paste0("Fitness evidence: ", fitness_evidence),
    paste0("Identity label: ", identity_label), # <-- NEW LINE
    paste0("Identity evidence: ", identity_evidence), # <-- NEW LINE
    paste0("Target genes expressed in this cell type: ", expressed_targets_line),
    paste0("Pathway background: ", pathway_background),
    pathway_lines,
    if (other_regulatory_genes_line != "") other_regulatory_genes_line,
    sep = "\n"
  )

  list(
    summary = summary_text,
    abundance = abundance_summary,
    abundance_code = abundance_code,
    abundance_severity = abundance_severity,
    fitness_label = fitness_label,
    fitness_severity = fitness_severity,
    fitness_evidence = fitness_evidence,
    identity_label = identity_label, # <-- NEW FIELD
    identity_evidence = identity_evidence, # <-- NEW FIELD
    present_above_thresh = cell_present_above_thresh,
    degs = cell_type_degs,
    goi_line = goi$goi_line,
    # Whether this cell has any REGULATOR (TF) genes of interest. goi_line now
    # also includes non-regulator DEGs, so it can no longer be used to gate the
    # LLM call (that would widen which cells get called); gate on this instead.
    has_regulator_goi = nrow(degs_of_interest) > 0,
    # Store only FDR-significant, non-empty-overlap pathways (same bar the LLM
    # prompt uses) instead of the raw fora dump -- the unfiltered result is
    # ~90%+ padj~1 / zero-overlap padding that bloats the impact table and any
    # downstream consumer without adding signal.
    pathways = if (!is.null(fora_res) && nrow(fora_res) > 0) {
      fora_res %>%
        dplyr::filter(!is.na(padj), padj < sig_p_val_thresh, overlap > 0) %>%
        dplyr::arrange(padj)
    } else {
      fora_res
    }
  )
}

build_lineage_context <- function(ct, parents, results, all_types = NULL) {
  # Only include parent summaries for parents in all_types (if provided)
  if (!is.null(all_types)) {
    parents <- intersect(parents, all_types)
  }

  # Build parent summaries. Pass the ancestor's LLM narrative (where nailed-down
  # literature and citations live), plus its genes-of-interest and enriched
  # pathways, so the descendant's call can reason about how an upstream defect
  # shapes this cell's own changes. This is CONTEXT only -- the descendant must
  # still name only its own <allowed_genes> in its output.
  parent_summaries <- unlist(lapply(parents, function(p) {
    pr <- results[[p]]
    if (is.null(pr)) return(NULL)
    # The distilled running summary is the bounded carrier of ancestral mechanism:
    # it is re-distilled at every called node and passed through unchanged by
    # empty nodes, so the salient root->...->here thread reaches this cell without
    # dumping the full ancestral chain into the prompt.
    narrative <- pr$lineage_context_summary
    if (is.null(narrative) || is.na(narrative)) narrative <- pr$llm_summary
    if (is.null(narrative) || is.na(narrative)) narrative <- pr$summary
    # Carry the IMMEDIATE parent's genes-of-interest across this one edge so the
    # model can connect an upstream regulator (down in the parent) to its target
    # (down in this cell). One hop only -- not accumulated up the whole lineage.
    goi_line <- if (!is.null(pr$goi_line) && !is.na(pr$goi_line) && nzchar(pr$goi_line) && pr$goi_line != "none") {
      paste0("\n  Immediate-parent genes of interest (look for upstream-regulator -> this-cell-target links): ", pr$goi_line)
    } else ""
    paste0("Ancestor (", p, ") [upstream context, do NOT report as this cell's own]:\n",
           narrative, goi_line)
  }), use.names = FALSE)

  # Build current cell type summary
  ct_context <- paste0(
    "Current Cell Type (", ct, "):\n",
    if (!is.null(results[[ct]]$llm_summary)) results[[ct]]$llm_summary else results[[ct]]$summary
  )

  # Combine parent summaries and current cell type summary
  context_text <- paste(
    c(parent_summaries, ct_context),
    collapse = "\n\n"
  )

  context_text
}

# Helper: Convert a list of Python DisruptedPathway objects to a tibble
py_disrupted_pathways_to_tibble <- function(x) {
  if (inherits(x, "data.frame")) {
    return(tibble::as_tibble(x))
  }
  if (is.null(x) || length(x) == 0) {
    return(tibble::tibble(
      name = character(),
      description = character(),
      dysregulated_genes = character()
    ))
  }
  # If reticulate hands back a python object (not yet converted), coerce again
  if (inherits(x, "python.builtin.object")) {
    x <- reticulate::py_to_r(x)
  }
  # If a data.frame/tibble already, just return it
  if (is.data.frame(x)) {
    return(tibble::as_tibble(x))
  }
  # If x is a character vector, not a list of lists
  if (is.atomic(x) && !is.list(x)) {
    return(tibble::tibble(
      name = as.character(x),
      description = NA_character_,
      dysregulated_genes = NA_character_
    ))
  }
  # If x is a single object with named fields, wrap it in a list
  if (!is.list(x)) {
    x <- list(x)
  }
  # Normalize any python objects inside the list
  x <- purrr::map(x, function(el) {
    if (inherits(el, "python.builtin.object")) reticulate::py_to_r(el) else el
  })
  # If x is a list of lists (the expected case)
  tryCatch(
    {
      tibble::tibble(
        name = purrr::map_chr(x, ~ .x$name %||% NA_character_),
        description = purrr::map_chr(x, ~ .x$description %||% NA_character_),
        dysregulated_genes = purrr::map_chr(x, ~ {
          dg <- .x$dysregulated_genes
          if (is.null(dg) || length(dg) == 0) return(NA_character_)
          if (is.list(dg)) dg <- unlist(dg)
          dg <- as.character(dg)
          dg <- dg[!is.na(dg) & nzchar(dg)]
          if (length(dg) == 0) NA_character_ else paste(dg, collapse = ", ")
        })
      ) %>%
        # Drop pathways the model named but could not ground in any gene from
        # <allowed_genes>. An ungrounded pathway carries no cell-specific
        # evidence; empirically these are the model anchoring on the
        # perturbation's canonical function (e.g. naming a myogenic pathway in a
        # non-muscle cell) rather than reading this cell's data. Demote-only.
        dplyr::filter(!is.na(.data$dysregulated_genes))
    },
    error = function(e) {
      preview <- paste(utils::capture.output(str(x, max.level = 2)), collapse = " ")
      message("LLM disrupted_pathways parse error; structure preview: ", preview)
      stop(e)
    }
  )
}

# --- Magnitude-arrow re-encoding of LLM-emitted gene lists -------------------
# The LLM returns gene tokens with a SINGLE direction arrow (copied from the
# <allowed_genes> block, e.g. "myod1 ↓"). To show fold-change MAGNITUDE in the
# rendered impact table -- in BOTH the report (goi_line) and the portal
# (dysregulated_genes) -- without adding per-surface logic, we rewrite those
# tokens once here, upstream: look each gene's shrunken log2FC up in this cell's
# DEG table, swap in magnitude arrows (see magnitude_arrow), and sort by |LFC|.

# Named vectors mapping lowercased gene symbol -> shrunken LFC and dysreg_type.
build_deg_arrow_lookup <- function(degs) {
  empty <- list(
    lfc = stats::setNames(numeric(0), character(0)),
    dir = stats::setNames(character(0), character(0))
  )
  if (is.null(degs) || !is.data.frame(degs) || nrow(degs) == 0 ||
    !"gene_short_name" %in% colnames(degs)) {
    return(empty)
  }
  d <- degs %>%
    dplyr::filter(!is.na(gene_short_name), nzchar(gene_short_name)) %>%
    dplyr::distinct(gene_short_name, .keep_all = TRUE)
  key <- tolower(d$gene_short_name)
  lfc <- if ("perturb_to_ctrl_shrunken_lfc" %in% colnames(d)) {
    as.numeric(d$perturb_to_ctrl_shrunken_lfc)
  } else {
    rep(NA_real_, nrow(d))
  }
  dir <- if ("dysreg_type" %in% colnames(d)) {
    as.character(d$dysreg_type)
  } else {
    rep(NA_character_, nrow(d))
  }
  list(lfc = stats::setNames(lfc, key), dir = stats::setNames(dir, key))
}

# Rewrite gene tokens (a character vector and/or comma-joined strings) with
# magnitude arrows, sorted by |LFC| descending (unmatched/NA last, alphabetical).
reencode_gene_arrows <- function(genes, lookup) {
  if (is.null(genes) || length(genes) == 0) {
    return(genes)
  }
  toks <- unlist(strsplit(as.character(genes), ","))
  toks <- trimws(toks)
  toks <- toks[!is.na(toks) & nzchar(toks)]
  if (length(toks) == 0) {
    return(character(0))
  }
  name <- trimws(sub("[↑↓]+$", "", toks))
  key <- tolower(name)
  lfc <- unname(lookup$lfc[key])
  dir <- unname(lookup$dir[key])
  out <- paste0(name, magnitude_arrow(lfc, dir))
  out[order(-abs(lfc), name)]
}

# Re-encode the dysregulated_genes string column of a disrupted-pathways tibble.
reencode_pathway_arrows <- function(tbl, lookup) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0 ||
    !"dysregulated_genes" %in% names(tbl)) {
    return(tbl)
  }
  tbl$dysregulated_genes <- vapply(
    tbl$dysregulated_genes,
    function(s) {
      if (is.na(s)) {
        return(NA_character_)
      }
      paste(reencode_gene_arrows(s, lookup), collapse = ", ")
    },
    character(1),
    USE.NAMES = FALSE
  )
  tbl
}

summarize_impact_in_lineage_context <- function(
  perturbation_description,
  target_gene_expression,
  dact_results,
  degs,
  ref_expression,
  ontologies,
  combined_psg,
  ai_notes_path,
  sig_p_val_thresh = 0.05,
  genes_of_interest = NULL,
  empirical_p_thresh = 1.0,
  llm_fun = NULL,
  max_lineage_depth = Inf,
  cell_types = NULL,
  primary_impact_summary = "",
  excluded_cell_types = NULL,
  abundance_phenotypes = NULL,
  fitness_phenotypes = NULL,
  identity_phenotypes = NULL,
  pre_cited_gene_claims = NULL,
  verbose = FALSE,
  ...
) {
  g <- if (class(combined_psg) == "cell_state_graph") combined_psg@graph else combined_psg

  # Determine which cell types to analyze
  if (is.null(cell_types)) {
    expressing_cell_types <- unique(target_gene_expression$cell_group)
    all_types <- unique(unlist(lapply(expressing_cell_types, function(ct) {
      if (ct %in% igraph::V(g)$name) {
        c(ct, get_descendants(ct, combined_psg))
      } else {
        ct
      }
    })))
  } else {
    all_types <- unique(unlist(lapply(cell_types, function(ct) {
      if (ct %in% igraph::V(g)$name) {
        c(ct, get_descendants(ct, combined_psg))
      } else {
        ct
      }
    })))
  }
  # Exclude specified cell types
  if (!is.null(excluded_cell_types)) {
    all_types <- setdiff(all_types, excluded_cell_types)
  }

  remaining <- all_types
  results <- list()
  processed <- character(0)
  max_iter <- length(all_types) * 2 # Prevent infinite loop
  iter <- 0

  if (!verbose) {
    pb <- txtProgressBar(min = 0, max = length(all_types), style = 3)
  }

  while (length(remaining) > 0 && iter < max_iter) {
    iter <- iter + 1
    progress <- FALSE

    ready <- purrr::keep(remaining, function(ct) {
      parents <- get_parents(combined_psg, ct)
      relevant_parents <- intersect(parents, all_types)
      all(relevant_parents %in% processed) || length(relevant_parents) == 0
    })
    if (length(ready) == 0) {
      ready <- remaining
    }
    for (ct in ready) {
      parents <- get_parents(combined_psg, ct)
      if (verbose) message(sprintf("[DEBUG] Processing cell type: %s (iteration %d)", ct, iter))
      results[[ct]] <- summarize_cell_type_impact(
        ct,
        perturbation_description,
        target_gene_expression,
        dact_results,
        degs,
        ref_expression,
        ontologies,
        ai_notes_path,
        sig_p_val_thresh,
        genes_of_interest,
        empirical_p_thresh = empirical_p_thresh,
        abundance_phenotypes = abundance_phenotypes,
        fitness_phenotypes = fitness_phenotypes,
        identity_phenotypes = identity_phenotypes
      )
      cell_impact_text <- build_lineage_context(ct, parents, results, all_types)
      if (!is.null(pre_cited_gene_claims) && nzchar(pre_cited_gene_claims)) {
        cell_impact_text <- paste0(
          "<pre_cited_claims>\n",
          pre_cited_gene_claims,
          "\n</pre_cited_claims>\n\n",
          cell_impact_text
        )
      }

      expresses_target <- !is.null(target_gene_expression) && ct %in% target_gene_expression$cell_group
      has_abundance <- !is.null(abundance_phenotypes) &&
        ct %in% abundance_phenotypes$cell_group &&
        !all(abundance_phenotypes %>% filter(cell_group == ct) %>% pull(abundance_code) %in% c("A0 No change"))

      has_fitness <- !is.null(fitness_phenotypes) &&
        ct %in% fitness_phenotypes$cell_group &&
        !all(fitness_phenotypes %>% filter(cell_group == ct) %>% pull(fitness_label) %in% c("F0 Normal", "F0 No significant phenotype"))

      has_identity <- !is.null(identity_phenotypes) &&
        ct %in% identity_phenotypes$cell_group &&
        !all(identity_phenotypes %>% filter(cell_group == ct) %>% pull(identity_label) %in% c("I0 Identity intact"))

      # A cell has something to explain if it shows its OWN transcriptional signal,
      # even without an abundance/fitness/identity phenotype label: its own
      # genes-of-interest (a key regulator moving is worth explaining regardless of
      # GO enrichment) or its own enriched pathways. We do NOT recycle an ancestor's
      # pathways as this cell's -- the ancestor is passed only as context (below).
      # Gate on regulator (TF) genes of interest only. goi_line was widened to
      # include non-regulator DEGs, so keying off it would call the LLM for cells
      # that previously fell through to the canned no-signal path. has_regulator_goi
      # preserves the original semantics (goi_line != "none" over regulators only).
      has_own_goi <- isTRUE(results[[ct]]$has_regulator_goi)
      has_own_enrichment <- !is.null(results[[ct]]$pathways) &&
        is.data.frame(results[[ct]]$pathways) && nrow(results[[ct]]$pathways) > 0

      # Distilled ancestral summary this cell forwards to its children. For a
      # called or excluded cell it is the cell's own concise_summary (set below);
      # a no-signal cell instead passes its parents' distilled summary straight
      # through, so the running narrative survives empty intermediate nodes.
      lineage_ctx_passthrough <- NULL

      if (!is.null(excluded_cell_types) && ct %in% excluded_cell_types) {
        if (verbose) message(sprintf("[DEBUG] Cell type %s is excluded.", ct))
        concise_summary <- NA_character_
        disrupted_pathways <- NULL
        other_dysregulated_genes <- NULL
      } else if (isTRUE(results[[ct]]$present_above_thresh) &&
                 (has_abundance || has_fitness || has_identity || has_own_goi || has_own_enrichment)) {
        if (verbose) message(sprintf("[DEBUG] Calling LLM for cell type: %s (expresses_target: %s, has_abundance: %s, has_fitness: %s, has_identity: %s, has_own_goi: %s, has_own_enrichment: %s)", ct, expresses_target, has_abundance, has_fitness, has_identity, has_own_goi, has_own_enrichment))
        allowed_degs_map <- build_allowed_degs_map(results[[ct]]$degs)
        if (verbose) {
          message(sprintf("[DEBUG] allowed_degs for %s: %d genes", ct, length(allowed_degs_map)))
        }
        if (length(allowed_degs_map) == 0) {
          # No callable genes -> skip the LLM and emit a short, factual line with
          # NO mechanistic interpretation (Amy's feedback: these near-loss /
          # no-gene cells were getting over-interpreted paragraphs). The cell's
          # phenotype is still shown in its impact-table row.
          if (verbose) message(sprintf("[DEBUG] No allowed genes for %s; canned no-gene summary (LLM skipped).", ct))
          concise_summary <- paste0(
            "No genes passed the significance threshold in this cell type, so no ",
            "cell-type-specific molecular mechanism is inferred here."
          )
          disrupted_pathways <- NULL
          other_dysregulated_genes <- NULL
        } else {
        llm_structured <- tryCatch(
          {
            if (!is.null(llm_fun)) {
              llm_res <- reticulate::py_to_r(
                llm_fun(
                  cell_type = ct,
                  cell_impact_text = cell_impact_text,
                  primary_effect_summary = primary_impact_summary,
                  allowed_degs = allowed_degs_map,
                  ...
                )
              )
              if (verbose) {
                message(sprintf("[DEBUG] LLM structured payload for %s:", ct))
                utils::str(llm_res, max.level = 2)
              }
              llm_res
            } else {
              NULL
            }
          },
          error = function(e) {
            if (verbose) message(sprintf("[ERROR] LLM call failed for cell type '%s': %s", ct, e$message))
            NULL
          }
        )
        if (verbose && !is.null(llm_structured)) {
          message(sprintf("[DEBUG] LLM structured keys for %s: %s", ct, paste(names(llm_structured), collapse = ", ")))
          dp <- llm_structured$disrupted_pathways
          message(sprintf("[DEBUG] LLM disrupted_pathways for %s: %s", ct, paste(utils::capture.output(str(dp, max.level = 2)), collapse = " ")))
        }
        concise_summary <- if (!is.null(llm_structured) && !is.null(llm_structured$concise_summary)) llm_structured$concise_summary else NA_character_
        disrupted_pathways <- if (!is.null(llm_structured)) llm_structured$disrupted_pathways else NULL
        other_dysregulated_genes <- if (!is.null(llm_structured)) llm_structured$other_dysregulated_genes else NULL
        # Defensive: ensure it's a character vector or NULL
        if (!is.null(other_dysregulated_genes) && !is.character(other_dysregulated_genes)) {
          other_dysregulated_genes <- as.character(other_dysregulated_genes)
        }
        # Re-encode single-arrow LLM gene tokens into magnitude arrows (1/2/3)
        # using this cell's DEG LFCs, so the report and portal render effect size
        # directly. disrupted_pathways is converted to a tibble here (idempotent
        # with the assembly-time conversion below).
        deg_lookup <- build_deg_arrow_lookup(results[[ct]]$degs)
        other_dysregulated_genes <- reencode_gene_arrows(other_dysregulated_genes, deg_lookup)
        disrupted_pathways <- tryCatch(
          reencode_pathway_arrows(py_disrupted_pathways_to_tibble(disrupted_pathways), deg_lookup),
          error = function(e) {
            if (verbose) message(sprintf("[WARN] arrow re-encode failed for %s: %s", ct, e$message))
            disrupted_pathways
          }
        )
        }
      } else {
        if (verbose) message(sprintf("[DEBUG] Skipping LLM for cell type: %s (no abundance, fitness, or identity phenotype)", ct))
        parent_phenotype_summaries <- sapply(parents, function(p) {
          parent_info <- results[[p]]
          parent_has_phenotype <- !is.null(parent_info$abundance) && !is.na(parent_info$abundance) && parent_info$abundance != "no_change"
          parent_has_fitness <- !is.null(parent_info$fitness_label) && !is.na(parent_info$fitness_label) && parent_info$fitness_label != "F0 Normal"
          parent_has_identity <- !is.null(parent_info$identity_label) && !is.na(parent_info$identity_label) && parent_info$identity_label != "I0 Identity intact"
          if (parent_has_phenotype || parent_has_fitness || parent_has_identity) {
            paste0(
              p, " of ", ct, " has phenotype: ",
              paste(
                c(
                  if (parent_has_phenotype) paste0("abundance: ", humanize_pheno_label(parent_info$abundance)),
                  if (parent_has_fitness) paste0("fitness: ", humanize_pheno_label(parent_info$fitness_label)),
                  if (parent_has_identity) paste0("identity: ", humanize_pheno_label(parent_info$identity_label))
                ),
                collapse = "; "
              )
            )
          } else {
            NULL
          }
        })
        parent_phenotype_summaries <- parent_phenotype_summaries[!is.na(parent_phenotype_summaries) & parent_phenotype_summaries != "" & parent_phenotype_summaries != "NULL"]

        # This branch is only reached when the cell has NO own signal of any kind
        # (no phenotype, no goi, no enrichment) -- there is nothing cell-autonomous
        # to explain. We record a pointer to any upstream (ancestral) phenotype as
        # context, but we do NOT fabricate a mechanism for this cell by copying the
        # ancestor's pathways/genes onto it -- that would pass off an ancestral
        # problem as a descendant problem. Its pathways/genes stay empty.
        if (length(parent_phenotype_summaries) == 0) {
          if (verbose) message(sprintf("[DEBUG] Cell type %s and its parents have no detectable phenotypes.", ct))
          concise_summary <- paste0(
            "No cell-autonomous signal for '", ct, "': no phenotype, genes-of-interest, or pathway enrichment in this cell type or its ancestors."
          )
        } else {
          if (verbose) message(sprintf("[DEBUG] Cell type %s: no own signal; noting upstream ancestral phenotype as context only (no pathways recycled).", ct))
          concise_summary <- paste0(
            "No cell-autonomous transcriptional signal for '", ct,
            "'. Downstream of ancestral phenotype(s): ",
            paste(parent_phenotype_summaries, collapse = " "),
            " See the ancestor(s) for the upstream mechanism."
          )
        }

        disrupted_pathways <- NULL
        other_dysregulated_genes <- NULL

        # Pass the parents' distilled running summary through unchanged (no
        # accretion), so this no-signal cell does not break the top-down summary
        # chain for its own descendants. Its displayed row still says "no own
        # change" (concise_summary above); this is only what it forwards.
        ups <- vapply(parents, function(p) {
          v <- results[[p]]$lineage_context_summary
          if (is.null(v) || is.na(v)) v <- results[[p]]$llm_summary
          if (is.null(v) || is.na(v)) NA_character_ else as.character(v)
        }, character(1))
        ups <- ups[!is.na(ups) & nzchar(ups)]
        lineage_ctx_passthrough <- if (length(ups) == 0) NA_character_ else paste(ups, collapse = "\n\n")
      }

      results[[ct]]$context <- cell_impact_text
      results[[ct]]$llm_summary <- ifelse(is.null(concise_summary) || concise_summary == "NULL", NA_character_, concise_summary)
      # Forwarded (distilled) ancestral summary: parents' pass-through for a
      # no-signal cell, else this cell's own concise_summary.
      results[[ct]]$lineage_context_summary <- if (!is.null(lineage_ctx_passthrough)) lineage_ctx_passthrough else results[[ct]]$llm_summary
      results[[ct]]$llm_disrupted_pathways <- disrupted_pathways
      results[[ct]]$llm_other_dysregulated_genes <- other_dysregulated_genes
      processed <- c(processed, ct)
      progress <- TRUE

      # Update progress bar or debug message for each cell processed
      if (!verbose) {
        setTxtProgressBar(pb, length(processed) + 1)
      } else {
        message(sprintf("[DEBUG] %d cell types remaining to process.", length(all_types) - length(processed)))
      }
    }
    remaining <- setdiff(remaining, ready)
    if (!progress) break
  }

  if (!verbose) close(pb)

  results <- tibble(
    cell_type = names(results),
    data = unname(results)
  ) %>% unnest_wider(data)

  # Drop cell types the abundance analysis judged not reliably present
  # (present_above_thresh = FALSE): this perturbation experiment cannot capture
  # them, so they get NO impact-table entry -- their large DEG sets are noise/
  # indirect on few cells, not attributable biology. They were still available as
  # context during processing, so any present descendants inherited the distilled
  # ancestral narrative through them.
  if ("present_above_thresh" %in% names(results)) {
    n_before <- nrow(results)
    results <- results %>% filter(present_above_thresh %in% TRUE)
    if (verbose && n_before > nrow(results)) {
      message(sprintf("[DEBUG] Dropped %d not-present-above-threshold cell types from the impact table.", n_before - nrow(results)))
    }
    results <- results %>% select(-present_above_thresh)
  }

  # Ensure llm_disrupted_pathways column exists before mutate
  if (!"llm_disrupted_pathways" %in% names(results)) {
    results$llm_disrupted_pathways <- vector("list", nrow(results))
  }

  results <- results %>%
    mutate(
      llm_disrupted_pathways = purrr::map(llm_disrupted_pathways, function(x) {
        tryCatch(
          py_disrupted_pathways_to_tibble(x),
          error = function(e) {
            warning(sprintf("Error processing llm_disrupted_pathways: %s", e$message))
            tibble::tibble(
              name = character(),
              description = character(),
              dysregulated_genes = character()
            )
          }
        )
      })
    )

  results
}

collect_cell_loss_explanations <- function(explanations_df) {
  explanations_df %>%
    rowwise() %>%
    mutate(
      pathway_explanations = {
        if (is.null(llm_disrupted_pathways) || !is.data.frame(llm_disrupted_pathways) || nrow(llm_disrupted_pathways) == 0) {
          "No disrupted pathways."
        } else {
          paste(
            apply(llm_disrupted_pathways, 1, function(pathway_row) {
              glue::glue(
                "- Name: {pathway_row[['name']]}\n  Description: {pathway_row[['description']]}\n  Dysregulated genes: {paste(as.character(pathway_row[['dysregulated_genes']]), collapse=', ')}"
              )
            }),
            collapse = "\n"
          )
        }
      },
      other_genes_explanation = if (is.null(llm_other_dysregulated_genes)) {
        "Other dysregulated genes: None"
      } else {
        paste(
          "Other dysregulated genes:",
          paste(as.character(llm_other_dysregulated_genes), collapse = ", ")
        )
      },
      explanation = paste(
        "Cell type:", cell_type,
        "\nSummary:", llm_summary,
        "\nDisrupted pathways:\n", pathway_explanations,
        "---------------------\n"
      )
    ) %>%
    ungroup() %>%
    pull(explanation)
}

# Helper to get cell type link
cell_type_link <- function(cell_type) {
  post_file <- cell_type_post_tbl %>%
    filter(cell_type == !!cell_type) %>%
    pull(post_file)
  if (length(post_file) > 0 && !is.na(post_file)) {
    gt::html(as.character(htmltools::a(href = post_link(post_file, prefix = "^content", base_url = base_url), cell_type)))
  } else {
    cell_type
  }
}

# Helper to get gene link
gene_link <- function(gene) {
  post_file <- gene_post_tbl %>%
    filter(gene == !!gene) %>%
    pull(post_file)
  if (length(post_file) > 0 && !is.na(post_file)) {
    gt::html(as.character(htmltools::a(href = post_link(post_file, prefix = "^content", base_url = base_url), gene)))
  } else {
    gene
  }
}


zscape_gt_perturbation_impact_table <- function(
  impact_table,
  cell_type_post_tbl,
  gene_post_tbl,
  cell_types_of_interest = NULL,
  base_url = "",
  show_only_with_pathways = TRUE
) {
  # Helper to get cell type link
  cell_type_link <- function(cell_type) {
    post_file <- cell_type_post_tbl %>%
      filter(cell_type == !!cell_type) %>%
      pull(post_file)
    if (length(post_file) > 0 && !is.na(post_file)) {
      gt::html(as.character(htmltools::a(href = post_link(post_file, prefix = "^content", base_url = base_url), cell_type)))
    } else {
      cell_type
    }
  }

  # Helper to get gene link
  gene_link <- function(gene) {
    post_file <- gene_post_tbl %>%
      filter(gene == !!gene) %>%
      pull(post_file)
    if (length(post_file) > 0 && !is.na(post_file)) {
      gt::html(as.character(htmltools::a(href = post_link(post_file, prefix = "^content", base_url = base_url), gene)))
    } else {
      gene
    }
  }

  # Optionally filter to cell types of interest
  if (!is.null(cell_types_of_interest)) {
    impact_table <- impact_table %>% filter(cell_type %in% cell_types_of_interest)
  }

  # Only proceed if there are valid pathway explanations
  if (!is.null(impact_table) && nrow(impact_table) > 0 &&
    "llm_disrupted_pathways" %in% colnames(impact_table) &&
    any(!purrr::map_lgl(impact_table$llm_disrupted_pathways, is.null))) {
    pathway_table <- impact_table %>%
      group_by(cell_type) %>%
      mutate(
        phenotype_summary = paste(
          na.omit(
            c(
              if (!is.na(abundance_code) && !(abundance_code %in% c("A0 No change"))) abundance_code else NULL,
              if (!is.na(fitness_label) && !(fitness_label %in% c("F0 No significant phenotype", "F0 Normal"))) fitness_label else NULL,
              if (!is.na(identity_label) && !(identity_label %in% c("I0 Identity intact"))) identity_label else NULL
            )
          ),
          collapse = "; "
        ),
        effect_type = paste(effect_type, phenotype_summary, sep = " - ")
      ) %>%
      ungroup() %>%
      dplyr::filter(!purrr::map_lgl(impact_table$llm_disrupted_pathways, is.null)) %>%
      tidyr::unnest(llm_disrupted_pathways) %>%
      select(
        cell_type,
        phenotype_summary,
        effect_type,
        pathway = name,
        description,
        dysregulated_genes
      ) %>%
      distinct() %>%
      ungroup() %>%
      arrange(effect_type, cell_type, pathway)

    pathway_table_linked <- pathway_table %>%
      arrange(cell_type, effect_type, pathway) %>%
      group_by(cell_type) %>%
      mutate(
        cell_type_row = row_number(),
        cell_type_display = ifelse(
          cell_type_row == 1,
          paste0(
            purrr::map_chr(cell_type, cell_type_link),
            "<br><span style='font-size:smaller;color:gray;'>", effect_type, "</span>"
          ),
          as.character(cell_type)
        ),
        dysregulated_genes = purrr::map(
          dysregulated_genes,
          ~ purrr::map(.x, gene_link)
        )
      ) %>%
      ungroup() %>%
      mutate(
        dysregulated_genes = purrr::map_chr(
          dysregulated_genes,
          ~ paste(purrr::map_chr(.x, as.character), collapse = ", ")
        )
      )

    # Render the gt table with markdown links
    pathway_table_linked %>%
      distinct() %>%
      arrange(cell_type, effect_type, pathway) %>%
      group_by(cell_type) %>%
      mutate(
        cell_type_display_row = row_number()
      ) %>%
      ungroup() %>%
      mutate(
        pathway = stringr::str_to_sentence(pathway),
        description = stringr::str_to_sentence(description),
        dysregulated_genes = stringr::str_replace_all(
          dysregulated_genes,
          "(↑)", "<span style='color:red;'>\\1</span>"
        ) %>%
          stringr::str_replace_all(
            "(↓)", "<span style='color:blue;'>\\1</span>"
          )
      ) %>%
      select(cell_type_display, pathway, description, dysregulated_genes, cell_type_display_row) %>%
      gt::gt() %>%
      gt::fmt_markdown(columns = c(cell_type_display, dysregulated_genes)) %>%
      gt::tab_style(
        style = gt::cell_text(color = "transparent"),
        locations = gt::cells_body(
          columns = "cell_type_display",
          rows = .data$cell_type_display_row != 1
        )
      ) %>%
      gt::tab_style(
        style = list(gt::cell_text(size = "smaller")),
        locations = gt::cells_body()
      ) %>%
      gt::tab_style(
        style = gt::cell_borders(sides = c("top", "bottom"), weight = px(0)),
        locations = gt::cells_body(
          rows = .data$cell_type_display_row != 1
        )
      ) %>%
      gt::cols_hide(columns = c("cell_type_display_row")) %>%
      gt::cols_width(
        cell_type_display ~ px(200),
        pathway ~ px(200),
        dysregulated_genes ~ px(200)
      ) %>%
      gt::cols_label(
        cell_type_display = "Cell type",
        pathway = "Pathway",
        description = "LLM explanation",
        dysregulated_genes = "Dysregulated genes"
      ) %>%
      gt::opt_interactive(use_search = TRUE, use_compact_mode = TRUE)
  } else {
    gt::gt(data.frame(Message = "No valid LLM pathways found in impact_table."))
  }
}


filter_and_pivot_dysregulated_genes <- function(
  impact_table,
  phenotype_types = c("abundance", "identity", "stress")
) {
  # Ensure columns exist and are list-columns
  if (!"llm_disrupted_pathways" %in% names(impact_table)) {
    impact_table$llm_disrupted_pathways <- vector("list", nrow(impact_table))
  }
  if (!"llm_other_dysregulated_genes" %in% names(impact_table)) {
    impact_table$llm_other_dysregulated_genes <- vector("list", nrow(impact_table))
  }

  filtered <- impact_table %>%
    filter(
      (
        ("abundance" %in% phenotype_types & !is.na(abundance_code) & abundance_code != "A0 No change") |
          ("identity" %in% phenotype_types & !is.na(identity_label) & identity_label != "I0 Identity intact") |
          ("fitness" %in% phenotype_types & !is.na(fitness_label) & !(fitness_label %in% c("F0 No significant phenotype", "F0 Normal"))) |
          ("stress" %in% phenotype_types & !is.na(fitness_label) & stringr::str_detect(fitness_label, regex("stress|F3", ignore_case = TRUE)))
      )
    )

  # Unnest pathways and preserve up/down status
  pathways_long <- filtered %>%
    select(cell_type, llm_disrupted_pathways) %>%
    tidyr::unnest(llm_disrupted_pathways, keep_empty = TRUE) %>%
    mutate(
      gene_list = if ("dysregulated_genes" %in% names(.)) stringr::str_split(dysregulated_genes, ",\\s*") else list(NA)
    ) %>%
    select(cell_type, pathway_name = name, pathway_description = description, gene_list) %>%
    tidyr::unnest(gene_list, keep_empty = TRUE) %>%
    mutate(
      gene = stringr::str_remove_all(gene_list, "[↑↓]"),
      direction = dplyr::case_when(
        stringr::str_detect(gene_list, "↑") ~ "up",
        stringr::str_detect(gene_list, "↓") ~ "down",
        TRUE ~ NA_character_
      )
    ) %>%
    select(cell_type, gene, direction, pathway_name, pathway_description) %>%
    dplyr::filter(!is.na(cell_type) & !is.na(gene))

  # Unnest other genes and preserve up/down status
  other_genes_long <- filtered %>%
    select(cell_type, llm_other_dysregulated_genes) %>%
    tidyr::unnest(llm_other_dysregulated_genes, keep_empty = TRUE) %>%
    mutate(
      gene = stringr::str_remove_all(llm_other_dysregulated_genes, "[↑↓]"),
      direction = dplyr::case_when(
        stringr::str_detect(llm_other_dysregulated_genes, "↑") ~ "up",
        stringr::str_detect(llm_other_dysregulated_genes, "↓") ~ "down",
        TRUE ~ NA_character_
      )
    ) %>%
    select(cell_type, gene, direction) %>%
    filter(!is.na(gene) & !is.na(cell_type))

  # Combine all gene/cell_type/pathway info
  all_genes_long <- dplyr::bind_rows(
    pathways_long %>% select(cell_type, gene, direction, pathway_name, pathway_description),
    other_genes_long %>% mutate(pathway_name = NA_character_, pathway_description = NA_character_)
  )

  # Summarize for each gene
  result <- all_genes_long %>%
    group_by(gene) %>%
    summarize(
      cell_types = list(unique(cell_type)),
      n_cell_types = length(unique(cell_type)),
      cell_types_down = list(unique(cell_type[direction == "down"])),
      n_cell_types_down = length(unique(cell_type[direction == "down"])),
      pathways = list(
        na.omit(unique(
          tibble::tibble(
            name = pathway_name,
            description = pathway_description
          )
        ))
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(n_cell_types_down))

  return(result)
}
