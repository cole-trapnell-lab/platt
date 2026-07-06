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

format_goi <- function(degs_of_interest) {
  if (
    is.null(degs_of_interest) ||
      !is.data.frame(degs_of_interest) ||
      nrow(degs_of_interest) == 0 ||
      !"gene_short_name" %in% colnames(degs_of_interest)
  ) {
    return(list(goi_vector_sorted = character(0), goi_line = "none"))
  }
  # Ensure columns exist
  if (!"dysreg_type" %in% colnames(degs_of_interest)) {
    degs_of_interest$dysreg_type <- NA_character_
  }
  if (!"perturb_to_ctrl_shrunken_lfc" %in% colnames(degs_of_interest)) {
    degs_of_interest$perturb_to_ctrl_shrunken_lfc <- NA_real_
  }
  goi_df <- degs_of_interest %>%
    mutate(
      arrow = dplyr::case_when(
        dysreg_type == "Overexpressed" ~ "\u2191",
        dysreg_type == "Underexpressed" ~ "\u2193",
        !is.na(perturb_to_ctrl_shrunken_lfc) & perturb_to_ctrl_shrunken_lfc > 0 ~ "\u2191",
        !is.na(perturb_to_ctrl_shrunken_lfc) & perturb_to_ctrl_shrunken_lfc < 0 ~ "\u2193",
        TRUE ~ ""
      ),
      lfc = perturb_to_ctrl_shrunken_lfc,
      gene_arrow = paste0(gene_short_name, arrow)
    )
  down <- goi_df %>%
    filter(arrow == "\u2193") %>%
    arrange(lfc) %>%
    pull(gene_arrow)
  up <- goi_df %>%
    filter(arrow == "\u2191") %>%
    arrange(desc(lfc)) %>%
    pull(gene_arrow)
  none <- goi_df %>%
    filter(arrow == "") %>%
    arrange(gene_short_name) %>%
    pull(gene_arrow)
  goi_vector_sorted <- c(down, up, none)
  goi_line <- paste(goi_vector_sorted, collapse = ", ")
  list(goi_vector_sorted = goi_vector_sorted, goi_line = goi_line)
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
      arrow = dplyr::case_when(
        dysreg_type == "Overexpressed" ~ "\u2191",
        dysreg_type == "Underexpressed" ~ "\u2193",
        !is.na(perturb_to_ctrl_shrunken_lfc) & perturb_to_ctrl_shrunken_lfc > 0 ~ "\u2191",
        !is.na(perturb_to_ctrl_shrunken_lfc) & perturb_to_ctrl_shrunken_lfc < 0 ~ "\u2193",
        TRUE ~ ""
      ),
      lfc = perturb_to_ctrl_shrunken_lfc,
      gene_arrow = paste0(gene_short_name, arrow)
    )
  down <- goi_df %>%
    filter(arrow == "\u2193") %>%
    arrange(lfc) %>%
    pull(gene_arrow)
  up <- goi_df %>%
    filter(arrow == "\u2191") %>%
    arrange(desc(lfc)) %>%
    pull(gene_arrow)
  none <- goi_df %>%
    filter(arrow == "") %>%
    arrange(gene_short_name) %>%
    pull(gene_arrow)
  goi_arrows_sorted <- c(down, up, none)
  paste(goi_arrows_sorted, collapse = ", ")
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
  empirical_fdr_thresh = 1.0,
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

  # 4. DEGs and key regulatory genes
  cell_type_degs <- degs %>%
    filter(cell_group == ct, perturb_to_ctrl_p_value < sig_p_val_thresh)
  degs_of_interest <- if (!is.null(genes_of_interest)) {
    cell_type_degs %>% filter(gene_short_name %in% genes_of_interest)
  } else {
    cell_type_degs
  }
  # Genes of interest are named, per-gene claims handed to the LLM (esp. TFs), so
  # gate them on the empirical-FDR control-sampling-artifact estimate. The FORA
  # foreground (deg_genes, below) is deliberately left on the relaxed nominal-p
  # set -- enrichment is a set-level test with its own pathway-level padj, and
  # gating individual genes would defeat it. NA (undecorated) rows are kept.
  if (!is.null(empirical_fdr_thresh) && empirical_fdr_thresh < 1 &&
      "empirical_fdr" %in% colnames(degs_of_interest)) {
    degs_of_interest <- degs_of_interest %>%
      filter(is.na(empirical_fdr) | empirical_fdr <= empirical_fdr_thresh)
  }
  goi <- format_goi(degs_of_interest)

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
    degs = cell_type_degs,
    goi_line = goi$goi_line,
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

  # Build parent summaries
  parent_summaries <- unlist(lapply(parents, function(p) {
    parent_context <- if (!is.null(results[[p]]$llm_summary)) {
      results[[p]]$llm_summary
    } else {
      results[[p]]$summary
    }
    paste0("Parent (", p, "):\n", parent_context)
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
          # Empty list / NULL — happens when the LLM correctly returns no genes
          # for a pathway (per the Phase-2 prompt). Treat as NA, not as an
          # error, so other pathways for this cell still survive.
          if (is.null(dg) || length(dg) == 0) return(NA_character_)
          if (is.list(dg)) dg <- unlist(dg)
          dg <- as.character(dg)
          dg <- dg[!is.na(dg) & nzchar(dg)]
          if (length(dg) == 0) NA_character_ else paste(dg, collapse = ", ")
        })
      )
    },
    error = function(e) {
      preview <- paste(utils::capture.output(str(x, max.level = 2)), collapse = " ")
      message("LLM disrupted_pathways parse error; structure preview: ", preview)
      stop(e)
    }
  )
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
  empirical_fdr_thresh = 1.0,
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
        empirical_fdr_thresh = empirical_fdr_thresh,
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

      if (!is.null(excluded_cell_types) && ct %in% excluded_cell_types) {
        if (verbose) message(sprintf("[DEBUG] Cell type %s is excluded.", ct))
        concise_summary <- NA_character_
        disrupted_pathways <- NULL
        other_dysregulated_genes <- NULL
      } else if (has_abundance || has_fitness || has_identity) {
        if (verbose) message(sprintf("[DEBUG] Calling LLM for cell type: %s (expresses_target: %s, has_abundance: %s, has_fitness: %s, has_identity: %s)", ct, expresses_target, has_abundance, has_fitness, has_identity))
        allowed_degs_map <- build_allowed_degs_map(results[[ct]]$degs)
        if (verbose) {
          message(sprintf("[DEBUG] allowed_degs for %s: %d genes", ct, length(allowed_degs_map)))
        }
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
                  if (parent_has_phenotype) paste0("abundance: ", parent_info$abundance),
                  if (parent_has_fitness) paste0("fitness: ", parent_info$fitness_label),
                  if (parent_has_identity) paste0("identity: ", parent_info$identity_label)
                ),
                collapse = "; "
              )
            )
          } else {
            NULL
          }
        })
        parent_phenotype_summaries <- parent_phenotype_summaries[!is.na(parent_phenotype_summaries) & parent_phenotype_summaries != "" & parent_phenotype_summaries != "NULL"]

        if (length(parent_phenotype_summaries) == 0) {
          if (verbose) message(sprintf("[DEBUG] Cell type %s and its parents have no detectable phenotypes.", ct))
          concise_summary <- paste0(
            "No LLM call for cell type '", ct, "': no detectable phenotypes in this cell type or its parents."
          )
        } else {
          if (verbose) message(sprintf("[DEBUG] Cell type %s: propagating parent phenotype info.", ct))
          concise_summary <- paste(parent_phenotype_summaries, collapse = "\n")
        }

        disrupted_pathways <- purrr::map(parents, ~ results[[.x]]$llm_disrupted_pathways) %>%
          purrr::compact() %>%
          purrr::flatten()
        other_dysregulated_genes <- purrr::map(parents, ~ results[[.x]]$llm_other_dysregulated_genes) %>%
          purrr::compact() %>%
          unlist()
      }

      results[[ct]]$context <- cell_impact_text
      results[[ct]]$llm_summary <- ifelse(is.null(concise_summary) || concise_summary == "NULL", NA_character_, concise_summary)
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
