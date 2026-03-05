#' Create an Empty Disrupted-Pathways Table
#'
#' Returns a zero-row tibble with the standardized columns used for storing
#' disrupted pathway summaries.
#'
#' @return A tibble with columns `name`, `description`, and
#'   `dysregulated_genes`.
empty_disrupted_pathways_tbl <- function() {
  tibble::tibble(
    name = character(),
    description = character(),
    dysregulated_genes = character()
  )
}

#' Normalize Overlap-Genes Input to a Character Vector
#'
#' Handles common overlap gene representations (`NULL`, vectors, nested lists)
#' and returns a flattened unique character vector.
#'
#' @param x Overlap gene input, typically from enrichment results.
#'
#' @return A character vector of gene identifiers.
extract_overlap_genes <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(character(0))
  }
  if (is.list(x)) {
    return(unique(unlist(x, use.names = FALSE)))
  }
  as.character(x)
}

#' Build Deterministic Pathway Summaries
#'
#' Converts pathway enrichment output into a compact deterministic table used in
#' no-LLM narrative outputs.
#'
#' @param pathways Data frame of pathway enrichment results. Expected columns
#'   include `pval`, `overlapGenes`, and optionally `ontology` and `pathway`.
#' @param sig_p_val_thresh Numeric p-value threshold used to retain pathways.
#' @param top_n_pathways Maximum number of pathways to keep per ontology (or in
#'   total when ontology is absent).
#'
#' @return A tibble with `name`, `description`, and `dysregulated_genes`.
build_deterministic_disrupted_pathways <- function(pathways, sig_p_val_thresh = 0.05, top_n_pathways = 5) {
  if (is.null(pathways) || !is.data.frame(pathways) || nrow(pathways) == 0) {
    return(empty_disrupted_pathways_tbl())
  }

  keep <- pathways %>%
    dplyr::filter(.data$pval < sig_p_val_thresh)

  if (nrow(keep) == 0) {
    return(empty_disrupted_pathways_tbl())
  }

  if ("ontology" %in% names(keep)) {
    keep <- keep %>%
      dplyr::group_by(.data$ontology) %>%
      dplyr::slice_head(n = top_n_pathways) %>%
      dplyr::ungroup()
  } else {
    keep <- keep %>%
      dplyr::slice_head(n = top_n_pathways)
  }

  keep %>%
    dplyr::mutate(
      name = if ("pathway" %in% names(keep)) {
        vapply(.data$pathway, humanize_pathway_name, FUN.VALUE = character(1))
      } else {
        NA_character_
      },
      description = paste0(
        if ("ontology" %in% names(keep)) paste0("[", .data$ontology, "] ") else "",
        "p=", signif(.data$pval, 2)
      ),
      dysregulated_genes = vapply(
        .data$overlapGenes,
        function(genes) paste(extract_overlap_genes(genes), collapse = ", "),
        FUN.VALUE = character(1)
      )
    ) %>%
    dplyr::select(.data$name, .data$description, .data$dysregulated_genes)
}

#' Summarize Cell-Type Impact Without LLM Notes
#'
#' Produces per-cell-type phenotype, DEG, and pathway summaries from
#' perturbation and ontology inputs.
#'
#' @param ct Character scalar cell type to summarize.
#' @param target_gene_expression Data frame of target-gene expression calls by
#'   cell group.
#' @param dact_results Differential abundance result table.
#' @param degs Differential expression table.
#' @param ontologies Named list of ontology gene-set tables with columns
#'   `gs_name` and `gene_short_name`.
#' @param sig_p_val_thresh Significance cutoff used for abundance/DEG/pathway
#'   calls.
#' @param genes_of_interest Optional character vector used to filter DEG output
#'   when generating `goi_line`.
#' @param power_thresh Power threshold for classifying non-significant abundance
#'   results as "no_change" vs underpowered.
#' @param top_n_pathways Maximum pathways retained per ontology in intermediate
#'   summaries.
#' @param abundance_phenotypes Optional precomputed abundance phenotype table.
#' @param fitness_phenotypes Optional precomputed fitness phenotype table.
#' @param identity_phenotypes Optional precomputed identity phenotype table.
#'
#' @return A named list containing abundance/phenotype fields, filtered DEGs,
#'   GOI summary text, and pathway enrichment results.
summarize_cell_type_impact_no_ai_notes <- function(
  ct,
  target_gene_expression,
  dact_results,
  degs,
  ontologies,
  sig_p_val_thresh = 0.05,
  genes_of_interest = NULL,
  power_thresh = 0.8,
  top_n_pathways = 5,
  abundance_phenotypes = NULL,
  fitness_phenotypes = NULL,
  identity_phenotypes = NULL
) {
  abundance_row <- if (!is.null(abundance_phenotypes)) {
    abundance_phenotypes %>% dplyr::filter(.data$cell_group == ct)
  } else {
    tibble::tibble()
  }
  abundance_code <- if (nrow(abundance_row) > 0) abundance_row$abundance_code[1] else NA_character_
  abundance_severity <- if (nrow(abundance_row) > 0) abundance_row$abundance_severity[1] else NA_character_

  dact_row <- dact_results %>% dplyr::filter(.data$cell_group == ct)
  abundance_summary <- if (!is.na(abundance_code)) {
    abundance_code
  } else if (nrow(dact_row) > 0 && all(c("delta_q_value", "delta_log_abund", "power") %in% colnames(dact_row))) {
    if (dact_row$delta_q_value[1] < sig_p_val_thresh) {
      if (dact_row$delta_log_abund[1] < 0) "depleted" else if (dact_row$delta_log_abund[1] > 0) "enriched" else "significant_no_direction"
    } else if (dact_row$power[1] >= power_thresh) {
      "no_change"
    } else {
      "unknown (underpowered)"
    }
  } else {
    NA_character_
  }

  fitness_rows <- if (!is.null(fitness_phenotypes)) {
    fitness_phenotypes %>% dplyr::filter(.data$cell_group == ct)
  } else {
    tibble::tibble()
  }
  fitness_label <- if (nrow(fitness_rows) > 0) paste(unique(fitness_rows$fitness_label), collapse = "; ") else NA_character_
  fitness_severity <- if (nrow(fitness_rows) > 0) paste(unique(fitness_rows$severity), collapse = "; ") else NA_character_
  fitness_evidence <- if (nrow(fitness_rows) > 0) paste(unique(fitness_rows$evidence), collapse = "; ") else NA_character_

  identity_rows <- if (!is.null(identity_phenotypes)) {
    identity_phenotypes %>% dplyr::filter(.data$cell_group == ct)
  } else {
    tibble::tibble()
  }
  identity_label <- if (nrow(identity_rows) > 0) paste(unique(identity_rows$identity_label), collapse = "; ") else NA_character_
  identity_evidence <- if (nrow(identity_rows) > 0) paste(unique(identity_rows$evidence), collapse = "; ") else NA_character_

  cell_type_degs <- degs %>%
    dplyr::filter(.data$cell_group == ct, .data$perturb_to_ctrl_p_value < sig_p_val_thresh)
  degs_of_interest <- if (!is.null(genes_of_interest)) {
    cell_type_degs %>% dplyr::filter(.data$gene_short_name %in% genes_of_interest)
  } else {
    cell_type_degs
  }
  goi <- format_goi(degs_of_interest)

  deg_genes <- unique(cell_type_degs$gene_short_name)
  fora_res_list <- list()
  for (ont in names(ontologies)) {
    gene_set <- ontologies[[ont]]
    if (!all(c("gs_name", "gene_short_name") %in% names(gene_set))) {
      next
    }
    pathway_gene_sets <- split(gene_set$gene_short_name, gene_set$gs_name)
    universe_genes <- unique(unlist(pathway_gene_sets, use.names = FALSE))
    if (length(universe_genes) == 0) {
      next
    }
    fora_res <- suppressWarnings(
      fgsea::fora(
        genes = deg_genes,
        pathways = pathway_gene_sets,
        universe = universe_genes
      )
    )
    fora_res_list[[ont]] <- fora_res %>%
      dplyr::arrange(.data$pval) %>%
      dplyr::mutate(ontology = ont)
  }
  fora_res <- dplyr::bind_rows(fora_res_list)

  expressed_targets <- character(0)
  if (!is.null(target_gene_expression)) {
    ancestor_expr <- target_gene_expression %>%
      dplyr::filter(.data$cell_group == ct)
    if (nrow(ancestor_expr) > 0) {
      expressed_targets <- unique(ancestor_expr$gene_short_name)
    }
  }
  expressed_targets_line <- if (length(expressed_targets) > 0) paste(expressed_targets, collapse = ", ") else "none"

  pathway_lines <- "Top perturbed pathways: none"
  if (!is.null(fora_res) && nrow(fora_res) > 0) {
    sig_pathways <- fora_res %>%
      dplyr::filter(.data$pval < sig_p_val_thresh) %>%
      dplyr::arrange(.data$pval) %>%
      dplyr::group_by(.data$ontology) %>%
      dplyr::slice_head(n = top_n_pathways) %>%
      dplyr::ungroup()
    if (nrow(sig_pathways) > 0) {
      pathway_lines <- paste0(
        "Top perturbed pathways:\n",
        paste(paste0(
          "- [", sig_pathways$ontology, "] ",
          vapply(sig_pathways$pathway, humanize_pathway_name, FUN.VALUE = character(1)),
          " (p=", signif(sig_pathways$pval, 2), ")"
        ), collapse = "\n")
      )
    }
  }

  list(
    abundance = abundance_summary,
    abundance_code = abundance_code,
    abundance_severity = abundance_severity,
    fitness_label = fitness_label,
    fitness_severity = fitness_severity,
    fitness_evidence = fitness_evidence,
    identity_label = identity_label,
    identity_evidence = identity_evidence,
    degs = cell_type_degs,
    goi_line = goi$goi_line,
    pathways = fora_res
  )
}

#' Summarize Lineage-Wide Impact Without LLM
#'
#' Traverses lineage-relevant cell types and compiles deterministic phenotype and
#' disrupted-pathway summaries for each included cell type.
#'
#' @param target_gene_expression Data frame of target-gene expression calls by
#'   cell group.
#' @param dact_results Differential abundance result table.
#' @param degs Differential expression table.
#' @param ref_expression Reference expression table (accepted for interface
#'   compatibility).
#' @param ontologies Named list of ontology gene-set tables.
#' @param combined_psg Cell-state graph object or graph used for lineage
#'   traversal.
#' @param sig_p_val_thresh Significance cutoff used in downstream summaries.
#' @param genes_of_interest Optional character vector used for DEG highlighting.
#' @param top_n_pathways Maximum pathways retained per ontology.
#' @param max_lineage_depth Optional maximum graph depth from roots for included
#'   cell types.
#' @param cell_types Optional explicit cell types to seed lineage traversal.
#' @param primary_impact_summary Reserved argument for interface compatibility.
#' @param excluded_cell_types Optional cell types to exclude from final output.
#' @param abundance_phenotypes Optional precomputed abundance phenotype table.
#' @param fitness_phenotypes Optional precomputed fitness phenotype table.
#' @param identity_phenotypes Optional precomputed identity phenotype table.
#' @param verbose Logical; if `TRUE`, prints debug messages.
#' @param ... Additional arguments accepted for interface compatibility.
#'
#' @return A tibble with one row per cell type and phenotype/pathway summary
#'   columns used by downstream reporting and plotting functions.
summarize_impact_in_lineage_context_no_llm <- function(
  target_gene_expression,
  dact_results,
  degs,
  ref_expression,
  ontologies,
  combined_psg,
  sig_p_val_thresh = 0.05,
  genes_of_interest = NULL,
  top_n_pathways = 5,
  max_lineage_depth = Inf,
  cell_types = NULL,
  primary_impact_summary = "",
  excluded_cell_types = NULL,
  abundance_phenotypes = NULL,
  fitness_phenotypes = NULL,
  identity_phenotypes = NULL,
  verbose = FALSE,
  ...
) {
  g <- if (class(combined_psg) == "cell_state_graph") combined_psg@graph else combined_psg

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

  if (!is.null(max_lineage_depth) && is.finite(max_lineage_depth) && max_lineage_depth >= 0) {
    roots <- unique(unlist(lapply(all_types, function(ct) get_roots(ct, combined_psg))))
    depth_ok <- vapply(all_types, function(ct) {
      if (length(roots) == 0) {
        return(TRUE)
      }
      min_depth <- suppressWarnings(min(unlist(lapply(roots, function(r) {
        d <- igraph::distances(coerce_state_graph(combined_psg), r, ct, mode = "out")
        as.numeric(d[1, 1])
      })), na.rm = TRUE))
      if (is.infinite(min_depth)) {
        return(FALSE)
      }
      min_depth <= max_lineage_depth
    }, FUN.VALUE = logical(1))
    all_types <- all_types[depth_ok]
  }

  if (!is.null(excluded_cell_types)) {
    all_types <- setdiff(all_types, excluded_cell_types)
  }

  remaining <- all_types
  results <- list()
  processed <- character(0)
  max_iter <- length(all_types) * 2
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

      if (verbose) {
        message(sprintf("[DEBUG] Processing cell type: %s (iteration %d)", ct, iter))
      }

      results[[ct]] <- summarize_cell_type_impact_no_ai_notes(
        ct,
        target_gene_expression,
        dact_results,
        degs,
        ontologies,
        sig_p_val_thresh,
        genes_of_interest,
        top_n_pathways = top_n_pathways,
        abundance_phenotypes = abundance_phenotypes,
        fitness_phenotypes = fitness_phenotypes,
        identity_phenotypes = identity_phenotypes
      )
      results[[ct]]$disrupted_pathways <- if (!is.null(excluded_cell_types) && ct %in% excluded_cell_types) {
        empty_disrupted_pathways_tbl()
      } else {
        build_deterministic_disrupted_pathways(
          pathways = results[[ct]]$pathways,
          sig_p_val_thresh = sig_p_val_thresh,
          top_n_pathways = top_n_pathways
        )
      }

      processed <- c(processed, ct)
      progress <- TRUE

      if (!verbose) {
        setTxtProgressBar(pb, length(processed) + 1)
      } else {
        message(sprintf("[DEBUG] %d cell types remaining to process.", length(all_types) - length(processed)))
      }
    }

    remaining <- setdiff(remaining, ready)
    if (!progress) {
      break
    }
  }

  if (!verbose) {
    close(pb)
  }

  out <- tibble::tibble(
    cell_type = names(results),
    data = unname(results)
  ) %>%
    tidyr::unnest_wider(.data$data)

  if (!"disrupted_pathways" %in% names(out)) {
    out$disrupted_pathways <- vector("list", nrow(out))
  }

  out <- out %>%
    dplyr::mutate(
      disrupted_pathways = purrr::map(.data$disrupted_pathways, function(x) {
        if (is.null(x) || !is.data.frame(x)) {
          empty_disrupted_pathways_tbl()
        } else {
          tibble::as_tibble(x)
        }
      })
    ) %>%
    dplyr::transmute(
      cell_type = .data$cell_type,
      abundance_code = .data$abundance_code,
      abundance_severity = .data$abundance_severity,
      identity_label = .data$identity_label,
      fitness_label = .data$fitness_label,
      fitness_evidence = .data$fitness_evidence,
      identity_evidence = .data$identity_evidence,
      disrupted_pathways = .data$disrupted_pathways
    )

  out
}
