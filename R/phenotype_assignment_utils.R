# Build a robust preranked vector from a DEG table
# deg_tbl columns (customize names via args):
#   gene, logFC
make_rank <- function(deg_tbl,
                      gene_col = "gene_short_name",
                      logFC_col = "perturb_to_ctrl_shrunken_lfc") {
    stopifnot(all(c(gene_col, logFC_col) %in% names(deg_tbl)))
    x <- deg_tbl %>%
        transmute(
            gene = .data[[gene_col]] %>% as.character(),
            logFC = as.numeric(.data[[logFC_col]])
        ) %>%
        distinct(gene, .keep_all = TRUE) %>%
        filter(is.finite(logFC))
    # Use shrunken logFC directly for ranking
    r <- x$logFC
    names(r) <- x$gene
    sort(r, decreasing = TRUE)
}

# Run fgsea on the preranked vector for a list of gene sets
run_fgsea_modules <- function(rank_vec,
                              gene_sets,
                              minSize = 10,
                              maxSize = 5000,
                              nperm = 10000) {
    # Keep only present genes per set
    gs <- lapply(gene_sets, function(v) intersect(unique(v), names(rank_vec)))
    keep <- lengths(gs) >= minSize
    if (!any(keep)) {
        return(tibble(path = character(), NES = numeric(), padj = numeric(), size = integer()))
    }
    gs <- gs[keep]
    suppressWarnings({
        res <- fgsea::fgsea(
            pathways = gs, stats = rank_vec,
            minSize = minSize, maxSize = maxSize, nperm = nperm
        )
    })
    res %>%
        as_tibble() %>%
        select(pathway, size, NES, padj) %>%
        arrange(padj, desc(abs(NES))) %>%
        rename(path = pathway)
}

# Classify F1–F4 from module-level enrichment
classify_fitness <- function(fgsea_res,
                             nes_mod = 1.5, # moderate threshold
                             nes_sev = 2.0, # severe threshold
                             q_cut = 0.05, # FDR cutoff
                             stress_modules = c("p53", "UPR_ER", "Oxidative_Stress", "Hypoxia", "Interferon"),
                             prolif_modules = c("Proliferation_S", "Proliferation_G2M"),
                             apoptosis_modules = c("Apoptosis"),
                             senescence_modules = c("Senescence")) {
    pick <- function(paths) fgsea_res %>% filter(path %in% paths)

    # Proliferation (F1): direction & strength
    pro <- pick(prolif_modules)
    pro_up <- pro %>% filter(NES >= nes_mod, padj <= q_cut)
    pro_down <- pro %>% filter(NES <= -nes_mod, padj <= q_cut)

    f1 <- NULL
    if (nrow(pro_up) > 0 && nrow(pro_down) == 0) f1 <- "F1 Proliferation increase"
    if (nrow(pro_down) > 0 && nrow(pro_up) == 0) f1 <- "F1 Proliferation decrease"
    if (nrow(pro_up) > 0 && nrow(pro_down) > 0) f1 <- "F1 Proliferation reprogrammed (mixed S/G2M)"

    # Apoptosis (F2)
    apo <- pick(apoptosis_modules) %>% filter(padj <= q_cut)
    f2 <- if (nrow(apo) && any(apo$NES >= nes_mod)) "F2 Apoptosis activation" else NULL

    # Stress/toxicity (F3)
    str <- pick(stress_modules) %>% filter(padj <= q_cut)
    f3 <- if (nrow(str) && any(str$NES >= nes_mod)) {
        pos <- str %>%
            filter(NES >= nes_mod) %>%
            arrange(desc(NES)) %>%
            pull(path)
        paste0("F3 Stress/toxicity (", paste(pos, collapse = ", "), ")")
    } else {
        NULL
    }

    # Senescence/quiescence (F4): senescence up + proliferation down
    sen <- pick(senescence_modules) %>% filter(padj <= q_cut, NES >= nes_mod)
    f4 <- if (nrow(sen) && nrow(pro_down) > 0) "F4 Senescence/quiescence shift" else NULL

    # Severity heuristic per label = max NES of its evidence
    sev_of <- function(paths, sign = "pos") {
        tab <- fgsea_res %>% filter(path %in% paths, padj <= q_cut)
        if (!nrow(tab)) {
            return(0)
        }
        nes <- if (sign == "neg") -tab$NES else tab$NES
        max(nes, na.rm = TRUE)
    }
    sev_map <- list(
        "F1 Proliferation increase" = sev_of(prolif_modules, "pos"),
        "F1 Proliferation decrease" = sev_of(prolif_modules, "neg"),
        "F1 Proliferation reprogrammed (mixed S/G2M)" = max(sev_of("Proliferation_S", "pos"), sev_of("Proliferation_G2M", "pos")),
        "F2 Apoptosis activation" = sev_of(apoptosis_modules, "pos"),
        "F3 Stress/toxicity" = sev_of(stress_modules, "pos"),
        "F4 Senescence/quiescence shift" = min(sev_of(senescence_modules, "pos"), sev_of(prolif_modules, "neg"))
    )

    labels <- c(f1, f2, f3, f4) %>%
        unique() %>%
        .[!is.null(.)]
    if (length(labels) == 0) labels <- "F0 No significant phenotype"

    out <- tibble(
        fitness_label = labels,
        severity = sapply(labels, function(L) {
            nes <- sev_map[[ifelse(L == "F3 Stress/toxicity", "F3 Stress/toxicity", L)]]
            if (is.null(nes) || is.na(nes)) {
                return("mild")
            }
            if (nes >= nes_sev) "severe" else if (nes >= nes_mod) "moderate" else "mild"
        }),
        evidence = sapply(labels, function(L) {
            if (str_starts(L, "F1")) {
                pro %>%
                    filter(padj <= q_cut) %>%
                    transmute(e = sprintf("%s NES=%.2f q=%.3g", path, NES, padj)) %>%
                    pull(e) %>%
                    paste(collapse = "; ")
            } else if (str_starts(L, "F2")) {
                apo %>%
                    transmute(e = sprintf("%s NES=%.2f q=%.3g", path, NES, padj)) %>%
                    pull(e) %>%
                    paste(collapse = "; ")
            } else if (str_starts(L, "F3")) {
                str %>%
                    transmute(e = sprintf("%s NES=%.2f q=%.3g", path, NES, padj)) %>%
                    pull(e) %>%
                    paste(collapse = "; ")
            } else if (str_starts(L, "F4")) {
                paste(
                    sen %>% transmute(e = sprintf("Senescence NES=%.2f q=%.3g", NES, padj)) %>% pull(e),
                    pro_down %>% transmute(e = sprintf("%s NES=%.2f q=%.3g", path, NES, padj)) %>% pull(e),
                    sep = " | "
                )
            } else {
                ""
            }
        })
    )
    out
}

# One convenience wrapper per cell type
score_fitness_from_deg <- function(deg_tbl,
                                   gene_sets,
                                   gene_col = "gene",
                                   logFC_col = "logFC",
                                   padj_col = "padj",
                                   minSize = 10,
                                   nperm = 10000,
                                   nes_mod = 1.5,
                                   nes_sev = 2.0,
                                   q_cut = 0.05) {
    r <- make_rank(deg_tbl, gene_col, logFC_col, padj_col)
    fg <- run_fgsea_modules(r, gene_sets, minSize = minSize, nperm = nperm)
    classes <- classify_fitness(fg,
        nes_mod = nes_mod, nes_sev = nes_sev, q_cut = q_cut
    )
    list(
        fgsea_table = fg,
        fitness_labels = classes
    )
}

# Build custom gene sets for fgsea modules from search terms
build_custom_modules <- function(gene_set_BP, module_definitions) {
    # Preprocess gene set names: convert underscores to spaces and to lower case
    gene_set_BP <- gene_set_BP %>%
        mutate(gs_name_clean = stringr::str_to_lower(stringr::str_replace_all(gs_name, "_", " ")))

    custom_sets <- purrr::map(module_definitions, function(terms) {
        # Convert terms to lower case for matching
        terms_clean <- stringr::str_to_lower(terms)
        matched_sets <- gene_set_BP %>%
            filter(
                purrr::map_lgl(gs_name_clean, function(gs) {
                    any(stringr::str_detect(gs, terms_clean))
                })
            )
        unique(matched_sets$gene_short_name)
    })
    names(custom_sets) <- names(module_definitions)
    custom_sets
}

load_deg_file <- function(deg_out_filename, deg_q_val_thresh = 1.0, cell_type_denylist = c()) {
    tryCatch(
        {
            deg_tbl <- data.table::fread(deg_out_filename)
            deg_tbl$cell_group <- stringr::str_trim(deg_tbl$cell_group)
            deg_tbl <- deg_tbl %>% filter(cell_group %in% cell_type_denylist == FALSE)
            # TODO: consider input validation with "problems", below.
            ## print (problems(deg_tbl))

            # FIXME: remove this filter, shouldn't be needed:
            # deg_tbl <- deg_tbl %>% mutate(perturb_to_ctrl_p_value = ifelse(abs(perturb_to_ctrl_shrunken_lfc) > 15, 1.0, perturb_to_ctrl_p_value))
            # deg_tbl <- deg_tbl %>% mutate(perturb_to_ctrl_shrunken_lfc = ifelse(abs(perturb_to_ctrl_shrunken_lfc) > 15, 0, perturb_to_ctrl_shrunken_lfc))

            deg_tbl <- deg_tbl %>% mutate(perturb_to_ctrl_q_value = p.adjust(perturb_to_ctrl_p_value))
            deg_tbl <- deg_tbl %>% filter(perturb_to_ctrl_q_value <= deg_q_val_thresh)


            deg_tbl
        },
        error = function(e) {
            print(e)
            NULL
        }
    )
}

assign_phenotypes <- function(contrast_tbls, fitness_gene_sets, identity_gene_sets, combined_psg, cell_type_denylist = NULL) {
    # Get all cell types across all perturbations
    all_cell_types <- unique(unlist(
        lapply(seq_len(nrow(contrast_tbls)), function(i) {
            perturb_record <- contrast_tbls[i, ]
            dact_tbl <- perturb_record$differential_cell_abundance[[1]]
            dact_tbl$cell_group
        })
    ))
    total_updates <- length(all_cell_types) * nrow(contrast_tbls)
    pb <- progress::progress_bar$new(
        format = "  Processing [:bar] :percent in :elapsed",
        total = total_updates,
        clear = FALSE,
        width = 60
    )

    purrr::map_dfr(seq_len(nrow(contrast_tbls)), function(i) {
        perturb_record <- contrast_tbls[i, ]
        perturb_name <- perturb_record$perturb_name
        perturb_group <- perturb_record$perturb_group
        perturb_time_window <- perturb_record$perturb_time_window[[1]]
        run <- perturb_record$run
        dact_tbl <- perturb_record$differential_cell_abundance[[1]]
        deg_filename <- perturb_record$differential_expression_filename[[1]]

        # Load DEGs for this perturbation
        deg_tbl <- load_deg_file(deg_filename)

        # Optionally filter out denylisted cell types
        dact_tbl <- filter_denylisted_cell_types(dact_tbl, deg_tbl, cell_type_denylist)

        # Use summarized differential cell abundance table
        dact_tbl <- perturb_record$summarized_differential_cell_abundance[[1]]

        assign_phenotypes_to_cell_types(
            dact_tbl, deg_tbl, fitness_gene_sets, identity_gene_sets,
            combined_psg,
            perturb_name, perturb_group, perturb_time_window, run,
            pb = pb # Pass the progress bar object
        )
    })
}

filter_denylisted_cell_types <- function(dact_tbl, deg_tbl, cell_type_denylist) {
    if (!is.null(cell_type_denylist)) {
        dact_tbl <- dact_tbl %>% filter(!cell_group %in% cell_type_denylist)
        deg_tbl <- deg_tbl %>% filter(!cell_group %in% cell_type_denylist)
    }
    dact_tbl
}


assign_phenotypes_to_cell_types <- function(
  dact_tbl, deg_tbl, gene_sets,
  identity_gene_sets,
  combined_psg,
  perturb_name, perturb_group, perturb_time_window, run,
  pb = NULL # Accept progress bar object
) {
    cell_types <- unique(dact_tbl$cell_group)
    results <- purrr::map_dfr(cell_types, function(ct) {
        if (!is.null(pb)) pb$tick()
        dact_row <- dact_tbl %>% filter(cell_group == ct)
        abundance_code <- if (nrow(dact_row) > 0) assign_abundance_code(dact_row$change_when_present, dact_row$change_when_present_q_val) else NA_character_
        abundance_severity <- if (nrow(dact_row) > 0) assign_abundance_severity(dact_row$change_when_present, dact_row$change_when_present_q_val) else NA_character_
        identity_labels <- assign_identity_maturation_labels(deg_tbl, ct, identity_gene_sets, combined_psg)
        fitness_labels <- assign_fitness_labels(deg_tbl, ct, gene_sets)
        tibble(
            cell_group = ct,
            perturb_group = perturb_group,
            perturb_name = perturb_name,
            run = run,
            time_window_start = perturb_time_window$start_time,
            time_window_end = perturb_time_window$stop_time,
            abundance_code = abundance_code,
            abundance_severity = abundance_severity,
            fitness_labels = list(fitness_labels),
            identity_labels = list(identity_labels)
        )
    })
    results
}

assign_abundance_code <- function(change_when_present, change_when_present_q_val) {
    case_when(
        is.na(change_when_present) | is.na(change_when_present_q_val) ~ "A0 No change",
        change_when_present >= 0.5 & change_when_present_q_val < 0.1 ~ "A1 Expansion",
        change_when_present <= -0.5 & change_when_present_q_val < 0.1 ~ "A2 Depletion",
        change_when_present <= -2.0 & change_when_present_q_val < 0.01 ~ "A3 Near-loss",
        TRUE ~ "A0 No change"
    )
}

assign_abundance_severity <- function(change_when_present, change_when_present_q_val) {
    case_when(
        abs(change_when_present) >= 2.0 & change_when_present_q_val < 0.01 ~ "severe",
        abs(change_when_present) >= 1.0 & change_when_present_q_val < 0.05 ~ "moderate",
        abs(change_when_present) >= 0.5 & change_when_present_q_val < 0.1 ~ "mild",
        TRUE ~ "none"
    )
}

assign_fitness_labels <- function(deg_tbl, ct, gene_sets) {
    degs_this_cell <- deg_tbl %>% filter(cell_group == ct)
    if (nrow(degs_this_cell) > 0) {
        rank_vec <- make_rank(degs_this_cell, gene_col = "gene_short_name", logFC_col = "perturb_to_ctrl_shrunken_lfc")
        fgsea_res <- run_fgsea_modules(rank_vec, gene_sets)
        classify_fitness(fgsea_res)
    } else {
        tibble(label = NA_character_, severity = NA_character_, evidence = NA_character_)
    }
}

assign_identity_maturation_labels <- function(deg_tbl, ct, identity_gene_sets, combined_psg,
                                              nes_mod = 1.5, nes_sev = 2.0, q_cut = 0.05) {
    degs_this_cell <- deg_tbl %>% filter(cell_group == ct)
    if (nrow(degs_this_cell) == 0) {
        return(tibble(identity_label = "I0 Identity intact", evidence = NA_character_))
    }

    # Get lineage info
    parents <- get_dir_parents(ct, combined_psg@graph)
    descendants <- get_descendants(ct, combined_psg@graph)
    roots <- get_roots(ct, combined_psg@graph)
    lineage_tree <- if (length(roots) > 0) {
        unique(unlist(lapply(roots, function(r) get_descendants(r, combined_psg@graph))))
    } else {
        character()
    }
    alt_fates <- setdiff(lineage_tree, c(ct, parents, descendants))

    # Filter identity_gene_sets to include only relevant cell types
    relevant_cell_types <- unique(c(ct, parents, descendants, alt_fates))
    filtered_identity_gene_sets <- identity_gene_sets %>%
        filter(cell_type %in% relevant_cell_types)

    # Build individual gene sets for relevant cell types
    gene_sets_list <- filtered_identity_gene_sets %>%
        group_by(cell_type, gene_set_name) %>%
        summarise(genes = list(unique(gene_short_name)), .groups = "drop") %>%
        mutate(set_id = paste(cell_type, gene_set_name, sep = "::")) %>%
        {
            setNames(.$genes, .$set_id)
        }

    # Prepare for fgsea
    rank_vec <- make_rank(degs_this_cell, gene_col = "gene_short_name", logFC_col = "perturb_to_ctrl_shrunken_lfc")
    fgsea_res <- run_fgsea_modules(rank_vec, gene_sets_list)

    # Annotate each result with its cell type and gene set
    fgsea_res <- fgsea_res %>%
        mutate(
            cell_type = sub("::.*", "", path),
            gene_set_name = sub(".*::", "", path)
        )

    evidence_list <- list()
    scores <- c(I1 = NA, I2 = NA, I4 = NA)
    evidence <- NA_character_
    label <- "I0 Identity intact"

    # I1: Maturation delay (upstream/parent sets upregulated)
    if (length(parents) > 0) {
        maturation_hits <- fgsea_res %>% filter(cell_type %in% parents, NES >= nes_mod, padj <= q_cut)
        if (nrow(maturation_hits) > 0) {
            scores["I1"] <- max(maturation_hits$NES)
            evidence_list$I1 <- paste0(
                "Upregulated parent/progenitor sets: ",
                paste(maturation_hits$gene_set_name, collapse = ", ")
            )
        }
    }

    # I2: Precocious maturation (downstream sets upregulated)
    if (length(descendants) > 0) {
        downstream_hits <- fgsea_res %>% filter(cell_type %in% descendants, NES >= nes_mod, padj <= q_cut)
        if (nrow(downstream_hits) > 0) {
            scores["I2"] <- max(downstream_hits$NES)
            evidence_list$I2 <- paste0(
                "Upregulated downstream/terminal sets: ",
                paste(downstream_hits$gene_set_name, collapse = ", ")
            )
        }
    }

    # I4: Fate switch/misspecification (alternative fate sets upregulated)
    if (length(alt_fates) > 0) {
        alt_hits <- fgsea_res %>% filter(cell_type %in% alt_fates, NES >= nes_mod, padj <= q_cut)
        if (nrow(alt_hits) > 0) {
            scores["I4"] <- max(alt_hits$NES)
            evidence_list$I4 <- paste0(
                "Upregulated alternative fate sets: ",
                paste(alt_hits$gene_set_name, collapse = ", ")
            )
        }
    }

    # Assign the strongest phenotype (highest NES among I1, I2, I4)
    valid_scores <- scores[!is.na(scores)]
    strongest <- if (length(valid_scores) > 0) names(valid_scores)[which.max(valid_scores)] else NA
    if (!is.na(strongest) && valid_scores[strongest] >= nes_mod) {
        label <- switch(strongest,
            I1 = "I1 Maturation delay",
            I2 = "I2 Precocious maturation",
            I4 = "I4 Fate switch / misspecification"
        )
        evidence <- evidence_list[[strongest]]
    }

    # I3: Program failure within identity (identity sets downregulated, no strong I1/I2/I4)
    identity_down <- fgsea_res %>% filter(cell_type == ct, NES <= -nes_mod, padj <= q_cut)
    if (nrow(identity_down) > 0 && label == "I0 Identity intact") {
        label <- "I3 Program failure within identity"
        evidence <- paste(
            "Identity sets downregulated (possible effector failure):",
            paste(identity_down$gene_set_name, collapse = ", ")
        )
    }

    # If no phenotype, check for identity set upregulation (I0)
    identity_up <- fgsea_res %>% filter(cell_type == ct, NES >= nes_mod, padj <= q_cut)
    if (nrow(identity_up) > 0 && label == "I0 Identity intact") {
        label <- "I0 Identity intact"
        evidence <- paste("Identity sets upregulated:", paste(identity_up$gene_set_name, collapse = ", "))
    }

    tibble(identity_label = label, evidence = evidence)
}

write_phenotype_outputs <- function(phenotype_tbl, base_dir) {
    phenotype_tbl %>%
        group_by(perturb_group, run, perturb_name) %>%
        group_walk(~ {
            # Output directory: <base_dir>/<perturb_group>/run_<run>/phenotypes/perturb_<perturb_name>/
            # out_dir <- file.path(base_dir, .y$perturb_group, paste0("run_", as.character(.y$run)), "phenotypes", paste0("perturb_", .y$perturb_name))

            out_dir <- file.path(base_dir, "phenotypes")
            if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

            # Write abundance info
            abundance_tsv <- file.path(out_dir, "abundance_phenotypes.tsv")
            .x %>%
                select(cell_group, abundance_code, abundance_severity) %>%
                readr::write_tsv(abundance_tsv)

            # Write fitness info (unnest fitness_labels if present)
            fitness_tsv <- file.path(out_dir, "fitness_phenotypes.tsv")
            if ("fitness_labels" %in% colnames(.x)) {
                .x %>%
                    select(cell_group, fitness_labels) %>%
                    unnest(fitness_labels) %>%
                    readr::write_tsv(fitness_tsv)
            }

            # Write identity/maturation info (unnest identity_labels if present)
            identity_tsv <- file.path(out_dir, "identity_phenotypes.tsv")
            if ("identity_labels" %in% colnames(.x)) {
                .x %>%
                    select(cell_group, identity_labels) %>%
                    unnest(identity_labels) %>%
                    readr::write_tsv(identity_tsv)
            }
        })
}
# Example usage:
# write_phenotype_outputs(phenotype_tbl, base_dir = "all_rna_outputs/v3.0.0/mclintock")

################
# Gene set construction support


# FIXME: Move this to zscapetools?
construct_identity_gene_sets <- function(ref_expression, gene_set_BP, gene_set_MF, gene_set_CC, cell_types, minSize = 15, maxSize = 500, padj_cutoff = 0.05) {
    # library(dplyr)
    # library(fgsea)
    # library(stringr)
    # library(progress)

    pb <- progress::progress_bar$new(
        format = "  Processing [:bar] :percent in :elapsed",
        total = length(cell_types),
        clear = FALSE,
        width = 60
    )

    identity_gene_sets <- purrr::map_dfr(cell_types, function(cell_type) {
        pb$tick()
        specificity_ranked_genes <- ref_expression %>%
            filter(cell_group == cell_type, fraction_expressing > 0.01)
        gene_ranks <- specificity_ranked_genes %>%
            arrange(desc(specificity)) %>%
            distinct(gene_short_name, .keep_all = TRUE) %>%
            select(gene_short_name, specificity) %>%
            deframe()
        # Run fgsea for BP, MF, CC
        fgsea_go_bp <- fgsea::fgsea(
            scoreType = "pos",
            pathways = split(gene_set_BP$gene_short_name, gene_set_BP$gs_name),
            stats = gene_ranks,
            minSize = minSize,
            maxSize = maxSize
        )
        fgsea_go_mf <- fgsea::fgsea(
            scoreType = "pos",
            pathways = split(gene_set_MF$gene_short_name, gene_set_MF$gs_name),
            stats = gene_ranks,
            minSize = minSize,
            maxSize = maxSize
        )
        fgsea_go_cc <- fgsea::fgsea(
            scoreType = "pos",
            pathways = split(gene_set_CC$gene_short_name, gene_set_CC$gs_name),
            stats = gene_ranks,
            minSize = minSize,
            maxSize = maxSize
        )
        # Collect significant pathways
        sig_bp <- fgsea_go_bp %>% filter(padj < padj_cutoff)
        if (nrow(sig_bp) > 0) sig_bp$ontology <- "BP"

        sig_mf <- fgsea_go_mf %>% filter(padj < padj_cutoff)
        if (nrow(sig_mf) > 0) sig_mf$ontology <- "MF"

        sig_cc <- fgsea_go_cc %>% filter(padj < padj_cutoff)
        if (nrow(sig_cc) > 0) sig_cc$ontology <- "CC"

        sig_pathways <- bind_rows(sig_bp, sig_mf, sig_cc)
        # Parse pathway names for display
        sig_pathways <- sig_pathways %>%
            mutate(
                display_name = pathway %>%
                    stringr::str_replace("^(GOBP_|GOMF_|GOCC_)", "") %>%
                    stringr::str_replace_all("_", " ") %>%
                    stringr::str_to_sentence()
            )
        # Return gene sets in tidy format
        if (nrow(sig_pathways) == 0) {
            return(tibble(cell_type = character(), gene_set_name = character(), gene_short_name = character()))
        } else {
            purrr::map_dfr(1:nrow(sig_pathways), function(i) {
                tibble(
                    cell_type = cell_type,
                    gene_set_name = sig_pathways$display_name[i],
                    gene_short_name = sig_pathways$leadingEdge[[i]]
                )
            })
        }
    })
    identity_gene_sets %>% tidyr::unnest(gene_short_name)
}
