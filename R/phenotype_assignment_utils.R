# Build a robust preranked vector from a DEG table
# deg_tbl columns (customize names via args):
#   gene, logFC, and (optionally) an empirical-p weight column
#
# When `weight_col` (default "empirical_p", written by the per-experiment
# empirical-FDR model) is present, genes are ranked by the MODEL-WEIGHTED
# statistic logFC * (1 - empirical_p): artifact-consistent genes (empirical_p
# near 1) are compressed toward the middle of the ranking, which halves the
# fitness/identity phenotype FDR without amplifying anything. Genes with no
# empirical_p (NA) are left un-compressed (weight 1). When the column is absent
# the ranking falls back to the raw shrunken logFC (backward compatible).
make_rank <- function(deg_tbl,
                      gene_col = "gene_short_name",
                      logFC_col = "perturb_to_ctrl_shrunken_lfc",
                      weight_col = "empirical_p") {
    stopifnot(all(c(gene_col, logFC_col) %in% names(deg_tbl)))
    use_weight <- !is.null(weight_col) && weight_col %in% names(deg_tbl)
    x <- deg_tbl %>%
        transmute(
            gene = .data[[gene_col]] %>% as.character(),
            logFC = as.numeric(.data[[logFC_col]]),
            p = if (use_weight) as.numeric(.data[[weight_col]]) else NA_real_
        ) %>%
        distinct(gene, .keep_all = TRUE) %>%
        filter(is.finite(logFC))
    if (use_weight) {
        # EXCLUDE gated (uncallable) genes -- empirical_p == 1 -- from the preranked list entirely,
        # rather than letting the (1 - empirical_p) weight collapse them to rank 0 in the MIDDLE. Under
        # heavy, asymmetric gating (e.g. the min-arm-cell gate demoting far more down- than up-calls),
        # a large mid-list block of zeros unbalances the ranking and manufactures spurious
        # one-directional GSEA enrichments (up-regulated fitness sets -- apoptosis/stress). Callable
        # genes keep the soft model weight (1 - empirical_p); NA weights are treated as un-demoted.
        x <- x %>% filter(!is.finite(.data$p) | .data$p < 1)
        x$w <- 1 - pmin(pmax(replace(x$p, !is.finite(x$p), 0), 0), 1)
    } else {
        x$w <- 1
    }
    # Rank by the model-weighted statistic (logFC * (1 - empirical_p)) over the CALLABLE genes; with
    # no weight column this reduces to the shrunken logFC over all genes.
    r <- x$logFC * x$w
    names(r) <- x$gene
    sort(r, decreasing = TRUE)
}

# Deterministic 31-bit hash of a string, used to derive a reproducible RNG seed
# from a call's own identity rather than from loop position. Double arithmetic
# (not integer) avoids 32-bit overflow; strength/collision-resistance doesn't
# matter here, only that identical keys always map to the same seed.
.fgsea_seed_from_key <- function(key) {
    h <- 5381
    for (b in utf8ToInt(enc2utf8(key))) h <- (h * 33 + b) %% 2147483647
    as.integer(h)
}

# Run fgsea on the preranked vector for a list of gene sets.
#
# Always dispatches to fgseaMultilevel (no `nperm` is ever forwarded to
# fgsea::fgsea() -- passing nperm there forces the deprecated fgseaSimple path,
# which floors padj resolution at ~1/(nperm+1) and emits a warning this
# function used to swallow via suppressWarnings). `nperm` is kept as the
# parameter name for backward compatibility with existing callers, but its
# value now maps onto fgseaMultilevel's `nPermSimple` (preliminary estimation
# permutations).
#
# Seeding: fgseaMultilevel is Monte Carlo, so identical inputs must get an
# identical seed to be reproducible. The seed is derived from `seed_key` --
# ideally the calling context's own identity (e.g. "<perturb_name>::<cell_group>::fitness"),
# passed in by the caller -- so results depend only on that call's own inputs,
# not on loop order, worker count, or serial vs. parallel execution. If no
# seed_key is supplied, one is derived from the gene ranking and gene sets
# themselves, which is still deterministic per-call but ties reproducibility
# to the data rather than to caller-supplied metadata.
run_fgsea_modules <- function(rank_vec,
                              gene_sets,
                              minSize = 10,
                              maxSize = 5000,
                              nperm = 1000,
                              seed_key = NULL) {
    # Keep only present genes per set
    gs <- lapply(gene_sets, function(v) intersect(unique(v), names(rank_vec)))
    keep <- lengths(gs) >= minSize
    if (!any(keep)) {
        return(tibble(path = character(), NES = numeric(), padj = numeric(), size = integer(), log2err = numeric(), leading_edge = character()))
    }
    gs <- gs[keep]
    if (is.null(seed_key)) {
        seed_key <- paste(c(names(rank_vec), names(gs)), collapse = "|")
    }
    seed <- .fgsea_seed_from_key(seed_key)
    res <- withr::with_seed(seed, {
        suppressWarnings({
            fgsea::fgsea(
                pathways = gs, stats = rank_vec,
                minSize = minSize, maxSize = maxSize, nPermSimple = nperm
            )
        })
    })
    res %>%
        as_tibble() %>%
        select(pathway, size, NES, padj, log2err, leadingEdge) %>%
        mutate(leading_edge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))) %>%
        select(pathway, size, NES, padj, log2err, leading_edge) %>%
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
                                   nperm = 1000,
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

log_ts <- function(...) {
    msg <- paste(..., collapse = " ")
    message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

get_phenotype_threads <- function(num_threads = NULL) {
    if (!is.null(num_threads) && is.numeric(num_threads) && num_threads > 0) {
        return(as.integer(num_threads))
    }
    1L
}

#' Assign abundance, fitness, and identity phenotype labels for every perturbation
#'
#' Iterates over the perturbations in `contrast_tbls`, loads each one's DEG
#' table, and classifies every cell type along three axes: abundance (`A*`
#' codes from the differential cell abundance test), fitness, and identity
#' (both from fgsea over the supplied gene sets).
#'
#' @param contrast_tbls A tibble with one row per perturbation, carrying the
#'   columns `perturb_name`, `perturb_group`, `run`, the list-columns
#'   `perturb_time_window` and `differential_expression_filename`, and the
#'   differential cell abundance list-column selected by `use_summarized_tbl`.
#' @param fitness_gene_sets Named list of gene sets used for the fitness fgsea.
#' @param identity_gene_sets Per-cell-type identity gene sets, as returned by
#'   [construct_identity_gene_sets()].
#' @param combined_psg Cell state graph supplying the lineage relationships
#'   (parents, descendants, alternative fates) used for identity labelling.
#' @param cell_type_denylist Character vector of cell types to drop before
#'   classification, or `NULL` to keep all of them.
#' @param num_threads Number of workers. `NULL` or a non-positive value runs
#'   serially.
#' @param use_summarized_tbl If `TRUE`, read the `summarized_differential_cell_abundance`
#'   column. If `FALSE`, derive the equivalent columns from
#'   `differential_cell_abundance` via [dacts_when_abundant()].
#' @param minSize,maxSize Gene set size bounds passed to fgsea.
#' @param nperm Permutation count passed to fgsea.
#' @param abundance_q_cut q-value cutoff for the abundance codes, passed to
#'   [assign_abundance_code()]. Defaults to `0.1`; pass `0.01` to reproduce the
#'   stricter `A3 Near-loss` gate used before this was configurable. Also passed
#'   to [assign_abundance_severity()], so that a cell type cannot come back as a
#'   non-call with a graded severity.
#'
#' @return A tibble with one row per perturbation and cell type, containing
#'   `cell_group`, the perturbation identifiers, the time window bounds,
#'   `abundance_code`, `abundance_severity`, and the `fitness_labels`,
#'   `fitness_fgsea`, `identity_labels` and `identity_fgsea` list-columns.
#'
#' @keywords internal
assign_phenotypes <- function(
  contrast_tbls,
  fitness_gene_sets,
  identity_gene_sets,
  combined_psg,
  cell_type_denylist = NULL,
  num_threads = NULL,
  use_summarized_tbl = TRUE,
  minSize = 10,
  maxSize = 5000,
  nperm = 1000,
  abundance_q_cut = 0.1
) {
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
    num_threads <- get_phenotype_threads(num_threads)
    log_ts(
        "assign_phenotypes start:",
        sprintf("perturbations=%d", nrow(contrast_tbls)),
        sprintf("total_cell_types=%d", length(all_cell_types)),
        sprintf("updates=%d", total_updates),
        sprintf("num_threads=%d", num_threads)
    )

    purrr::map_dfr(seq_len(nrow(contrast_tbls)), function(i) {
        perturb_record <- contrast_tbls[i, ]
        perturb_name <- perturb_record$perturb_name
        perturb_group <- perturb_record$perturb_group
        perturb_time_window <- perturb_record$perturb_time_window[[1]]
        run <- perturb_record$run

        deg_filename <- perturb_record$differential_expression_filename[[1]]

        # Load DEGs for this perturbation
        deg_tbl <- load_deg_file(deg_filename)
        deg_tbl <- deg_tbl %>% filter(present_above_thresh)

        # Use summarized differential cell abundance table
        if (use_summarized_tbl) {
            dact_tbl <- perturb_record$summarized_differential_cell_abundance[[1]]
            dact_tbl <- filter_denylisted_cell_types(dact_tbl, deg_tbl, cell_type_denylist)
            cell_type_universe <- unique(dact_tbl$cell_group)
        } else {
            # Optionally filter out denylisted cell types
            dact_tbl <- perturb_record$differential_cell_abundance[[1]]
            dact_tbl <- filter_denylisted_cell_types(dact_tbl, deg_tbl, cell_type_denylist)
            # The cell types this perturbation could have had something to say
            # about, captured BEFORE dacts_when_abundant() drops any.
            #
            # That filter requires present_above_thresh, which is a property of
            # the WILD-TYPE reference (fit_wt_model.R: percent_max_abund >= 0.1),
            # not of the perturbation: it marks states the experiment's timepoint
            # window could not assess. A cell type whose every row fails it
            # disappeared from the impact table entirely -- no row, no code, no
            # reason -- so "we could not assess this state here" was
            # indistinguishable from "this state was not in the experiment".
            # Measured on GENE6: 585 of 5,203 cell-type/perturbation pairs
            # (11.2%) vanished this way, every one of them for
            # present_above_thresh and none for a missing delta_log_abund.
            #
            # They are kept in the universe so they surface as "AN Not assessed"
            # rather than as silence.
            cell_type_universe <- unique(dact_tbl$cell_group)
            dact_tbl <- dacts_when_abundant(dact_tbl, percent_max_thresh = 0)
            # dact_tbl <- dact_tbl %>%
            #     mutate(timepoint_x = as.numeric(timepoint_x)) %>%
            #     group_by(cell_group) %>%
            #     slice_max(percent_max_abund, with_ties = F)
            dact_tbl <- dact_tbl %>% mutate(change_when_present = delta_log_abund, change_when_present_q_val = delta_q_value)
        }

        log_ts(
            sprintf(
                "perturbation %s (%d/%d): %d cell types",
                perturb_name, i, nrow(contrast_tbls), length(unique(dact_tbl$cell_group))
            )
        )

        assign_phenotypes_to_cell_types(
            dact_tbl, deg_tbl, fitness_gene_sets, identity_gene_sets,
            combined_psg,
            perturb_name, perturb_group, perturb_time_window, run,
            all_cell_types = cell_type_universe,
            pb = pb, # Pass the progress bar object
            log_fn = log_ts,
            num_threads = num_threads,
            minSize = minSize,
            maxSize = maxSize,
            nperm = nperm,
            abundance_q_cut = abundance_q_cut
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


#' Pick one abundance contrast per cell type, where the cell type was abundant
#'
#' Reduces a decorated Hooke abundance contrast table to one row per cell type:
#' the most significant timepoint among those where the cell type was actually
#' present in the wild-type reference.
#'
#' This is the canonical implementation. It was duplicated verbatim in
#' zscapetools for some time; that copy now calls this one. Do not reintroduce a
#' second copy -- one of the two was exported and the other was not, so bare
#' calls resolved by `library()` order.
#'
#' @section Rows are dropped, not flagged:
#' Cell types are **removed entirely** when every one of their rows has
#' `present_above_thresh` `FALSE` or `NA`, or a missing `delta_log_abund`. That
#' matters more than it looks: `present_above_thresh` is a property of the
#' WILD-TYPE reference, set in mcclintock's `fit_wt_model.R` as
#' `percent_max_abund >= 0.1`, so it marks states the experiment's timepoint
#' window could not assess -- not states the perturbation removed.
#'
#' On a 15-perturbation GENE6 run, 585 of 5,203 cell-type/perturbation pairs
#' (11.2%) were dropped this way, every one of them for `present_above_thresh`
#' and none for a missing `delta_log_abund`.
#'
#' Callers that build a per-cell-type table from this output therefore lose
#' those cell types silently. To report them instead, capture the cell-type
#' universe *before* calling this and pass it to
#' [assign_phenotypes_to_cell_types()] as `all_cell_types`, which surfaces them
#' as `"AN Not assessed"`.
#'
#' @param differential_cell_abundance A tibble of abundance contrasts carrying
#'   `cell_group`, `timepoint_x`, `delta_log_abund`, `delta_q_value`,
#'   `percent_max_abund` and `present_above_thresh`. The latter two are added
#'   downstream of Hooke by the wild-type reference join, so a raw
#'   `compare_abundances()` result will not have them.
#' @param percent_max_thresh Rows below this fraction of the cell type's peak
#'   abundance are neutralised rather than dropped: `delta_log_abund` is set to
#'   `0` and `delta_q_value` to `1`. Defaults to `0.10`. Pass `0` to neutralise
#'   nothing.
#' @param with_ties Passed to [dplyr::slice_min()]. `FALSE` (the default) keeps
#'   exactly one row per cell type even when several tie on `delta_q_value`.
#' @return A tibble with at most one row per `cell_group`. Cell types with no
#'   qualifying row are absent, not `NA` -- see the section above.
#' @examples
#' \dontrun{
#' # one row per cell type, neutralising nothing
#' dacts <- dacts_when_abundant(differential_cell_abundance, percent_max_thresh = 0)
#' }
#' @export
dacts_when_abundant <- function(differential_cell_abundance, percent_max_thresh = 0.10, with_ties = FALSE) {
    perturb_table_at_when_abundant <- differential_cell_abundance %>%
        mutate(timepoint_x = as.numeric(timepoint_x)) %>%
        filter(is.na(delta_log_abund) == FALSE &
            is.na(present_above_thresh) == FALSE &
            present_above_thresh) %>%
        # filter(present_above_thresh) %>%
        mutate(delta_log_abund = ifelse(percent_max_abund < percent_max_thresh, 0, delta_log_abund)) %>%
        mutate(delta_q_value = ifelse(delta_log_abund == 0, 1, delta_q_value)) %>%
        group_by(cell_group) %>%
        slice_min(delta_q_value, with_ties = with_ties) %>%
        ungroup()

    return(perturb_table_at_when_abundant)
}

#' Assign phenotype labels for every cell type within one perturbation
#'
#' Worker behind [assign_phenotypes()]; classifies each cell type in
#' `dact_tbl` along the abundance, fitness, and identity axes.
#'
#' @param dact_tbl Differential cell abundance table for this perturbation,
#'   carrying `cell_group`, `change_when_present` and
#'   `change_when_present_q_val`.
#' @param deg_tbl Differential expression table for this perturbation.
#' @param gene_sets Named list of gene sets used for the fitness fgsea.
#' @param identity_gene_sets Per-cell-type identity gene sets, as returned by
#'   [construct_identity_gene_sets()].
#' @param combined_psg Cell state graph supplying lineage relationships.
#' @param perturb_name,perturb_group,run Identifiers copied onto every output
#'   row.
#' @param perturb_time_window A list or data frame with `start_time` and
#'   `stop_time`, recorded as the output's time window bounds.
#' @param pb Optional progress bar object; ticked once per cell type.
#' @param log_fn Optional logging function called with progress messages.
#' @param num_threads Number of workers. `NULL` or a non-positive value runs
#'   serially.
#' @param minSize,maxSize Gene set size bounds passed to fgsea.
#' @param nperm Permutation count passed to fgsea.
#' @param abundance_q_cut q-value cutoff for the abundance codes, passed to both
#'   [assign_abundance_code()] and [assign_abundance_severity()] so the two stay
#'   consistent. Defaults to `0.1`.
#'
#' @return A tibble with one row per cell type; see [assign_phenotypes()] for
#'   the columns.
#'
#' @keywords internal
assign_phenotypes_to_cell_types <- function(
  dact_tbl, deg_tbl, gene_sets,
  identity_gene_sets,
  combined_psg,
  perturb_name, perturb_group, perturb_time_window, run,
  pb = NULL, # Accept progress bar object
  log_fn = NULL,
  num_threads = NULL,
  minSize = 10,
  maxSize = 5000,
  nperm = 1000,
  abundance_q_cut = 0.1,
  all_cell_types = NULL
) {
    # Defaults to the filtered table's own cell types, which is the pre-0.0.3
    # behaviour. Callers that know the wider universe pass it, so states the
    # experiment could not assess get a row saying so instead of vanishing.
    cell_types <- if (is.null(all_cell_types)) unique(dact_tbl$cell_group) else all_cell_types
    num_threads <- get_phenotype_threads(num_threads)
    worker_fn <- function(i) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BiocParallel::register(BiocParallel::SerialParam())
        }
        ct <- cell_types[[i]]
        if (!is.null(log_fn)) {
            log_fn(sprintf(
                "perturbation %s: cell type %s (%d/%d) start",
                perturb_name, ct, i, length(cell_types)
            ))
        }
        ct_start <- Sys.time()
        if (!is.null(pb)) pb$tick()
        dact_row <- dact_tbl %>% filter(cell_group == ct)
        # `$` on an absent column yields NULL, so this degrades to
        # "AU Undetermined" on tables that predate the SE/df columns.
        abundance_se <- dact_row$change_when_present_se
        abundance_df <- dact_row$change_when_present_tvalue_df
        # No surviving abundance row: the contrast never assessed this cell type
        # here. That is exactly "AN Not assessed", not an absent value -- and it
        # is why this branch is now reachable at all.
        abundance_code <- if (nrow(dact_row) > 0) assign_abundance_code(dact_row$change_when_present, dact_row$change_when_present_q_val, q_cut = abundance_q_cut, se = abundance_se, df = abundance_df) else "AN Not assessed"
        abundance_mdfc80_val <- if (nrow(dact_row) > 0) abundance_mdfc80(if (is.null(abundance_se)) NA_real_ else abundance_se, if (is.null(abundance_df)) NA_real_ else abundance_df, alpha = abundance_q_cut)[1] else NA_real_
        # q_cut must match the one the code cascade used, or a non-call can come
        # back graded.
        abundance_severity <- if (nrow(dact_row) > 0) assign_abundance_severity(dact_row$change_when_present, dact_row$change_when_present_q_val, q_cut = abundance_q_cut) else NA_character_
        identity_assignment <- assign_identity_maturation_labels(deg_tbl, ct, identity_gene_sets, combined_psg, minSize = minSize, nperm = nperm, maxSize = maxSize, perturb_name = perturb_name)
        fitness_assignment <- assign_fitness_labels(deg_tbl, ct, gene_sets, minSize = minSize, nperm = nperm, maxSize = maxSize, perturb_name = perturb_name)
        if (!is.null(log_fn)) {
            elapsed <- as.numeric(difftime(Sys.time(), ct_start, units = "secs"))
            log_fn(sprintf(
                "perturbation %s: cell type %s done (%.1fs)",
                perturb_name, ct, elapsed
            ))
        }
        tibble(
            cell_group = ct,
            perturb_group = perturb_group,
            perturb_name = perturb_name,
            run = run,
            time_window_start = perturb_time_window$start_time,
            time_window_end = perturb_time_window$stop_time,
            abundance_code = abundance_code,
            abundance_severity = abundance_severity,
            # Carry the evidence behind the code, so a negative result can be
            # evaluated without going back to the Hooke table.
            abundance_lfc = if (nrow(dact_row) > 0) dact_row$change_when_present[1] else NA_real_,
            abundance_se = if (is.null(abundance_se) || length(abundance_se) == 0) NA_real_ else as.numeric(abundance_se)[1],
            abundance_df = if (is.null(abundance_df) || length(abundance_df) == 0) NA_real_ else as.numeric(abundance_df)[1],
            abundance_mdfc80 = abundance_mdfc80_val,
            fitness_labels = list(fitness_assignment$labels),
            fitness_fgsea = list(fitness_assignment$fgsea_res),
            identity_labels = list(identity_assignment$labels),
            identity_fgsea = list(identity_assignment$fgsea_res)
        )
    }
    results <- if (num_threads > 1) {
        old_max <- getOption("future.globals.maxSize")
        options(future.globals.maxSize = 10 * 1024^3)
        on.exit(options(future.globals.maxSize = old_max), add = TRUE)
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        plan_strategy <- if (future::supportsMulticore()) {
            future::multicore
        } else {
            future::multisession
        }
        future::plan(plan_strategy, workers = num_threads)
        future.apply::future_lapply(seq_along(cell_types), worker_fn, future.seed = TRUE) %>%
            dplyr::bind_rows()
    } else {
        purrr::map_dfr(seq_along(cell_types), worker_fn)
    }
    results
}

#' Abundance codes that do not assert a phenotype
#'
#' `"A0 No change"` is a positive claim: the contrast was powered to detect a
#' change of the declared margin and saw none. `"AU Undetermined"` and
#' `"AN Not assessed"` are the absence of a claim. None of the three is a
#' phenotype, so anything asking "did this cell type change?" must treat all
#' three as no.
#'
#' Exported so that consumers outside platt -- zscape_portal in particular --
#' key off one definition rather than repeating the strings.
#'
#' @export
NON_CALLED_ABUNDANCE_CODES <- c("A0 No change", "AU Undetermined", "AN Not assessed")

# Minimum detectable fold change at `power`, on the scale `se` is measured on.
#
# Effect-independent by construction: a function of the standard error, alpha,
# the residual degrees of freedom and the requested power, never of the observed
# change. Returns NA where no detection limit is defined, so a caller can tell
# "we could not have seen it" apart from "there was nothing to see".
#
# Mirrors hooke::compare_abundances()'s `mdfc80`. Kept local because the
# summarised abundance table carries its own SE and df
# (`change_when_present_se`, `change_when_present_tvalue_df`), which are a
# weighted mean across the timepoints a cell type was present.
abundance_mdfc80 <- function(se, df, alpha = 0.1, power = 0.8) {
    se <- as.numeric(se)
    df <- as.numeric(df)
    if (length(se) == 0 || length(df) == 0) return(numeric(0))
    # Recycle to a common length BEFORE subsetting. Indexing a length-1 `df`
    # with a longer logical yields NA past the first element, which would make
    # a scalar df silently poison every row but the first.
    n <- max(length(se), length(df))
    se <- rep_len(se, n)
    df <- rep_len(df, n)
    usable <- is.finite(se) & se > 0 & is.finite(df) & df > 0
    out <- rep(NA_real_, length(se))
    if (any(usable)) {
        out[usable] <- exp(
            (qt(1 - alpha / 2, df[usable]) + qt(power, df[usable])) * se[usable]
        )
    }
    out
}

# Assign a four-state abundance code.
#
# Before this was a two-state cascade whose `TRUE ~` fallthrough was
# "A0 No change". That bucket absorbed, indistinguishably: a well-measured
# genuine null; a cell type whose standard error was so large it could not have
# caught a 68-fold change; and a degenerate fit reporting no estimate at all.
# Because the impact table carried only the binned code, "A0 No change" was
# unfalsifiable at the point of consumption.
#
# The four states, and what separates them:
#
#   A1/A2/A3         called   -- significant, as before
#   A0 No change     resolved -- not significant, AND mdfc80 <= margin, i.e. we
#                                were powered to see a change of that size and
#                                did not
#   AU Undetermined  not resolved -- not significant, and mdfc80 > margin (or
#                                unknown), i.e. we could not have seen such a
#                                change even if it were there
#   AN Not assessed  no usable contrast at all
#
# `margin_fold_change` is a FOLD CHANGE, and `mdfc80` is a fold change, so they
# are compared directly. (The design note wrote this as `mdfc80 < log(2)`, which
# is a units error: mdfc80 is always >= 1, so nothing would ever have resolved.)
#
# WHERE THE MARGIN COMES FROM. It is not a free parameter and it is not a round
# number. A "resolved null" claims we were powered to see a change we would have
# called a phenotype. The smallest change this function will call is `lfc_cut`
# on the log scale, so the margin that makes that claim true is exactly
# exp(lfc_cut) -- 1.649-fold at the default 0.5. It is derived from an existing
# declared threshold in the same way `.EFDR_MIN_EXPECTED = 3` is derived from
# alpha (Poisson P(0 | 3) ~ 5%), rather than being invented alongside it.
#
# A larger margin is not merely conservative, it is WRONG in the lenient
# direction. Set margin = 2 while calling phenotypes at 1.649, and every cell
# type with mdfc80 in (1.649, 2.0] is labelled a resolved null even though its
# detection limit exceeds the smallest change that would have counted -- the
# claim "no phenotype, and we would have caught one" is false for exactly those
# rows. That window is a ~39% band of standard-error space at every df the
# screens run at, so it is not a corner case.
#
# Because the margin is derived, tuning `lfc_cut` keeps the verdict coherent
# automatically. Override `margin_fold_change` only to answer a different
# question ("were we powered for a 2-fold change?"), not to set policy.
#
# When `se`/`df` are unavailable -- older tables that never carried them -- every
# non-significant row becomes "AU Undetermined". That is deliberate: without a
# standard error we cannot certify a null, and saying so is the point.
assign_abundance_code <- function(change_when_present, change_when_present_q_val,
                                  q_cut = 0.1,
                                  se = NULL, df = NULL,
                                  lfc_cut = 0.5,
                                  near_loss_cut = 2.0,
                                  margin_fold_change = exp(lfc_cut),
                                  power = 0.8) {
    n <- length(change_when_present)
    fill <- function(x) {
        if (is.null(x) || length(x) == 0) rep(NA_real_, n) else rep_len(as.numeric(x), n)
    }
    se <- fill(se)
    df <- fill(df)
    mdfc80 <- abundance_mdfc80(se, df, alpha = q_cut, power = power)
    resolved <- !is.na(mdfc80) & mdfc80 <= margin_fold_change

    # PRESENT BUT INVALID -> AN; ABSENT -> AU. Both inputs follow the same rule.
    #
    # A standard error that is present but zero or non-finite is a degenerate
    # fit -- the model reported no uncertainty at all. Residual df that is
    # present but <= 0 or non-finite is the same thing from the other side: an
    # overparameterised fit, where n - k - 1 went underwater. Neither was
    # really tested, so neither is "tested and inconclusive".
    #
    # An ABSENT se or df is a different situation: the row may be perfectly
    # fine, the table just never carried the column. We cannot certify a null
    # without it, so those fall through to AU below rather than claiming the
    # fit was broken.
    #
    # Without the df half of this, a row with a good se and an invalid df
    # produces mdfc80 = NA and lands in AU -- inconclusive -- when the fit
    # behind it never supported a test at all. No production row currently
    # does this (df runs 7-89 across v3.1.0, never NA, never <= 0), so this
    # closes a path rather than fixing an observed miscall.
    degenerate <- (!is.na(se) & (se <= 0 | !is.finite(se))) |
                  (!is.na(df) & (df <= 0 | !is.finite(df)))

    # ORDER IS LOAD-BEARING throughout: case_when takes the first match.
    case_when(
        is.na(change_when_present) | is.na(change_when_present_q_val) ~ "AN Not assessed",
        degenerate ~ "AN Not assessed",
        change_when_present >= lfc_cut & change_when_present_q_val < q_cut ~ "A1 Expansion",
        # A3 MUST be tested before A2, and its threshold is the more extreme of
        # the two on purpose. Read in isolation the next two lines look
        # backwards -- A2's -lfc_cut (-0.5) is a weaker cutoff than A3's
        # -near_loss_cut (-2.0), so every A3 row also satisfies A2. It is the
        # ORDER that separates them, not the thresholds: a -2.5 loss matches
        # A3 first and never reaches A2, while a -0.8 loss fails A3 and falls
        # to A2. Swap these two lines and A3 becomes unreachable.
        change_when_present <= -near_loss_cut & change_when_present_q_val < q_cut ~ "A3 Near-loss",
        change_when_present <= -lfc_cut & change_when_present_q_val < q_cut ~ "A2 Depletion",
        # Only non-calls reach here: anything significant with a real effect
        # size exited above.
        resolved ~ "A0 No change",
        TRUE ~ "AU Undetermined"
    )
}

# Graded severity for an abundance call.
#
# The "mild" tier is deliberately the SAME condition as being called A1/A2 at
# all, so that `severity == "none"` and a non-call always agree. That invariant
# used to hold by coincidence: this function hardcoded 0.5 and 0.1, which
# happened to equal assign_abundance_code()'s defaults. Once those became
# tunable arguments the coincidence broke -- setting lfc_cut = 0.8 produced
# cell types reported as "A0 No change" with "mild" severity, and q_cut was
# already user-facing, so that half was reachable in production.
#
# The shared cutoffs are therefore threaded, with the same defaults, so
# behaviour is unchanged unless a caller tunes them -- at which point both
# functions move together. `moderate_lfc_cut` and the stricter q levels are
# severity's own grading and have no counterpart in the code cascade.
assign_abundance_severity <- function(change_when_present, change_when_present_q_val,
                                      q_cut = 0.1,
                                      lfc_cut = 0.5,
                                      near_loss_cut = 2.0,
                                      moderate_lfc_cut = 1.0,
                                      severe_q_cut = 0.01,
                                      moderate_q_cut = 0.05) {
    # Every tier must IMPLY the calling condition, or a tier can fire on a row
    # the cascade rejected. Threading the shared cutoffs is not enough on its
    # own: severity's private q levels are stricter than q_cut by default, but a
    # caller tightening q_cut to 0.01 would leave `moderate` at 0.05 and grade
    # rows that were never called. So each tier is clamped to be no looser than
    # the calling condition on either axis.
    severe_q   <- min(severe_q_cut, q_cut)
    moderate_q <- min(moderate_q_cut, q_cut)
    moderate_l <- max(moderate_lfc_cut, lfc_cut)
    severe_l   <- max(near_loss_cut, moderate_l)

    case_when(
        abs(change_when_present) >= severe_l & change_when_present_q_val < severe_q ~ "severe",
        abs(change_when_present) >= moderate_l & change_when_present_q_val < moderate_q ~ "moderate",
        # same condition as A1/A2 in assign_abundance_code()
        abs(change_when_present) >= lfc_cut & change_when_present_q_val < q_cut ~ "mild",
        TRUE ~ "none"
    )
}

assign_fitness_labels <- function(deg_tbl, ct, gene_sets,
                                  minSize = 10,
                                  maxSize = 5000,
                                  nperm = 1000,
                                  perturb_name = NULL) {
    degs_this_cell <- deg_tbl %>% filter(cell_group == ct)
    if (nrow(degs_this_cell) > 0) {
        rank_vec <- make_rank(degs_this_cell, gene_col = "gene_short_name", logFC_col = "perturb_to_ctrl_shrunken_lfc")
        seed_key <- paste(perturb_name, ct, "fitness", sep = "::")
        fgsea_res <- run_fgsea_modules(rank_vec, gene_sets, minSize = minSize, maxSize = maxSize, nperm = nperm, seed_key = seed_key)
        list(
            labels = classify_fitness(fgsea_res),
            fgsea_res = fgsea_res
        )
    } else {
        list(
            labels = tibble(fitness_label = NA_character_, severity = NA_character_, evidence = NA_character_),
            fgsea_res = tibble()
        )
    }
}

assign_identity_maturation_labels <- function(deg_tbl, ct, identity_gene_sets, combined_psg,
                                              nes_mod = 1.5, nes_sev = 2.0, q_cut = 0.05,
                                              minSize = 10,
                                              maxSize = 5000,
                                              nperm = 1000,
                                              perturb_name = NULL) {
    degs_this_cell <- deg_tbl %>% filter(cell_group == ct)
    if (nrow(degs_this_cell) == 0) {
        return(list(
            labels = tibble(identity_label = "I0 Identity intact", evidence = NA_character_),
            fgsea_res = tibble()
        ))
    }

    # Get lineage info
    g <- coerce_state_graph(combined_psg)
    parents <- get_parents(g, ct)
    descendants <- get_descendants(ct, g)
    roots <- get_roots(ct, g)
    # lineage_tree <- if (length(roots) > 0) {
    #     unique(unlist(lapply(roots, function(r) get_descendants(r, g))))
    # } else {
    #     character()
    # }
    lineage_tree <- if (length(parents) > 0) {
        unique(unlist(lapply(parents, function(r) get_descendants(r, g))))
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
    seed_key <- paste(perturb_name, ct, "identity", sep = "::")
    fgsea_res <- run_fgsea_modules(rank_vec, gene_sets_list, minSize = minSize, maxSize = maxSize, nperm = nperm, seed_key = seed_key)

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

    list(
        labels = tibble(identity_label = label, evidence = evidence),
        fgsea_res = fgsea_res
    )
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

            fitness_fgsea_tsv <- file.path(out_dir, "fitness_fgsea.tsv")
            if ("fitness_fgsea" %in% colnames(.x)) {
                .x %>%
                    select(cell_group, fitness_fgsea) %>%
                    unnest(fitness_fgsea) %>%
                    readr::write_tsv(fitness_fgsea_tsv)
            }

            # Write identity/maturation info (unnest identity_labels if present)
            identity_tsv <- file.path(out_dir, "identity_phenotypes.tsv")
            if ("identity_labels" %in% colnames(.x)) {
                .x %>%
                    select(cell_group, identity_labels) %>%
                    unnest(identity_labels) %>%
                    readr::write_tsv(identity_tsv)
            }

            identity_fgsea_tsv <- file.path(out_dir, "identity_fgsea.tsv")
            if ("identity_fgsea" %in% colnames(.x)) {
                .x %>%
                    select(cell_group, identity_fgsea) %>%
                    unnest(identity_fgsea) %>%
                    readr::write_tsv(identity_fgsea_tsv)
            }
        })
}
# Example usage:
# write_phenotype_outputs(phenotype_tbl, base_dir = "all_rna_outputs/v3.0.0/mclintock")

################
# Gene set construction support


# Preranked (by specificity) gene vector for ONE cell type.
#
# Kept as a standalone helper so the cell-type subset is unit-testable without
# running fgsea, and so the subset is taken with base-R indexing rather than
# inside a dplyr data mask. That last point is the whole reason this function
# exists: callers may hand us a ref_expression carrying its own `cell_type`
# COLUMN (sulston's make_gene_sets created one, equal to cell_group, from
# 2026-01 onward). Inside filter(), such a column SHADOWS a function argument
# of the same name, silently turning `cell_group == cell_type` into a row-wise
# self-comparison that is TRUE for every row -- so no subset is taken, every
# gene keeps its maximum specificity across all cell types, and every cell type
# receives the same "identity" gene sets. Base-R indexing has no data mask and
# is immune to whatever columns the caller supplies.
.rank_genes_by_specificity <- function(ref_expression,
                                       cell_type,
                                       min_fraction_expressing = 0.01) {
    required <- c("cell_group", "gene_short_name", "fraction_expressing", "specificity")
    missing_cols <- setdiff(required, names(ref_expression))
    if (length(missing_cols) > 0) {
        stop(
            "ref_expression is missing required column(s): ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }
    keep <- ref_expression$cell_group == cell_type &
        ref_expression$fraction_expressing > min_fraction_expressing
    keep[is.na(keep)] <- FALSE
    rows <- ref_expression[keep, c("gene_short_name", "specificity"), drop = FALSE]
    rows <- rows[order(rows$specificity, decreasing = TRUE), , drop = FALSE]
    rows <- rows[!duplicated(rows$gene_short_name), , drop = FALSE]
    stats::setNames(rows$specificity, rows$gene_short_name)
}

# FIXME: Move this to zscapetools?
construct_identity_gene_sets <- function(ref_expression, gene_set_BP, gene_set_MF, gene_set_CC, cell_types, minSize = 15, maxSize = 500, padj_cutoff = 0.05, min_fraction_expressing = 0.01) {
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

    identity_gene_sets <- purrr::map_dfr(cell_types, function(this_cell_type) {
        pb$tick()
        # NOTE: the loop variable is deliberately NOT named `cell_type`. A
        # ref_expression column of that name would shadow it inside dplyr verbs.
        # See .rank_genes_by_specificity() for the full story.
        gene_ranks <- .rank_genes_by_specificity(
            ref_expression, this_cell_type,
            min_fraction_expressing = min_fraction_expressing
        )
        if (length(gene_ranks) < minSize) {
            return(tibble(cell_type = character(), gene_set_name = character(), gene_short_name = character()))
        }
        # Run fgsea for BP, MF, CC.
        #
        # fgseaMultilevel is Monte Carlo, so an unseeded run makes *which*
        # pathways clear padj_cutoff depend on RNG state -- i.e. on loop order,
        # worker count, and serial vs. parallel execution. Seed from this cell
        # type's own identity (same approach as run_fgsea_modules) so results
        # depend only on the inputs. Determinism additionally assumes the caller
        # has registered BiocParallel::SerialParam().
        go_res <- withr::with_seed(
            .fgsea_seed_from_key(paste(this_cell_type, "identity_go", sep = "::")),
            {
                suppressWarnings(list(
                    bp = fgsea::fgsea(
                        scoreType = "pos",
                        pathways = split(gene_set_BP$gene_short_name, gene_set_BP$gs_name),
                        stats = gene_ranks,
                        minSize = minSize,
                        maxSize = maxSize
                    ),
                    mf = fgsea::fgsea(
                        scoreType = "pos",
                        pathways = split(gene_set_MF$gene_short_name, gene_set_MF$gs_name),
                        stats = gene_ranks,
                        minSize = minSize,
                        maxSize = maxSize
                    ),
                    cc = fgsea::fgsea(
                        scoreType = "pos",
                        pathways = split(gene_set_CC$gene_short_name, gene_set_CC$gs_name),
                        stats = gene_ranks,
                        minSize = minSize,
                        maxSize = maxSize
                    )
                ))
            }
        )
        fgsea_go_bp <- go_res$bp
        fgsea_go_mf <- go_res$mf
        fgsea_go_cc <- go_res$cc
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
                    cell_type = this_cell_type,
                    gene_set_name = sig_pathways$display_name[i],
                    gene_short_name = sig_pathways$leadingEdge[[i]]
                )
            })
        }
    })
    identity_gene_sets %>% tidyr::unnest(gene_short_name)
}
