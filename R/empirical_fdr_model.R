# Empirical-FDR model for the control-sampling artifact.
#
# In control-heavy perturbation designs the deep control arm detects lowly
# expressed genes the shallow perturbation arm misses, producing an excess of
# low-expression, on-in-control (down) DEG calls that are sampling artifacts,
# not biology. This module estimates a per-experiment empirical null for the
# standardized DEG statistic z = shrunken_lfc / shrunken_lfc_se by control-split
# permutation, learns the null tail as a function of (expression, direction,
# sampling log-ratio), and assigns a per-gene empirical FDR to real DEG calls.
#
# Built ONCE PER EXPERIMENT on a cheap stratified subsample, across a small grid
# of control:perturbation log-ratios bracketing the experiment's perturbations.
#
# @importFrom dplyr filter mutate group_by ungroup summarise transmute select
#   arrange if_else n cur_group across all_of bind_rows
# @importFrom purrr map map_dfr set_names
# @importFrom tidyr replace_na
# @importFrom rlang .data
# @importFrom magrittr %>%
NULL

.EFDR_BREAKS <- c(-Inf, seq(-8, 0.5, 0.5), Inf)

#' Build a control-split empirical null for the DEG artifact.
#'
#' Relabels control samples as a pseudo-perturbation at several split sizes, runs
#' the standard within-state DEG contrast on a stratified subsample of cell types
#' (and, optionally, genes), and returns the pooled null gene-tests tagged with
#' their sampling log-ratio. Training data for [train_efdr_model()].
#'
#' @param cds A cell_data_set containing the experiment's control cells.
#' @param sample_group,cell_group colData columns for replicate and cell type.
#' @param control_ids,perturbation_col control values and perturbation column.
#' @param n_perturb_grid Sizes of the pseudo-perturbation (control samples
#'   relabelled). If NULL, chosen from the control count to give control-heavy
#'   log-ratios bracketing real perturbations.
#' @param n_cell_types,min_cells_per_type Stratified cell-type panel controls.
#' @param gene_subsample Optional cap on genes fit per cell type (expression-spanning).
#' @param nperm,cores,max_simultaneous_genes,seed,write_dir,verbose Pass-through / control.
#' @return A tibble with columns cell_group, log_mean_expression, z, log_ratio.
#' @export
build_empirical_null <- function(cds,
                                 sample_group = "embryo_ID",
                                 cell_group = "cell_type",
                                 control_ids = c("ctrl-inj"),
                                 perturbation_col = "perturbation",
                                 n_perturb_grid = NULL,
                                 n_cell_types = 80L,
                                 min_cells_per_type = 100L,
                                 gene_subsample = NULL,
                                 nperm = 10000L,
                                 cores = 1L,
                                 seed = 42L,
                                 max_simultaneous_genes = 2000L,
                                 write_dir = tempfile("efdr_null_"),
                                 verbose = TRUE) {
  stopifnot(methods::is(cds, "cell_data_set"))
  cd <- SummarizedExperiment::colData(cds)
  cds <- cds[, as.character(cd[[perturbation_col]]) %in% control_ids]
  cd <- SummarizedExperiment::colData(cds)
  embs <- unique(as.character(cd[[sample_group]]))
  n_ctrl <- length(embs)

  if (is.null(n_perturb_grid)) {
    n_perturb_grid <- unique(pmax(3L, round(c(0.15, 0.20, 0.25, 0.30, 0.40) * n_ctrl)))
  }
  n_perturb_grid <- n_perturb_grid[n_perturb_grid > 0 & n_perturb_grid < n_ctrl]

  # abundance-stratified cell-type panel
  ok <- tibble::tibble(cell_group = as.character(cd[[cell_group]])) %>%
    dplyr::count(.data$cell_group, name = "n") %>%
    dplyr::filter(.data$n >= min_cells_per_type) %>%
    dplyr::arrange(.data$n) %>%
    dplyr::pull(.data$cell_group)
  if (length(ok) > n_cell_types) ok <- ok[unique(round(seq(1, length(ok), length.out = n_cell_types)))]

  keep_genes <- rownames(cds)
  if (!is.null(gene_subsample) && gene_subsample < nrow(cds)) {
    set.seed(seed)
    tot <- Matrix::rowSums(monocle3::exprs(cds)[, sample(ncol(cds), min(2000L, ncol(cds)))])
    keep_genes <- rownames(cds)[order(tot)][unique(round(seq(1, nrow(cds), length.out = gene_subsample)))]
  }

  dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
  purrr::map_dfr(n_perturb_grid, function(N) {
    set.seed(seed + N)
    null_pert <- sample(embs, N)
    sub <- cds[keep_genes, ]
    SummarizedExperiment::colData(sub)[[perturbation_col]] <-
      dplyr::if_else(as.character(SummarizedExperiment::colData(sub)[[sample_group]]) %in% null_pert,
                     "NULLPERT", "ctrl-inj")
    sub <- sub[, as.character(SummarizedExperiment::colData(sub)[[cell_group]]) %in% ok]
    ccs <- hooke::new_cell_count_set(sub, sample_group = sample_group, cell_group = cell_group)
    wd <- file.path(write_dir, paste0("N", N)); dir.create(wd, showWarnings = FALSE)
    if (verbose) message(sprintf("[efdr null] N=%d  log2ratio=%.2f  cts=%d genes=%d",
                                 N, log2(N / (n_ctrl - N)), length(ok), length(keep_genes)))
    compare_genes_within_state_graph(ccs, perturbation_col = perturbation_col,
      control_ids = "ctrl-inj", perturbations = "NULLPERT",
      nuisance_model_formula_str = "0", cell_groups = ok, cores = cores,
      write_dir = wd, max_simultaneous_genes = max_simultaneous_genes,
      filter_mode = "by_background_count", detection_min_samples = 1,
      background_bottom_frac = 0.25, background_quantile_p = 0.99,
      background_count_floor = 2, nperm = nperm)
    .read_null_dir(wd) %>% dplyr::mutate(log_ratio = log2(N / (n_ctrl - N)))
  })
}

.read_null_dir <- function(wd) {
  files <- list.files(wd, "_within_node_degs.csv$", full.names = TRUE)
  if (!length(files)) return(tibble::tibble())
  purrr::map_dfr(files, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE))) %>%
    dplyr::filter(is.finite(.data$perturb_to_ctrl_shrunken_lfc),
                  .data$perturb_to_ctrl_shrunken_lfc_se > 0,
                  is.finite(.data$log_mean_expression)) %>%
    dplyr::transmute(.data$cell_group, .data$log_mean_expression,
                     z = .data$perturb_to_ctrl_shrunken_lfc / .data$perturb_to_ctrl_shrunken_lfc_se) %>%
    dplyr::filter(is.finite(.data$z))
}

#' Train the empirical-FDR null model from control-split null draws.
#' @param null Output of [build_empirical_null()] (log_mean_expression, z, log_ratio).
#' @param expr_breaks Expression bin edges.
#' @return An `efdr_model` object (a tibble of per-stratum sorted |z| references).
#' @export
train_efdr_model <- function(null, expr_breaks = .EFDR_BREAKS) {
  ref <- null %>%
    dplyr::filter(is.finite(.data$log_mean_expression), is.finite(.data$z)) %>%
    dplyr::mutate(eb = cut(.data$log_mean_expression, expr_breaks),
                  dir = dplyr::if_else(.data$z < 0, "dn", "up"),
                  a = abs(.data$z)) %>%
    dplyr::group_by(.data$eb, .data$dir, .data$log_ratio) %>%
    dplyr::summarise(ref = list(sort(.data$a)), .groups = "drop")
  structure(list(ref = ref, ratios = sort(unique(ref$log_ratio)), expr_breaks = expr_breaks),
            class = "efdr_model")
}

# vectorised right-tail probability of |z| values `a` against a sorted reference
.efdr_tailp <- function(a, ref) {
  if (is.null(ref) || !length(ref)) return(rep(0.5, length(a)))
  (length(ref) - findInterval(a, ref, left.open = TRUE) + 1) / (length(ref) + 1)
}

#' Annotate a DEG table with per-gene empirical p and BH FDR at a given log-ratio.
#'
#' @param model An `efdr_model` from [train_efdr_model()].
#' @param deg_tbl DEG table with `log_mean_expression` and either `z` or
#'   `perturb_to_ctrl_shrunken_lfc` + `perturb_to_ctrl_shrunken_lfc_se`.
#' @param log_ratio The perturbation's log2(n_perturb / n_control).
#' @param group_col Column to BH-adjust within (default `cell_group`).
#' @return `deg_tbl` with added `empirical_p` and `empirical_fdr`.
#' @export
annotate_empirical_fdr_model <- function(model, deg_tbl, log_ratio, group_col = "cell_group") {
  R <- model$ratios
  br <- if (log_ratio <= R[1]) list(lo = R[1], hi = R[1], w = 1) else
        if (log_ratio >= R[length(R)]) list(lo = R[length(R)], hi = R[length(R)], w = 1) else {
          hi <- R[which(R >= log_ratio)[1]]; lo <- R[max(which(R <= log_ratio))]
          list(lo = lo, hi = hi, w = if (hi == lo) 1 else (hi - log_ratio) / (hi - lo))
        }
  ref_at <- function(r) model$ref %>% dplyr::filter(.data$log_ratio == r) %>%
    { purrr::set_names(.$ref, paste(.$eb, .$dir)) }
  rlo <- ref_at(br$lo); rhi <- if (br$hi == br$lo) rlo else ref_at(br$hi)

  deg_tbl %>%
    dplyr::mutate(
      z = if ("z" %in% names(.)) .data$z
          else .data$perturb_to_ctrl_shrunken_lfc / .data$perturb_to_ctrl_shrunken_lfc_se,
      eb = cut(.data$log_mean_expression, model$expr_breaks),
      dir = dplyr::if_else(.data$z < 0, "dn", "up"),
      .a = abs(.data$z)) %>%
    dplyr::group_by(.data$eb, .data$dir) %>%
    # one findInterval per (eb,dir) group over the whole group's |z| -- fully vectorised
    dplyr::mutate(empirical_p = {
      key <- paste(dplyr::cur_group()$eb, dplyr::cur_group()$dir)
      br$w * .efdr_tailp(.data$.a, rlo[[key]]) +
        (1 - br$w) * .efdr_tailp(.data$.a, rhi[[key]])
    }) %>%
    dplyr::group_by(dplyr::across(all_of(group_col))) %>%
    dplyr::mutate(empirical_fdr = stats::p.adjust(.data$empirical_p, "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-"eb", -"dir", -".a")
}
