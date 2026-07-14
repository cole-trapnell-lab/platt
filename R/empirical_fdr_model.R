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

.EFDR_BREAKS <- c(-Inf, seq(-8, 0.5, 0.5), Inf)   # retained for back-compat / bin midpoints
# The null tail is a rare, upper-tail phenomenon (the lower quantiles are a flat,
# shot-noise-dominated bulk), so the tau grid is concentrated above 0.9 with a
# couple of low anchors for interpolation.
.EFDR_TAUS   <- c(0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)
.EFDR_EGRID  <- seq(-9, 1, 0.25)      # expression grid the tail surface is evaluated on
.EFDR_PEGRID <- seq(0, 100, 10)       # %embryos grid (control detection rate)

# --- The empirical-null policy (one set of baked constants; not pipeline-configurable) ---
.EFDR_C_FULL <- 10    # cells an embryo must carry to count as ONE full perturbation replicate;
                      # below this it contributes fractionally (cells / c_full). Drives the
                      # cell-count-aware degrees-of-freedom of the within-state contrast.
.EFDR_MIN_PB <- 2L    # a cell type needs at least this many perturbation embryos to have a
                      # valid two-group contrast; below it, calls are gated (empirical_p = 1).

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
    # Per-pseudo-group detection for THIS split: pct_emb = ctrl-inj (losing) arm,
    # pct_emb_pert = NULLPERT (gaining) arm. The up-null is conditioned on the
    # gaining-arm detection so the few-cell blip distribution is captured; the
    # down-null on the losing (control) arm, as before.
    det <- efdr_detection_rate(sub, sample_group = sample_group, cell_group = cell_group,
                               control_ids = "ctrl-inj", perturbation_col = perturbation_col)
    # Per-split perturbation-arm replicate support (NULLPERT pseudobulks per cell type):
    # carried so the model can calibrate the honest degrees-of-freedom vs support (the
    # trend-dispersion shrinkage buys back df, so the effective df exceeds n_pert_pb - 1;
    # the null is the reference for how much).
    supp <- efdr_perturbation_support(sub, sample_group = sample_group, cell_group = cell_group,
                                      control_ids = "ctrl-inj", perturbation_col = perturbation_col)
    # supp carries per-cell-type pseudo-arm cell counts (n_ctrl_cells, n_pert_cells), so the
    # null z is SE-floored on the same footing as the real DEG.
    .read_null_dir(wd, ct_cells = supp[, c("cell_group", "n_ctrl_cells", "n_pert_cells")]) %>%
      dplyr::left_join(det, by = c("gene_short_name", "cell_group")) %>%
      dplyr::left_join(supp[, c("cell_group", "n_pert_pb", "n_eff_pb")], by = "cell_group") %>%
      dplyr::mutate(log_ratio = log2(N / (n_ctrl - N)))
  })
}

# Theoretical lower bound on the shrunken-LFC standard error. The trend-dispersion
# shrinkage can manufacture an insane |z| on near-zero-count genes (a gene detected in <1%
# of cells getting a tighter SE than a broadly-expressed one), and because the control-split
# null is built from the same GLM its tail inherits the same insanity -- so the empirical
# model can only demote so far. This caps the SE at the information the counts can actually
# support: recover each arm's mean expression from the pooled mean, the fold change and the
# arm cell counts, then apply the Poisson/delta-method floor se >= sqrt(1/K_ctrl + 1/K_pert)
# on the effective counts K = N * mu. Well-powered calls are untouched (their real SE already
# exceeds the floor); it is applied identically to the null and the real DEG so the statistic
# and the distribution it is judged against stay on the same footing. LFC is natural-log.
.efdr_se_floor <- function(log_mean_expression, lfc, n_ctrl, n_pert, kappa = 1) {
  mu <- exp(log_mean_expression)
  r  <- exp(lfc)                                        # mu_pert / mu_ctrl
  denom <- n_ctrl + n_pert * r
  mu_ctrl <- ifelse(denom > 0, mu * (n_ctrl + n_pert) / denom, mu)
  K_ctrl <- pmax(n_ctrl * mu_ctrl, 0)
  K_pert <- pmax(n_pert * mu_ctrl * r, 0)
  sqrt(1 / (K_ctrl + kappa) + 1 / (K_pert + kappa))
}

# Apply the SE floor to a table carrying lfc/se/log_mean_expression + per-cell-type arm cell
# counts (n_ctrl_cells, n_pert_cells), returning the floored standardized statistic.
.efdr_floored_z <- function(lfc, se, log_mean_expression, n_ctrl_cells, n_pert_cells) {
  ok <- is.finite(n_ctrl_cells) & is.finite(n_pert_cells)
  se_fl <- se
  se_fl[ok] <- pmax(se[ok], .efdr_se_floor(log_mean_expression[ok], lfc[ok],
                                           n_ctrl_cells[ok], n_pert_cells[ok]))
  lfc / se_fl
}

.read_null_dir <- function(wd, ct_cells = NULL) {
  files <- list.files(wd, "_within_node_degs.csv$", full.names = TRUE)
  if (!length(files)) return(tibble::tibble())
  d <- purrr::map_dfr(files, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE))) %>%
    dplyr::filter(is.finite(.data$perturb_to_ctrl_shrunken_lfc),
                  .data$perturb_to_ctrl_shrunken_lfc_se > 0,
                  is.finite(.data$log_mean_expression))
  # Floor the SE with the split's per-cell-type pseudo-arm cell counts (same treatment the
  # real DEG gets in annotate); fall back to the raw z if counts are unavailable.
  if (!is.null(ct_cells)) {
    d <- dplyr::left_join(d, ct_cells, by = "cell_group")
    z <- .efdr_floored_z(d$perturb_to_ctrl_shrunken_lfc, d$perturb_to_ctrl_shrunken_lfc_se,
                         d$log_mean_expression, d$n_ctrl_cells, d$n_pert_cells)
  } else {
    z <- d$perturb_to_ctrl_shrunken_lfc / d$perturb_to_ctrl_shrunken_lfc_se
  }
  tibble::tibble(gene_short_name = d$gene_short_name, cell_group = d$cell_group,
                 log_mean_expression = d$log_mean_expression, z = z) %>%
    dplyr::filter(is.finite(.data$z))
}

#' Per-(gene, cell-type) control detection rate (%embryos with detectable expression).
#'
#' The empirical null's second covariate. A gene carried by few embryos has
#' between-embryo variance the mean-dispersion trend under-prices (the DESeq-style
#' shrinkage used in the DEG fit), which inflates |z|; %embryos is a scale-free,
#' per-experiment proxy for that dispersion excess. Computed from control cells
#' only (the null is control-split), so the same table applies to the null draws
#' and to the real DEG table.
#'
#' @param cds cell_data_set (control cells are selected internally).
#' @param sample_group,cell_group,control_ids,perturbation_col colData columns/values.
#' @return tibble(gene_short_name, cell_group, pct_emb) with pct_emb in [0, 100].
#' @export
efdr_detection_rate <- function(cds,
                                sample_group = "embryo_ID",
                                cell_group = "cell_type",
                                control_ids = c("ctrl-inj"),
                                perturbation_col = "perturbation") {
  # Per-(gene, cell) detection breadth (% of embryos with >=1 cell detecting), computed
  # SEPARATELY for the control and perturbation arms. The two are used direction-specifically
  # downstream: control detection conditions the DOWN null (on-in-control artifact), perturbation
  # detection conditions the UP null (off-in-control, few-perturbation-cell blip artifact). They
  # are never used together, avoiding colinearity between the two (highly correlated) metrics.
  .rate <- function(cds_sub, col_name) {
    cd <- SummarizedExperiment::colData(cds_sub)
    M <- monocle3::exprs(cds_sub)
    gsym <- SummarizedExperiment::rowData(cds_sub)$gene_short_name
    if (is.null(gsym)) gsym <- rownames(cds_sub)
    ct <- as.character(cd[[cell_group]]); emb <- as.character(cd[[sample_group]])
    out <- purrr::map_dfr(unique(ct), function(g) {
      idx <- which(ct == g); if (length(idx) < 2L) return(NULL)
      e <- factor(emb[idx])
      B <- methods::as(M[, idx, drop = FALSE], "dgCMatrix"); B@x[] <- 1  # binarise detection (BPCells -> sparse)
      E <- Matrix::sparseMatrix(i = seq_along(e), j = as.integer(e), x = 1,
                                dims = c(length(e), nlevels(e)))
      detc <- as.matrix(B %*% E)
      tibble::tibble(gene_short_name = gsym, cell_group = g,
                     rate = 100 * Matrix::rowMeans(detc > 0))
    })
    if (is.null(out) || nrow(out) == 0) {
      return(tibble::tibble(gene_short_name = character(), cell_group = character()) %>%
               dplyr::mutate(!!col_name := numeric()))
    }
    out %>%
      dplyr::group_by(.data$gene_short_name, .data$cell_group) %>%
      dplyr::summarise(rate = max(.data$rate), .groups = "drop") %>%
      dplyr::rename(!!col_name := "rate")
  }
  cd0 <- SummarizedExperiment::colData(cds)
  is_ctrl <- as.character(cd0[[perturbation_col]]) %in% control_ids
  ctrl <- .rate(cds[, is_ctrl, drop = FALSE], "pct_emb")                 # control detection (kept name for back-compat)
  pert <- .rate(cds[, !is_ctrl, drop = FALSE], "pct_emb_pert")           # perturbation detection
  dplyr::full_join(ctrl, pert, by = c("gene_short_name", "cell_group")) %>%
    dplyr::mutate(pct_emb = ifelse(is.na(.data$pct_emb), 0, .data$pct_emb),
                  pct_emb_pert = ifelse(is.na(.data$pct_emb_pert), 0, .data$pct_emb_pert))
}

#' Per-cell-type perturbation-arm replicate support.
#'
#' Counts, per cell type, the number of **perturbation** embryos (pseudobulks)
#' contributing at least `min_cells` cells to the within-state contrast. This is
#' the replicate count the pseudobulk GLM actually has on the perturbation side,
#' and it is the honest degrees-of-freedom for the contrast: when a cell type is
#' depleted in the perturbation, the deep control arm supplies a small standard
#' error as if the contrast were well-powered, inflating |z| and manufacturing
#' spurious (predominantly down) calls. [annotate_empirical_fdr_model()] uses this
#' to floor the ordinary p by a t-distribution with `n_pert_pb - 1` df, and to
#' gate cell types with <2 perturbation pseudobulks (no valid two-group contrast).
#'
#' @param cds contrast cell_data_set.
#' @param sample_group,cell_group,control_ids,perturbation_col colData columns/values.
#' @param min_cells Minimum cells for an embryo to count as a replicate. Default 1:
#'   the replicate unit is the embryo (an embryo either has the cell type or not);
#'   cells-per-embryo feeds the per-pseudobulk variance, not the replicate count, and
#'   a higher bar would wrongly drop cell types that are present across many embryos
#'   but sparse per embryo (common in shallow captures). Matches the within-state DEG,
#'   which fits on any embryo with cells (`min_cells_per_pseudobulk = NULL`).
#' @return tibble(cell_group, n_pert_pb, n_ctrl_pb, n_eff_pb, n_pert_cells, n_ctrl_cells).
#' @export
efdr_perturbation_support <- function(cds,
                                      sample_group = "embryo_ID",
                                      cell_group = "cell_type",
                                      control_ids = c("ctrl-inj"),
                                      perturbation_col = "perturbation") {
  cd <- SummarizedExperiment::colData(cds)
  efdr_support_from_coldata(
    data.frame(emb  = as.character(cd[[sample_group]]),
               pert = as.character(cd[[perturbation_col]]),
               ct   = as.character(cd[[cell_group]]),
               stringsAsFactors = FALSE),
    sample_group = "emb", cell_group = "ct",
    control_ids = control_ids, perturbation_col = "pert")
}

#' Per-cell-type perturbation-arm support, from a plain coldata table.
#'
#' The single source of truth for every support quantity the empirical-FDR model needs,
#' computed from a (embryo, perturbation, cell_type) table so it can be driven from either a
#' cell_data_set's colData ([efdr_perturbation_support()]) or a lightweight coldata TSV (the
#' mcclintock efdr stage) without duplicating the definitions across repos:
#'   - n_pert_pb   : perturbation embryos with the cell type (the gate's replicate count).
#'   - n_ctrl_pb   : control embryos with the cell type.
#'   - n_eff_pb    : cell-weighted effective perturbation replicates, sum min(1, cells/C_FULL)
#'                   -- an embryo is a full replicate only once it carries `.EFDR_C_FULL`
#'                   cells; drives the degrees-of-freedom correction.
#'   - n_pert_cells / n_ctrl_cells : total cells per arm, feeding the theoretical SE floor.
#' @param coldata data.frame/table with embryo, perturbation and cell-type columns.
#' @param sample_group,cell_group,control_ids,perturbation_col column names / control labels.
#' @export
efdr_support_from_coldata <- function(coldata,
                                      sample_group = "embryo_ID",
                                      cell_group = "cell_type",
                                      control_ids = c("ctrl-inj"),
                                      perturbation_col = "perturbation") {
  df <- tibble::tibble(ct = as.character(coldata[[cell_group]]),
                       emb = as.character(coldata[[sample_group]]),
                       is_ctrl = as.character(coldata[[perturbation_col]]) %in% control_ids)
  df %>%
    dplyr::count(.data$ct, .data$emb, .data$is_ctrl, name = "n") %>%
    dplyr::group_by(.data$ct) %>%
    dplyr::summarise(
      n_pert_pb    = dplyr::n_distinct(.data$emb[!.data$is_ctrl]),
      n_ctrl_pb    = dplyr::n_distinct(.data$emb[.data$is_ctrl]),
      n_eff_pb     = sum(pmin(1, .data$n[!.data$is_ctrl] / .EFDR_C_FULL)),
      n_pert_cells = sum(.data$n[!.data$is_ctrl]),
      n_ctrl_cells = sum(.data$n[.data$is_ctrl]),
      .groups = "drop") %>%
    dplyr::rename(cell_group = "ct")
}

#' Train the empirical-FDR null model from control-split null draws.
#'
#' Models the null upper tail of |z| as a smooth function of **expression and
#' control %embryos**, per (direction, sampling log-ratio):
#'   - **monotone-decreasing in expression** (scam `bs="mpd"`): the dispersion
#'     mis-pricing that inflates |z| worsens monotonically as expression -> 0, so
#'     the tail must not droop back toward zero in the sparse extreme-low region.
#'   - **decreasing in %embryos** (scam `bs="mpd"`): broadly-detected genes have
#'     well-estimated dispersion (light tail); narrowly-detected genes carry the
#'     under-priced between-embryo variance (heavy tail).
#' Direction is kept separate (the control-depth artifact is one-directional).
#'
#' @param null Output of [build_empirical_null()] (gene_short_name, cell_group,
#'   log_mean_expression, z, log_ratio).
#' @param detection Output of [efdr_detection_rate()] (gene_short_name, cell_group,
#'   pct_emb); joined to the null to supply the %embryos covariate.
#' @param taus Quantile levels (tail-concentrated) the surface is evaluated at.
#' @param egrid,pegrid Expression and %embryos grids the surface is evaluated on.
#' @param bin_width Expression bin width for the per-cell empirical quantiles.
#' @param pe_bin %embryos bin width for the per-cell empirical quantiles.
#' @return An `efdr_model`: per-(direction, log-ratio) 2-D tail surfaces on the
#'   (egrid x pegrid) grid, one column per tau.
#' @export
train_efdr_model <- function(null, detection = NULL, taus = .EFDR_TAUS,
                             egrid = .EFDR_EGRID, pegrid = .EFDR_PEGRID,
                             bin_width = 0.5, pe_bin = 10) {
  # The null carries per-pseudo-group detection (pct_emb = losing/control arm,
  # pct_emb_pert = gaining/perturbation arm). Prefer those; only join a supplied
  # `detection` for columns the null lacks (back-compat with older nulls).
  if (!is.null(detection)) {
    join_cols <- setdiff(intersect(c("pct_emb", "pct_emb_pert"), names(detection)), names(null))
    if (length(join_cols)) {
      null <- dplyr::left_join(null, detection[, c("gene_short_name", "cell_group", join_cols)],
                               by = c("gene_short_name", "cell_group"))
    }
  }
  null <- null %>%
    dplyr::filter(is.finite(.data$log_mean_expression), is.finite(.data$z)) %>%
    dplyr::mutate(dir = dplyr::if_else(.data$z < 0, "dn", "up"), a = abs(.data$z))
  ratios <- sort(unique(null$log_ratio))
  GR <- expand.grid(eb = egrid, pe = pegrid)          # eb varies fastest
  surf <- list()
  for (r in ratios) for (d in c("dn", "up")) {
    key <- paste(d, r)
    # Direction-specific detection covariate: the arm that "has" the gene in a
    # call of this direction. up -> gaining (perturbation) arm; dn -> control arm.
    sub <- null %>% dplyr::filter(.data$log_ratio == r, .data$dir == d)
    sub$pct_emb <- if (d == "up") sub$pct_emb_pert else sub$pct_emb
    use_pe <- mean(is.finite(sub$pct_emb)) > 0.5 &&
              length(unique(round(sub$pct_emb[is.finite(sub$pct_emb)] / pe_bin))) >= 3L
    pb <- sub %>%
      dplyr::mutate(eb = round(.data$log_mean_expression / bin_width) * bin_width,
                    pe = round(.data$pct_emb / pe_bin) * pe_bin) %>%
      { if (use_pe) dplyr::group_by(., .data$eb, .data$pe) else dplyr::group_by(., .data$eb) } %>%
      dplyr::summarise(q = list(stats::quantile(.data$a, taus, names = FALSE)),
                       n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(.data$n >= 30L)
    if (nrow(pb) < 8L) { surf[[key]] <- NULL; next }   # too sparse -> caller defers to ashr
    qmat <- do.call(rbind, pb$q)                       # cells x taus
    Q <- vapply(seq_along(taus), function(i) {
      if (use_pe) {
        df <- data.frame(q = qmat[, i], eb = pb$eb, pe = pb$pe, n = pb$n)
        m <- tryCatch(
          scam::scam(q ~ s(eb, bs = "mpd") + s(pe, bs = "mpd"), weights = n, data = df),
          error = function(e1) tryCatch(
            scam::scam(q ~ s(eb, bs = "mpd"), weights = n, data = df),
            error = function(e2) mgcv::gam(q ~ s(eb), weights = n, data = df)))
        pmax(as.numeric(stats::predict(m, newdata = GR)), 0)
      } else {
        # no usable %embryos -> 1-D monotone in expression, tiled across the pe grid
        df <- data.frame(q = qmat[, i], eb = pb$eb, n = pb$n)
        m <- tryCatch(scam::scam(q ~ s(eb, bs = "mpd"), weights = n, data = df),
                      error = function(e2) mgcv::gam(q ~ s(eb), weights = n, data = df))
        v <- pmax(as.numeric(stats::predict(m, newdata = data.frame(eb = egrid))), 0)
        rep(v, times = length(pegrid))                 # eb fastest in GR -> tile per pe block
      }
    }, numeric(nrow(GR)))                               # nrow(GR) x taus
    Q <- t(apply(Q, 1, cummax))                        # monotone in tau within each grid cell
    surf[[key]] <- Q
  }
  df_calib <- .efdr_calibrate_df(null)
  if (is.null(df_calib))
    stop("efdr: could not calibrate degrees-of-freedom from the null (too few draws). ",
         "A valid efdr_model requires the support->df calibration.", call. = FALSE)
  structure(list(surf = surf, taus = taus, egrid = egrid, pegrid = pegrid, ratios = ratios,
                 df_calib = df_calib),
            class = "efdr_model")
}

# Calibrate the honest effective degrees-of-freedom of the within-state contrast as a
# function of the cell-weighted effective perturbation support (n_eff_pb). Naive Welch says
# the df is (replicates - 1), but the trend-dispersion shrinkage borrows dispersion across
# genes and buys real df back, so the effective df is higher. The null is the reference: at
# moderate expression (isolating the sample-size effect from the low-expression artifact the
# surface already prices), match the null |z| right-tail exceedance at `z0` to a
# t-distribution and solve for its df. Returns a monotone (support -> df) lookup; NULL only
# if the null has too few draws to calibrate (train_efdr_model then errors).
.efdr_calibrate_df <- function(null, z0 = 3.5, min_n = 500L) {
  d <- null[is.finite(null$n_eff_pb) & is.finite(null$z) &
              is.finite(null$log_mean_expression) & null$log_mean_expression > -1, ]
  if (nrow(d) < min_n) return(NULL)
  a <- abs(d$z); k <- d$n_eff_pb
  kb <- round(k)                                        # bin continuous effective support to integer levels
  solve_df <- function(p) {                              # 2*pt(-z0, df) is decreasing in df
    if (!is.finite(p) || p <= 2 * stats::pnorm(-z0)) return(Inf)  # normal already conservative
    g <- function(l) 2 * stats::pt(-z0, exp(l)) - p
    if (g(log(0.5)) < 0) return(0.5)
    if (g(log(1e4)) > 0) return(Inf)
    exp(stats::uniroot(g, c(log(0.5), log(1e4)))$root)
  }
  levs <- sort(unique(kb[kb >= 2]))
  rows <- lapply(levs, function(kk) {
    ak <- a[kb == kk]; if (length(ak) < min_n) return(NULL)
    data.frame(support = kk, df_eff = solve_df(mean(ak >= z0)), n = length(ak))
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) < 2L) return(NULL)
  # Enforce df non-decreasing in support with isotonic regression (pool-adjacent-violators),
  # weighted by null mass; Inf (no correction) is capped at `df_cap` for the fit, above which
  # the t-distribution is indistinguishable from normal so the correction is negligible.
  df_cap <- 60
  y <- pmin(out$df_eff, df_cap)
  out$df_eff <- stats::isoreg(out$support, y)$yf
  out[, c("support", "df_eff")]
}

# right-tail probability of |z| = `a` at (expression `e`, %embryos `pe`) against a
# 2-D tail surface. NA where the surface is absent (too-sparse stratum) or pe is
# missing -> caller defers to the ordinary test.
.efdr_tail_2d <- function(Q, taus, egrid, pegrid, e, pe, a) {
  if (is.null(Q)) return(rep(NA_real_, length(a)))
  ne <- length(egrid); nt <- length(taus)
  ie <- pmin(pmax(findInterval(e, egrid), 1L), ne)
  ip <- pmin(pmax(findInterval(pe, pegrid), 1L), length(pegrid))
  ri <- (ip - 1L) * ne + ie                            # row in GR (eb fastest)
  vapply(seq_along(a), function(k) {
    if (is.na(ri[k]) || is.na(a[k])) return(NA_real_)
    qi <- Q[ri[k], ]; ak <- a[k]
    if (ak <= qi[1])  return(1 - taus[1] * 0.5)
    if (ak >= qi[nt]) return(1 - taus[nt])
    j <- max(which(qi <= ak))
    1 - (taus[j] + (ak - qi[j]) / (qi[j + 1] - qi[j]) * (taus[j + 1] - taus[j]))
  }, numeric(1))
}

#' Annotate a DEG table with per-gene empirical p and BH FDR at a given log-ratio.
#'
#' @param model An `efdr_model` from [train_efdr_model()].
#' @param deg_tbl DEG table with `gene_short_name`, `cell_group`,
#'   `log_mean_expression`, and either `z` or `perturb_to_ctrl_shrunken_lfc` +
#'   `perturb_to_ctrl_shrunken_lfc_se`.
#' @param detection Output of [efdr_detection_rate()] (gene_short_name, cell_group,
#'   pct_emb); supplies the %embryos covariate. Genes with no match fall back to
#'   the ordinary test (still floored below).
#' @param log_ratio The perturbation's log2(n_perturb / n_control).
#' @param group_col Column to BH-adjust within (default `cell_group`).
#' @return `deg_tbl` with added `empirical_p` and `empirical_fdr`.
#' @details Two coupled pieces:
#'   * The 2-D null tail (monotone in expression, smooth in %embryos) *demotes*
#'     the dispersion-under-priced artifact calls (narrow / low-expression).
#'   * An **ashr floor** — `empirical_p = max(tail, perturb_to_ctrl_p_value)` — is
#'     applied **always** (including where the tail is unavailable). It is
#'     load-bearing: without it, conditioning the null tighter lets trivially small
#'     |z| pass. So the empirical step may only *demote* a call ashr made, never
#'     promote one.
#' @export
annotate_empirical_fdr_model <- function(model, deg_tbl, log_ratio, detection = NULL,
                                         support = NULL, min_pseudobulks = .EFDR_MIN_PB,
                                         group_col = "cell_group") {
  R <- model$ratios
  br <- if (log_ratio <= R[1]) list(lo = R[1], hi = R[1], w = 1) else
        if (log_ratio >= R[length(R)]) list(lo = R[length(R)], hi = R[length(R)], w = 1) else {
          hi <- R[which(R >= log_ratio)[1]]; lo <- R[max(which(R <= log_ratio))]
          list(lo = lo, hi = hi, w = if (hi == lo) 1 else (hi - log_ratio) / (hi - lo))
        }
  # Drop any pre-existing detection/support columns before (re)joining: decorating an
  # already-decorated DEG table (e.g. re-running the efdr stage) would otherwise collide
  # on `pct_emb`/`pct_emb_pert`/`n_pert_pb`, producing `.x`/`.y` and leaving no plain
  # `pct_emb` -- which silently NA's the covariate and skips the tail. Re-derive them here.
  if (!is.null(detection)) {
    deg_tbl <- deg_tbl[, setdiff(names(deg_tbl), c("pct_emb", "pct_emb_pert")), drop = FALSE]
    deg_tbl <- dplyr::left_join(deg_tbl, detection, by = c("gene_short_name", "cell_group"))
  }
  if (!is.null(support)) {
    supp_cols <- intersect(c("n_pert_pb", "n_eff_pb", "n_ctrl_cells", "n_pert_cells"), names(support))
    deg_tbl <- deg_tbl[, setdiff(names(deg_tbl), supp_cols), drop = FALSE]
    deg_tbl <- dplyr::left_join(deg_tbl, support[, c("cell_group", supp_cols)], by = "cell_group")
  }
  if (!"n_ctrl_cells" %in% names(deg_tbl)) deg_tbl$n_ctrl_cells <- NA_real_
  if (!"n_pert_cells" %in% names(deg_tbl)) deg_tbl$n_pert_cells <- NA_real_
  # Map effective perturbation-arm support (n_eff_pb) -> honest effective df. Within the
  # range the null covers, use the null-calibrated (support -> df) lookup (accounts for the
  # shrinkage's df buy-back). BELOW the smallest calibrated support the control-split null
  # has no coverage (its cell types are well-populated), so df is ramped analytically from a
  # 0.5 floor at <=1 effective replicate up to the smallest calibrated df -- this is what
  # demotes cell types that are thin *in cells* (n_eff -> 0), whose over-shrunk SE would
  # otherwise manufacture a confident call.
  cal <- model$df_calib
  s0 <- cal$support[1]; d0 <- cal$df_eff[1]              # smallest calibrated support and its df
  dfmap <- function(n) {
    out <- rep(NA_real_, length(n)); ok <- is.finite(n); nn <- n[ok]
    idx <- findInterval(nn, cal$support)                 # within/above range: step lookup
    val <- cal$df_eff[pmin(pmax(idx, 1L), nrow(cal))]
    below <- nn < s0                                     # extend below the calibrated floor
    val[below] <- if (s0 > 1) pmax(0.5, 0.5 + (d0 - 0.5) * (nn[below] - 1) / (s0 - 1)) else 0.5
    out[ok] <- val
    out
  }

  out <- deg_tbl %>%
    dplyr::mutate(
      # SE-floored standardized statistic: cap the shrunken SE at the Poisson floor the counts
      # support, so a near-zero gene can't get an insane |z| (same floor the null z got, so
      # they stay on the same footing). Falls back to the raw z only if lfc/se/expression are
      # absent (a caller that passed only `z`).
      z = if (all(c("perturb_to_ctrl_shrunken_lfc", "perturb_to_ctrl_shrunken_lfc_se",
                    "log_mean_expression") %in% names(.)))
            .efdr_floored_z(.data$perturb_to_ctrl_shrunken_lfc, .data$perturb_to_ctrl_shrunken_lfc_se,
                            .data$log_mean_expression, .data$n_ctrl_cells, .data$n_pert_cells)
          else if ("z" %in% names(.)) .data$z
          else .data$perturb_to_ctrl_shrunken_lfc / .data$perturb_to_ctrl_shrunken_lfc_se,
      .dir = dplyr::if_else(.data$z < 0, "dn", "up"),
      .a = abs(.data$z),
      # Ordinary reference p (the floor): the trusted ashr shrunken p if present, else
      # normal(z). When perturbation-arm support is known, the ordinary test is made honest
      # about the KO-arm degrees of freedom: a normal/ashr p assumes ~infinite df, but a
      # contrast backed by few perturbation pseudobulks has finite (null-calibrated) df, so
      # |z| is re-scored against a t-distribution. This demotes the cell-abundance confound
      # (depleted cell types where the deep control arm lends unearned power) at the same
      # decorate-in-place layer, with no DEG re-fit; it is a pure function of KO support, so
      # it does not touch expression and cannot demote well-supported calls.
      .df = dfmap(.data$n_eff_pb),
      .p_df = dplyr::if_else(is.finite(.data$.df),
                             2 * stats::pt(-.data$.a, df = pmax(.data$.df, 0.5)), NA_real_),
      .ord_p = if ("perturb_to_ctrl_p_value" %in% names(.)) .data$perturb_to_ctrl_p_value
               else 2 * (1 - stats::pnorm(.data$.a)),
      # df-honesty is demote-only: never let the ordinary p be more confident than the t-test
      .ord_p = dplyr::if_else(is.na(.data$.p_df), .data$.ord_p, pmax(.data$.ord_p, .data$.p_df)))

  # 2-D empirical tail, interpolated across the two bracketing log-ratios.
  emp_tail <- rep(NA_real_, nrow(out))
  for (d in c("dn", "up")) {
    idx <- which(out$.dir == d); if (!length(idx)) next
    # Direction-specific detection: up-calls scored against the perturbation-arm
    # (gaining) detection surface; down-calls against the control-arm surface.
    pe_src <- if (d == "up") out$pct_emb_pert else out$pct_emb
    e <- out$log_mean_expression[idx]; pe <- pe_src[idx]; a <- out$.a[idx]
    tl <- .efdr_tail_2d(model$surf[[paste(d, br$lo)]], model$taus, model$egrid, model$pegrid, e, pe, a)
    th <- if (br$hi == br$lo) tl else
          .efdr_tail_2d(model$surf[[paste(d, br$hi)]], model$taus, model$egrid, model$pegrid, e, pe, a)
    emp_tail[idx] <- dplyr::coalesce(br$w * tl + (1 - br$w) * th, tl, th)
  }
  # ashr floor, ALWAYS: demote-only. Where the tail is NA (missing surface or pct_emb),
  # empirical_p is the ordinary p, so the floor still holds (no fallback hole).
  emp_tail <- pmin(pmax(emp_tail, 0), 1)
  out$empirical_p <- dplyr::if_else(is.na(emp_tail), out$.ord_p, pmax(emp_tail, out$.ord_p))
  # Support gate: a cell type with fewer than `min_pseudobulks` perturbation pseudobulks
  # has no valid two-group contrast (df < 1); its calls are uncallable, not significant.
  gate <- is.finite(out$n_pert_pb) & out$n_pert_pb < min_pseudobulks
  out$empirical_p[gate] <- 1

  out %>%
    dplyr::group_by(dplyr::across(all_of(group_col))) %>%
    dplyr::mutate(empirical_fdr = stats::p.adjust(.data$empirical_p, "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-".dir", -".a", -".ord_p", -".df", -".p_df")
}
