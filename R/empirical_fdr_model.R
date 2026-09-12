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
.EFDR_TAUS   <- c(0.5, 0.75, 0.85, 0.9, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999)
                                      # densified around the demotion boundary (~0.95-0.98) so the
                                      # tail p interpolates finely where the call/demote decision sits
.EFDR_EGRID  <- seq(-9, 1, 0.25)      # expression grid the tail surface is evaluated on
.EFDR_PEGRID <- seq(0, 100, 10)       # %embryos grid (legacy detection axis; unused by the arm-matched surface)
.EFDR_LCGRID <- seq(0, 3.6, 0.15)     # log10(thin-arm cell count + 1) grid the tail is conditioned on
                                      # (~1 cell -> ~4000 cells); the arm-size axis that decides
                                      # whether a large |z| is a real depletion (deep arm) or a
                                      # thin-arm sampling artifact.

# --- The empirical-null policy (one set of baked constants; not pipeline-configurable) ---
.EFDR_C_FULL <- 10    # cells an embryo must carry to count as ONE full perturbation replicate;
                      # below this it contributes fractionally (cells / c_full). Drives the
                      # cell-count-aware degrees-of-freedom of the within-state contrast.
.EFDR_MIN_PB <- 2L    # a cell type needs at least this many perturbation embryos to have a
                      # valid two-group contrast; below it, calls are gated (empirical_p = 1).
.EFDR_MIN_ARM_UMI <- 1  # a (gene, cell type) needs >= this many total UMIs in AT LEAST ONE arm
                      # to be callable; with both arms empty the gene is absent from the cell
                      # type and its |z| is shrinkage noise (gated). Tests max(K_ctrl, K_pert),
                      # so a well-counted control arm + empty perturbation arm (real complete
                      # depletion) is KEPT -- only genuinely-absent genes are gated.
.EFDR_MIN_SRC_DET <- 0.02  # a LOSS (down call) needs the CONTROL (source) arm to detect the gene in
                      # at least this fraction of its cells to be callable. A depletion presupposes a
                      # robustly-expressed baseline; a "loss" of a gene detected in <2% of control
                      # cells is trend-dispersion shrinkage noise (huge |z| on near-zero counts) the
                      # control-split null cannot reach at moderate arm sizes. Data-derived, not an
                      # a-priori expression floor: real depletions sit at >=12% control detection,
                      # Amy's confirmed 'too low' down-artifacts all at <=1.1%, with a clean ~10x gap;
                      # 2% sits in that gap. DIRECTIONAL (control arm, down only) -- NOT applied to up
                      # calls, because real GAINS can be sparse (phox2a is a real up-call at 0.7%
                      # perturbation detection), so no analogous gap exists on the up side.
                      # OFF by default; see .EFDR_DEFAULT_GATES.
# --- Is the thinner arm powered to see THIS gene? (per call, both directions) -----------------
# A reduced-or-absent signal in the thinner arm is only callable if that arm held enough cells to
# have detected the gene in the first place. If a gene is detected in a fraction p of cells, an arm
# of n cells is expected to show it in n*p of them; at n*p = .EFDR_MIN_EXPECTED the Poisson chance
# of seeing none is exp(-3) ~ 5%, so below that "we saw none" is not distinguishable from "we did
# not sample it". The test is therefore expected detections >= .EFDR_MIN_EXPECTED, where
#
#     expected = min(n_ctrl_cells, n_pert_cells) * max(det_ctrl, det_pert)
#
# p is read as the MAX of the two arms' detection rates: whichever arm carries the gene supplies
# the rate (that arm is well sampled, so the rate is well estimated), and the thinner arm supplies
# only its cell COUNT, which is exact. Taking the max rather than the other arm's rate is what
# keeps a gain from a zero baseline testable -- the reference arm's rate would be 0 there.
#
# This REPLACES a flat per-experiment cell floor K = ceil(3 / p75(det_ctrl)). That floor asked only
# "how many cells are in this arm", so it discarded every gene in a small cell type regardless of
# expression, while waving through trace genes in large ones. It also inherited two underivable
# constants: the .75 percentile, and the notion of a single reference "well-detected" gene standing
# in for all of them. The per-call form keeps only .EFDR_MIN_EXPECTED, the one constant with a
# derivation. The K machinery is retained below solely so the old policy stays reproducible via
# `gates` for back-comparison.
.EFDR_MIN_EXPECTED  <- 3     # expected detections needed in the thinner arm (Poisson P(0)=0.05)
.EFDR_DET_PCTL      <- 0.75  # legacy: percentile locating the retired flat cell floor K
.EFDR_MIN_ARM_CELLS <- 25L   # legacy: fallback for K if the detection distribution is degenerate

# Every uncallability gate that exists, in the order it is applied. Every one pins empirical_p
# to 1. Which SUBSET runs is chosen by `gates` (see annotate_empirical_fdr_model); the default
# is .EFDR_DEFAULT_GATES below. Names are the stable public identifiers -- the config, the CLI
# and the diagnostics table all use these strings.
.EFDR_GATES <- c("not_replicated",         # detected in <min_replicate_embryos embryos in BOTH arms
                 "too_few_embryos",        # <2 perturbation embryos: no valid two-group contrast
                 "not_detected_anywhere",  # ~0 UMI in BOTH arms: any |z| is shrinkage noise
                 "gain_from_zero",         # UP call whose CONTROL arm has ~0 UMI
                 "gain_no_counts",         # UP call whose PERTURBATION (gaining) arm has ~0 UMI
                 "loss_from_trace",        # DOWN call detected in <MIN_SRC_DET of control cells
                 "too_few_cells",          # thinner arm not powered to see THIS gene (expected < 3)
                 "flat_cell_floor")        # legacy: thinner arm below the per-experiment floor K

# The gates applied unless a caller says otherwise -- three rules.
#
#   too_few_embryos  <2 perturbation embryos: no two-group contrast exists.
#   not_replicated   the gene was seen in <min_replicate_embryos distinct embryos in BOTH arms,
#                    i.e. it was never observed more than once anywhere.
#   too_few_cells    the thinner arm was not powered to see this gene (see .EFDR_MIN_EXPECTED).
#
# The previous policy's other five gates are off by default, each subsumed:
#   not_detected_anywhere  a gene in >=2 embryos of either arm necessarily has counts somewhere.
#   gain_from_zero         an absolute "control arm has 0 UMI" rule that could not tell one stray
#                          cell from sixteen, so it removed real ectopic activation on low-depth
#                          runs; the expected-detections test answers that question directly.
#   gain_no_counts         measured to remove ZERO otherwise-significant calls in all six test
#                          experiments once too_few_cells applies everywhere (it fires on up to
#                          10% of calls but never on one that would have survived).
#   flat_cell_floor        the retired per-experiment cell floor K; see .EFDR_MIN_EXPECTED above.
#   loss_from_trace        a flat 2% control-detection floor. Removed calls well replicated across
#                          20+ embryos purely for being lowly expressed, and 64-82% of what it
#                          uniquely removed is caught by expected detections on a derived quantity.
#
# All five remain selectable, so the previous policy is reproducible for back-comparison:
#   gates = c("too_few_embryos", "not_detected_anywhere", "gain_from_zero", "gain_no_counts",
#             "loss_from_trace", "flat_cell_floor")
.EFDR_DEFAULT_GATES <- c("too_few_embryos", "not_replicated", "too_few_cells")

.EFDR_MIN_REP_EMB <- 2L

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
    # Thin-extended SPARSE grid. The arm-matched null must SPAN arm thinness -- from as thin
    # as a valid two-embryo pseudobulk contrast allows (~3% of control embryos) up to a thick
    # anchor -- not sample densely in the narrow high band the old grid used. Cell-type rarity
    # turns the 3% embryo split into single-digit-cell arms (the artifact regime); the thick
    # point anchors the well-sampled / real-call region. Three points {thin, mid, thick}
    # reconstruct the tail surface as well as the old five (validated on GENE6) while dropping
    # the most expensive (thickest) splits. Going below ~3% embryos buys nothing: we condition
    # on per-cell-type arm CELL count, which already reaches 1-3 cells for rare cell types at
    # the 3% point, and <2 embryos cannot form a contrast (the >=2-pseudobulk gate).
    n_perturb_grid <- unique(pmax(3L, round(c(0.03, 0.08, 0.40) * n_ctrl)))
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
    # Per-cell-type embryo support: n_pert_pb still drives the <2-embryo gate.
    supp <- efdr_perturbation_support(sub, sample_group = sample_group, cell_group = cell_group,
                                      control_ids = "ctrl-inj", perturbation_col = perturbation_col)
    # Per-(gene, cell_type) pseudo-arm UMI counts -> K_eff (retained as a diagnostic).
    null_counts <- .efdr_arm_counts(sub, sample_group = sample_group, cell_group = cell_group,
                                    control_ids = "ctrl-inj", perturbation_col = perturbation_col)
    # Per-(cell_type) pseudo-arm CELL counts -- the arm-size axis the tail is conditioned on.
    # The pseudo-KO ("NULLPERT") arm is the thinner arm here (a minority split of controls),
    # so its cell count spans the thin regime the artifacts live in.
    .cg <- as.character(SummarizedExperiment::colData(sub)[[cell_group]])
    .pc <- as.character(SummarizedExperiment::colData(sub)[[perturbation_col]])
    arm_cells <- data.frame(
      cell_group   = names(tapply(.pc == "NULLPERT", .cg, sum)),
      n_pert_cells = as.integer(tapply(.pc == "NULLPERT", .cg, sum)),
      n_ctrl_cells = as.integer(tapply(.pc == "ctrl-inj", .cg, sum)),
      stringsAsFactors = FALSE, row.names = NULL)
    .read_null_dir(wd, counts = null_counts) %>%
      dplyr::left_join(det, by = c("gene_short_name", "cell_group")) %>%
      dplyr::left_join(supp[, c("cell_group", "n_pert_pb")], by = "cell_group") %>%
      dplyr::left_join(arm_cells, by = "cell_group") %>%
      dplyr::mutate(log_ratio = log2(N / (n_ctrl - N)))
  })
}

# Effective count backing a within-state contrast: the harmonic combination of the per-arm
# total UMI counts, K_eff = 1/(1/K_ctrl + 1/K_pert), the reciprocal of the Poisson log-ratio
# variance (a near-empty arm drives K_eff -> ~0; a well-counted contrast gives K_eff in the
# hundreds). RETAINED AS A DIAGNOSTIC ONLY: the count-calibrated t that once scored z against
# a K_eff-calibrated t-distribution has been REMOVED in favour of the arm-matched tail (which
# conditions on the thin-arm CELL count) plus the directional detection/absence gates. K_eff
# is still computed and stored on the null/DEG so a demotion can be inspected against the arm
# counts, but it no longer drives any p-value. `kappa` guards the K = 0 singularity.
.efdr_keff <- function(K_ctrl, K_pert, kappa = 1) {
  1 / (1 / (K_ctrl + kappa) + 1 / (K_pert + kappa))
}

# Per-(gene, cell_type, arm) total UMI counts from a cds, as a long table
# (id, gene_short_name, cell_group, K_ctrl, K_pert). One grouped sum of the count matrix over
# the cell_type x arm indicator -- the cheap one-off during a load we already do. Feeds the
# count-calibrated t (K_eff = .efdr_keff(K_ctrl, K_pert)).
#
# KEYED ON `id` (= cds rownames = the DEG's unique feature id), NOT gene_short_name: multi-copy
# families (histones etc.) share a gene_short_name across many Ensembl ids, so a gene_short_name
# key is one-to-many and a downstream join on it silently duplicates DEG rows. `id` is unique.
.efdr_arm_counts <- function(cds, sample_group = "embryo_ID", cell_group = "cell_type",
                             control_ids = c("ctrl-inj"), perturbation_col = "perturbation") {
  cd <- SummarizedExperiment::colData(cds)
  ct <- as.character(cd[[cell_group]])
  arm <- ifelse(as.character(cd[[perturbation_col]]) %in% control_ids, "ctrl", "pert")
  grp <- factor(paste(ct, arm, sep = "\r"))
  ind <- Matrix::sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                              dims = c(length(grp), nlevels(grp)))
  # Materialize once: a lazy BPCells matrix's `>` comparison doesn't reliably survive a
  # subsequent matmul + as.matrix() (silently comes out NA), the same class of problem
  # efdr_detection_rate()'s .rate() helper already coerces around above. Doing it once here
  # also avoids materializing exprs(cds) twice (once for S, once for D).
  # drop0: a materialized BPCells matrix can carry EXPLICIT ZEROS in its sparsity pattern
  # (4.65% of stored entries in GENE8's contrast_cds; 0% in SGSeq2's). The binarisation below
  # sets every STORED entry to 1, so without dropping them a stored zero counts as a cell
  # expressing the gene -- inflating n_expr_* and hence det_ctrl, which drives the
  # low-source-detection gate, and n_emb_*. Symptom in a decorated table: K_ctrl < n_expr_ctrl,
  # arithmetically impossible, in 19.6% of rows of GENE8 gata5 as shipped.
  M <- Matrix::drop0(methods::as(monocle3::exprs(cds), "dgCMatrix"))
  S <- as.matrix(M %*% ind)                              # genes x (cell_type|arm) raw UMI totals
  # cells expressing the gene per (gene, cell_type|arm) -- binarized matmul, same grouping. Feeds
  # the per-arm DETECTION (expressing cells / arm cells) used by the down-call source-detection gate.
  Db <- M; Db@x[] <- 1                                   # binarise (stored entries are always > 0 for counts)
  D <- as.matrix(Db %*% ind)                             # genes x (cell_type|arm) cells expressing
  # Distinct EMBRYOS detecting the gene per (gene, cell_type|arm) -- the replicate-level count the
  # `not_replicated` gate needs. Accumulated one embryo at a time and kept sparse: the direct
  # route (grouping by cell_type|arm|embryo in one matmul) is genes x ~C*2*E, which densifies to
  # GBs on a full run. Per embryo the product is genes x (cell_type|arm), binarised to "this embryo
  # saw it" and summed, so the accumulator never leaves the size of D.
  smp <- as.character(cd[[sample_group]])
  Nemb <- matrix(0L, nrow(M), nlevels(grp))              # genes x (cell_type|arm), dense int
  for (s in unique(smp)) {
    k <- which(smp == s)
    Gs <- Matrix::sparseMatrix(i = seq_along(k), j = as.integer(grp[k]), x = 1,
                               dims = c(length(k), nlevels(grp)))
    Xs <- methods::as(Db[, k, drop = FALSE] %*% Gs, "TsparseMatrix")
    # +1 only on the (gene, arm) entries this embryo touched: O(nnz) per embryo rather than a
    # whole-matrix add, which is what makes ~500 embryos affordable.
    ij <- cbind(Xs@i + 1L, Xs@j + 1L)
    Nemb[ij] <- Nemb[ij] + 1L
  }
  ids <- rownames(cds)
  gsym <- SummarizedExperiment::rowData(cds)$gene_short_name
  if (is.null(gsym)) gsym <- ids
  key <- strsplit(levels(grp), "\r", fixed = TRUE)
  ct_lev <- vapply(key, `[`, "", 1); arm_lev <- vapply(key, `[`, "", 2)
  # Build (id, cell_type) x {K_ctrl, K_pert, n_expr_ctrl, n_expr_pert} per cell type (no long/pivot).
  res <- lapply(unique(ct_lev), function(ctype) {
    jc <- which(ct_lev == ctype & arm_lev == "ctrl")
    jp <- which(ct_lev == ctype & arm_lev == "pert")
    data.frame(id = ids, gene_short_name = gsym, cell_group = ctype,
               K_ctrl = if (length(jc)) S[, jc] else 0,
               K_pert = if (length(jp)) S[, jp] else 0,
               n_expr_ctrl = if (length(jc)) D[, jc] else 0,
               n_expr_pert = if (length(jp)) D[, jp] else 0,
               n_emb_ctrl = if (length(jc)) Nemb[, jc] else 0,
               n_emb_pert = if (length(jp)) Nemb[, jp] else 0, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

# Join per-(gene, cell) arm counts onto a table `x` on the unique feature id when both carry
# it, else fall back to gene_short_name (deduplicating the counts by summing so the fallback
# can never be one-to-many). Shared by the null build and the real-DEG annotation so K_eff is
# computed the same way on both sides.
.efdr_join_counts <- function(x, counts) {
  has_det <- all(c("n_expr_ctrl", "n_expr_pert") %in% names(counts))
  has_emb <- all(c("n_emb_ctrl", "n_emb_pert") %in% names(counts))
  if ("id" %in% names(x) && "id" %in% names(counts)) {
    cols <- c("id", "cell_group", "K_ctrl", "K_pert",
              if (has_det) c("n_expr_ctrl", "n_expr_pert"),
              if (has_emb) c("n_emb_ctrl", "n_emb_pert"))
    dplyr::left_join(x, counts[, cols], by = c("id", "cell_group"))
  } else {
    # gene_short_name fallback: UMI and cell counts SUM across a symbol's ids, but embryo counts
    # cannot -- the same embryo can detect two ids of one symbol, so summing double-counts the
    # replicate. max() is the honest bound available without the per-embryo detail (it is the
    # count for the single best-detected id, so it never overstates the number of distinct embryos).
    cnt <- counts %>%
      dplyr::group_by(.data$gene_short_name, .data$cell_group) %>%
      dplyr::summarise(K_ctrl = sum(.data$K_ctrl), K_pert = sum(.data$K_pert),
                       n_expr_ctrl = if (has_det) sum(.data$n_expr_ctrl) else NA_real_,
                       n_expr_pert = if (has_det) sum(.data$n_expr_pert) else NA_real_,
                       n_emb_ctrl  = if (has_emb) max(.data$n_emb_ctrl)  else NA_real_,
                       n_emb_pert  = if (has_emb) max(.data$n_emb_pert)  else NA_real_,
                       .groups = "drop")
    dplyr::left_join(x, cnt, by = c("gene_short_name", "cell_group"))
  }
}

.read_null_dir <- function(wd, counts = NULL) {
  files <- list.files(wd, "_within_node_degs.csv$", full.names = TRUE)
  if (!length(files)) return(tibble::tibble())
  d <- purrr::map_dfr(files, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE))) %>%
    dplyr::filter(is.finite(.data$perturb_to_ctrl_shrunken_lfc),
                  .data$perturb_to_ctrl_shrunken_lfc_se > 0,
                  is.finite(.data$log_mean_expression))
  d$z <- d$perturb_to_ctrl_shrunken_lfc / d$perturb_to_ctrl_shrunken_lfc_se   # raw magnitude z
  # Per-(gene, cell_type) effective count, so the count-calibrated t is calibrated against the
  # same quantity the real DEG is scored on.
  if (!is.null(counts)) {
    d <- .efdr_join_counts(d, counts)                   # id-keyed (falls back to summed symbol)
    d$K_eff <- .efdr_keff(d$K_ctrl, d$K_pert)
  } else {
    d$K_eff <- NA_real_
  }
  tibble::tibble(gene_short_name = d$gene_short_name, cell_group = d$cell_group,
                 log_mean_expression = d$log_mean_expression, z = d$z, K_eff = d$K_eff) %>%
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
#' Models the null upper tail of |z| as a LOCAL-EMPIRICAL surface over **expression
#' and the thin-arm cell count**, per direction. Each grid cell's tail quantiles are
#' estimated from the null draws in a local window around it (widened adaptively until
#' `min_n` draws are in; still-empty cells filled from the nearest evaluated cell;
#' quantiles made monotone in tau by cummax). The two axes:
#'   - **expression**: the dispersion mis-pricing that inflates |z| worsens as
#'     expression -> 0, so the tail stays heavy in the sparse extreme-low region
#'     rather than drooping back toward zero.
#'   - **thin-arm cell count**: a thin arm hits exact/near-zero counts which, times
#'     the trend-dispersion-shrunk SE, manufactures large |z| under no real effect.
#'     This is the axis that separates a real depletion (huge |z|, deep arm) from a
#'     thin-arm sampling artifact -- the previous %embryos/log-ratio surface could
#'     not, because it never conditioned on the actual arm size.
#' A local-empirical estimate is used rather than a global shape-constrained smooth
#' (an earlier scam `bs="mpd"` surface), which under-fit the low-expression/thin-arm
#' corner by borrowing from lighter, better-sampled neighbours.
#' Draws are POOLED across the sampling grid: the grid exists only to span the
#' arm-size axis, which is then read directly off `n_pert_cells`. Direction is kept
#' separate. The thin arm in a control-split null is always the pseudo-KO
#' ("NULLPERT") minority arm, so its cell count (`n_pert_cells`) is the covariate;
#' `dn` = the gene is lower in that thin arm (z<0), `up` = higher.
#'
#' @param null Output of [build_empirical_null()] (gene_short_name, cell_group,
#'   log_mean_expression, z, n_pert_cells).
#' @param detection Unused (kept for back-compat with callers).
#' @param taus Quantile levels (tail-concentrated) the surface is evaluated at.
#' @param egrid,lcgrid Expression and log10(thin-arm cells + 1) grids evaluated on.
#' @param we,wc Half-widths of the local window (expression, log10-cells) the tail
#'   quantiles are estimated over; widened adaptively until `min_n` draws are in.
#' @param min_n Minimum null draws required in a window before its quantiles are used.
#' @return An `efdr_model`: per-direction 2-D tail surfaces on (egrid x lcgrid).
#' @export
train_efdr_model <- function(null, detection = NULL, taus = .EFDR_TAUS,
                             egrid = .EFDR_EGRID, lcgrid = .EFDR_LCGRID,
                             we = 0.75, wc = 0.4, min_n = 40L, max_we = 2.5) {
  if (!"n_pert_cells" %in% names(null))
    stop("efdr: null lacks `n_pert_cells` (arm cell count). Rebuild with the current ",
         "build_empirical_null(), which records per-cell-type pseudo-arm cell counts.",
         call. = FALSE)
  null <- null %>%
    dplyr::filter(is.finite(.data$log_mean_expression), is.finite(.data$z),
                  is.finite(.data$n_pert_cells), .data$n_pert_cells > 0) %>%
    dplyr::mutate(dir = dplyr::if_else(.data$z < 0, "dn", "up"), a = abs(.data$z),
                  lc = log10(.data$n_pert_cells + 1))          # thin-arm size axis
  GR <- expand.grid(eb = egrid, pe = lcgrid)          # eb varies fastest; `pe` slot holds log-cells
  # aspect ratio to make grid-space "distance" comparable across the two axes when NN-filling.
  asp <- diff(range(egrid)) / max(diff(range(lcgrid)), 1e-6)
  surf <- list()
  for (d in c("dn", "up")) {
    sub <- null %>% dplyr::filter(.data$dir == d)
    if (nrow(sub) < 8L * min_n) { surf[[d]] <- NULL; next }   # too sparse -> caller defers to ashr
    e <- sub$log_mean_expression; l <- sub$lc; a <- sub$a
    ord <- order(e); e <- e[ord]; l <- l[ord]; a <- a[ord]    # sort by expr for fast window slicing
    # LOCAL EMPIRICAL tail on the grid (a kNN baked to a lookup). The scam monotone-smooth
    # under-fit the (low-expression, thin-arm) corner where the artifacts live -- borrowing
    # from lighter neighbours -- so estimate each grid cell from its own local window instead.
    Q <- matrix(NA_real_, nrow(GR), length(taus))
    for (r in seq_len(nrow(GR))) {
      eb <- GR$eb[r]; pe <- GR$pe[r]; ww <- we; wl <- wc; keep <- integer(0)
      repeat {
        lo <- findInterval(eb - ww, e) + 1L; hi <- findInterval(eb + ww, e)
        if (hi >= lo) { j <- lo:hi; keep <- j[abs(l[j] - pe) <= wl] }
        if (length(keep) >= min_n || ww >= max_we) break
        ww <- ww * 1.5; wl <- wl * 1.5
      }
      if (length(keep) >= 30L) Q[r, ] <- stats::quantile(a[keep], taus, names = FALSE)
    }
    # fill any still-empty grid cell from the nearest evaluated cell (grid-space NN), so the
    # extreme corners are covered rather than left to caller-side extrapolation.
    miss <- which(is.na(Q[, 1])); have <- which(!is.na(Q[, 1]))
    if (length(miss) && length(have)) for (m in miss) {
      dd <- (GR$eb[have] - GR$eb[m])^2 + ((GR$pe[have] - GR$pe[m]) * asp)^2
      Q[m, ] <- Q[have[which.min(dd)], ]
    }
    Q <- t(apply(Q, 1, cummax))                        # monotone in tau within each grid cell
    surf[[d]] <- Q
  }
  structure(list(surf = surf, taus = taus, egrid = egrid, pegrid = lcgrid),
            class = "efdr_model")
}

# right-tail probability of |z| = `a` at (expression `e`, arm-size `pe`) against a
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

#' Annotate a DEG table with per-gene empirical p and BH FDR.
#'
#' @param model An `efdr_model` from [train_efdr_model()].
#' @param deg_tbl DEG table with `gene_short_name`, `cell_group`,
#'   `log_mean_expression`, and either `z` or `perturb_to_ctrl_shrunken_lfc` +
#'   `perturb_to_ctrl_shrunken_lfc_se`.
#' @param arm_cells Per-cell-type arm CELL counts (data.frame with `cell_group`,
#'   `n_ctrl_cells`, `n_pert_cells`). Supplies the arm-size covariate the tail is
#'   conditioned on. Genes in cell types with no match fall back to ashr only.
#' @param counts Per-(gene, cell_type) arm UMI totals (`id`/`gene_short_name`,
#'   `cell_group`, `K_ctrl`, `K_pert`) for the gene-absent and control-absent gates.
#' @param support Per-cell-type embryo counts (`cell_group`, `n_pert_pb`) for the
#'   <2-pseudobulk gate.
#' @param detection Deprecated / unused (the %embryos covariate is retired).
#' @param log_ratio Deprecated / ignored (the tail is conditioned on the actual
#'   per-cell-type arm cell count, not a single global sampling ratio).
#' @param gates Character vector of the uncallability gates to APPLY, from
#'   `.EFDR_GATES`. Defaults to `.EFDR_DEFAULT_GATES` = `too_few_embryos`,
#'   `not_replicated`, `too_few_cells`. The previous policy's other five gates
#'   (`not_detected_anywhere`, `gain_from_zero`, `gain_no_counts`,
#'   `loss_from_trace`, and `flat_cell_floor`, the retired per-experiment cell
#'   floor K) are off by default because each is subsumed; pass them explicitly
#'   to reproduce that policy for back-comparison. Gates not listed are still
#'   computed and counted in the `efdr_stats` attribute -- they just stop pinning
#'   `empirical_p` to 1, so those calls fall back to the two statistical floors
#'   (`max(p_ashr, p_tail)`). Unknown names error rather than silently no-op.
#' @param min_replicate_embryos Distinct embryos that must detect the gene in at
#'   least one arm for the `not_replicated` gate to pass it (default
#'   `.EFDR_MIN_REP_EMB` = 2). Needs `counts` to carry `n_emb_ctrl` / `n_emb_pert`
#'   (produced by `.efdr_arm_counts`); without them nothing is gated on this axis.
#' @param group_col Column to BH-adjust within (default `cell_group`).
#' @return `deg_tbl` with added `empirical_p` and `empirical_fdr`.
#' @details The tail null is conditioned on **expression x thin-arm cell count**
#'   (not %embryos / global ratio), so a real depletion (huge |z|, deep arm) is
#'   kept while a thin-arm sampling artifact of the same |z| is demoted. Each call
#'   is matched on `min(n_ctrl_cells, n_pert_cells)` (the thinner arm) and the sign
#'   of the effect **in that thin arm**, which generalizes to perturbation-heavy
#'   designs (where the control arm is the thin one). `empirical_p = max(p_ashr,
#'   p_tail)`: the empirical step may only *demote* a call ashr made.
#' @export
annotate_empirical_fdr_model <- function(model, deg_tbl, log_ratio = NULL, detection = NULL,
                                         support = NULL, counts = NULL, arm_cells = NULL,
                                         min_pseudobulks = .EFDR_MIN_PB,
                                         gates = .EFDR_DEFAULT_GATES,
                                         min_replicate_embryos = .EFDR_MIN_REP_EMB,
                                         group_col = "cell_group") {
  gates <- setdiff(as.character(gates), c("", NA_character_))
  if (length(bad <- setdiff(gates, .EFDR_GATES)))
    stop("unknown gate(s): ", paste(bad, collapse = ", "),
         ". Valid: ", paste(.EFDR_GATES, collapse = ", "))
  # Gates are always COMPUTED (the diagnostics report membership either way); `on()` decides
  # whether a gate is APPLIED to empirical_p.
  on <- function(g) g %in% gates
  message(sprintf("[efdr] gates applied: %s (min_replicate_embryos=%d)",
                  paste(gates, collapse = ", "), as.integer(min_replicate_embryos)))
  # (Re)join the per-cell-type inputs, dropping any stale copies first so re-decoration
  # of an already-decorated table doesn't collide into `.x`/`.y` and silently NA a covariate.
  if (!is.null(support)) {                                # n_pert_pb -> <2-embryo gate
    deg_tbl <- deg_tbl[, setdiff(names(deg_tbl), "n_pert_pb"), drop = FALSE]
    deg_tbl <- dplyr::left_join(deg_tbl, support[, c("cell_group", "n_pert_pb")], by = "cell_group")
  }
  if (!is.null(counts)) {                                 # per-(gene, cell_type) UMI totals -> gates
    # Also drop n_expr_ctrl/n_expr_pert (the per-arm detection counts .efdr_join_counts adds
    # when `counts` carries them): re-decorating an already-arm-matched-decorated table
    # otherwise collides into n_expr_ctrl.x/.y (dplyr's duplicate-column suffixing), and
    # det_ctrl's later reference to the plain column name silently resolves to nothing.
    deg_tbl <- deg_tbl[, setdiff(names(deg_tbl), c("K_ctrl", "K_pert", "n_expr_ctrl", "n_expr_pert",
                                                   "n_emb_ctrl", "n_emb_pert")), drop = FALSE]
    deg_tbl <- .efdr_join_counts(deg_tbl, counts)         # id-keyed (falls back to summed symbol)
  }
  if (!is.null(arm_cells)) {                              # per-cell-type arm CELL counts -> tail axis
    deg_tbl <- deg_tbl[, setdiff(names(deg_tbl), c("n_ctrl_cells", "n_pert_cells")), drop = FALSE]
    deg_tbl <- dplyr::left_join(deg_tbl, arm_cells[, c("cell_group", "n_ctrl_cells", "n_pert_cells")],
                                by = "cell_group")
  }
  for (col in c("n_pert_pb", "K_ctrl", "K_pert", "n_ctrl_cells", "n_pert_cells",
                "n_expr_ctrl", "n_expr_pert", "n_emb_ctrl", "n_emb_pert"))
    if (!col %in% names(deg_tbl)) deg_tbl[[col]] <- NA_real_

  out <- deg_tbl %>%
    dplyr::mutate(
      z = if ("z" %in% names(.)) .data$z
          else .data$perturb_to_ctrl_shrunken_lfc / .data$perturb_to_ctrl_shrunken_lfc_se,
      .a = abs(.data$z),
      # thin arm = the smaller of the two arms; its cell count is the tail covariate.
      thin_arm_cells = pmin(.data$n_ctrl_cells, .data$n_pert_cells),
      .thin_is_ko = is.finite(.data$n_pert_cells) & is.finite(.data$n_ctrl_cells) &
                    .data$n_pert_cells <= .data$n_ctrl_cells,
      # direction RELATIVE TO THE THIN ARM (matches the null's convention, whose thin arm is the
      # pseudo-KO): dn = gene lower in the thin arm. For a control-heavy design thin=KO so this is
      # sign(z); for a perturbation-heavy design thin=control and it flips -- which is exactly the
      # generalization (the artifact becomes apparent GAINS when the control arm is the thin one).
      .dir = dplyr::if_else((.data$.thin_is_ko & .data$z < 0) |
                            (!.data$.thin_is_ko & .data$z > 0), "dn", "up"),
      # p_ashr: ordinary shrunken significance (demote-only baseline; empirical can only demote).
      p_ashr = if ("perturb_to_ctrl_p_value" %in% names(.)) .data$perturb_to_ctrl_p_value
               else 2 * (1 - stats::pnorm(.data$.a)))

  # empirical tail, conditioned on (expression x thin-arm cell count), per thin-arm direction.
  emp_tail <- rep(NA_real_, nrow(out))
  for (d in c("dn", "up")) {
    idx <- which(out$.dir == d & is.finite(out$thin_arm_cells) & out$thin_arm_cells > 0)
    if (!length(idx)) next
    e  <- out$log_mean_expression[idx]
    lc <- log10(out$thin_arm_cells[idx] + 1)             # arm-size axis (matches training)
    emp_tail[idx] <- .efdr_tail_2d(model$surf[[d]], model$taus, model$egrid, model$pegrid, e, lc, out$.a[idx])
  }
  out$p_tail <- pmin(pmax(emp_tail, 0), 1)               # thin-arm sampling-artifact floor
  # empirical_p = max of the two demote-only floors: ashr (ordinary significance) and tail
  # (thin-arm artifact). A missing tail drops out of the max (no fallback hole).
  out$empirical_p <- pmax(out$p_ashr, dplyr::coalesce(out$p_tail, 0))

  # REPLICABILITY (not_replicated): the gene must have been seen in at least
  # `min_replicate_embryos` distinct embryos of AT LEAST ONE arm. At ~500 UMI/cell most genes
  # are absent from any given cell, so catching the same gene in two separate specimens of one
  # arm and none of the other is not a coincidence the sampling produces easily -- which is
  # exactly what the arm-matched null prices, and why this is a better instrument than the
  # absolute zero-UMI rules it replaces. Direction-agnostic and arm-agnostic on purpose: it
  # asks "is this gene really measured in this cell type at all", not "which way did it move".
  # NA counts (no `counts` input) -> not gated, matching every other count-driven gate.
  out$max_emb_det <- pmax(out$n_emb_ctrl, out$n_emb_pert, na.rm = TRUE)
  # WHY a call was pinned, recorded alongside WHETHER it was. Until now a gated
  # call and a tested-and-clearly-null call were both empirical_p == 1 in the
  # published table, so "we could not call this gene here" and "we tested this
  # gene here and it was null" were literally the same number. Only aggregate
  # counts survived, in the efdr_stats attribute.
  #
  # This records the reason and changes nothing else: empirical_p, empirical_fdr
  # and every gate's behaviour are untouched, and a gate that is off records
  # nothing. NA means the row was not gated.
  #
  # Several gates can fire on one row, so reasons accumulate in evaluation order,
  # ";"-separated. Note that expected_detections already tells you about the
  # count-driven gates; it is the others -- not_replicated above all -- that are
  # otherwise invisible, because a row can pass every count test and still be
  # pinned.
  out$gate_reason <- NA_character_
  note_gate <- function(reason, hit, label) {
    hit <- hit & !is.na(hit)
    ifelse(hit, ifelse(is.na(reason), label, paste(reason, label, sep = ";")), reason)
  }

  not_replicated_gate <- is.finite(out$max_emb_det) & out$max_emb_det < min_replicate_embryos
  if (on("not_replicated")) out$empirical_p[not_replicated_gate] <- 1
  if (on("not_replicated")) out$gate_reason <- note_gate(out$gate_reason, not_replicated_gate, "not_replicated")
  # Gate 1 (<2 perturbation embryos): no valid two-group contrast; calls uncallable.
  gate <- is.finite(out$n_pert_pb) & out$n_pert_pb < min_pseudobulks
  if (on("too_few_embryos")) out$empirical_p[gate] <- 1
  if (on("too_few_embryos")) out$gate_reason <- note_gate(out$gate_reason, gate, "too_few_embryos")
  # Gate 2 (gene absent from the cell type): ~zero UMIs in BOTH arms -> not expressed here at all;
  # any |z| is pure shrinkage noise no floor can price. Gate only when the LARGER arm is empty
  # (max(K)), so a deep control arm with an empty perturbation arm = genuine COMPLETE DEPLETION is
  # KEPT. Not an expression floor: "the gene has no counts in either arm of this cell type".
  absent_gate <- is.finite(out$K_ctrl) & is.finite(out$K_pert) &
    pmax(out$K_ctrl, out$K_pert) < .EFDR_MIN_ARM_UMI
  if (on("not_detected_anywhere")) out$empirical_p[absent_gate] <- 1
  if (on("not_detected_anywhere")) out$gate_reason <- note_gate(out$gate_reason, absent_gate, "not_detected_anywhere")
  # Gate 3 (control-absent UP-call): the mirror of the gene-absent gate for gains. An UP call
  # (higher in perturbation) whose CONTROL arm carries ~zero UMIs is an "increase" over a baseline
  # we never detected -- uncallable, and the control-split null cannot price it (control-vs-control
  # is 0 in both arms). A DOWN call with an empty perturbation arm is the opposite case (a real
  # depletion) and is NOT gated. Directional on the control (reference) arm.
  up_absent_gate <- (out$z > 0) & is.finite(out$K_ctrl) & (out$K_ctrl < .EFDR_MIN_ARM_UMI)
  if (on("gain_from_zero")) out$empirical_p[up_absent_gate] <- 1
  if (on("gain_from_zero")) out$gate_reason <- note_gate(out$gate_reason, up_absent_gate, "gain_from_zero")
  # Gate 3c (empty-gaining-arm UP-call): an UP call (higher in perturbation) whose PERTURBATION
  # (gaining) arm carries ~zero UMIs cannot be a real gain -- the gene is absent from the very arm
  # it is called up in, so the shrunken-LFC sign is noise (a near-zero control arm plus size-factor
  # shrinkage can emit z>0 with K_pert=0). This is the true up-mirror of the low-source-detection
  # LOSS gate on the GAINING side, and the counterpart of the gene-absent gate: gate 2 keeps a gene
  # present in EITHER arm (max(K)), which lets a K_pert=0 up-call survive on its control counts, so
  # this gate closes that hole on the up side. Deliberately a ZERO-count floor, not the 2% detection
  # floor used for losses: real GAINS can be sparse (phox2a is real at 0.7% perturbation detection,
  # K_pert>0), so only a truly empty gaining arm is gated. Directional on the perturbation arm.
  up_empty_gain_gate <- (out$z > 0) & is.finite(out$K_pert) & (out$K_pert < .EFDR_MIN_ARM_UMI)
  if (on("gain_no_counts")) out$empirical_p[up_empty_gain_gate] <- 1
  if (on("gain_no_counts")) out$gate_reason <- note_gate(out$gate_reason, up_empty_gain_gate, "gain_no_counts")
  # Gate 3b (low source-detection LOSS): a DOWN call whose CONTROL (source) arm detects the gene in
  # fewer than MIN_SRC_DET of its cells has no robustly-expressed baseline to lose -- its large |z|
  # is trend-dispersion shrinkage noise the control-split null cannot reach at moderate arm sizes
  # (the residual the count-t used to catch). Directional and DOWN-ONLY: real gains can be sparse
  # (a real up-call can sit below this in the perturbation arm), so no symmetric up gate. Data-derived
  # floor sitting in the ~10x gap between real depletions (>=12% control detection) and the confirmed
  # low-expression artifacts (<=1.1%). NA detection -> not gated (no data to judge).
  out$det_ctrl <- out$n_expr_ctrl / out$n_ctrl_cells
  low_src_det_gate <- (out$z < 0) & is.finite(out$det_ctrl) & (out$det_ctrl < .EFDR_MIN_SRC_DET)
  if (on("loss_from_trace")) out$empirical_p[low_src_det_gate] <- 1
  if (on("loss_from_trace")) out$gate_reason <- note_gate(out$gate_reason, low_src_det_gate, "loss_from_trace")
  # too_few_cells (BOTH directions): was the thinner arm powered to see THIS gene? Expected
  # detections = (cells in the thinner arm) x (the better-estimated of the two arms' detection
  # rates); below .EFDR_MIN_EXPECTED, "we saw none / fewer" cannot be told from "we did not sample
  # it" (derivation at the constant). A complete absence in the thin arm is a real depletion only
  # if that arm could have detected the gene -- brd1b in a cell type depleted to 8 perturbation
  # cells, none expressing, is a too-few-cells call, not a loss. Being per call, it keeps a
  # well-detected gene callable in a small cell type and gates a trace gene even in a large one;
  # the flat floor it replaces could do neither. Symmetric on purpose: it also covers thin-arm
  # GAINS, which the control-split null structurally cannot generate. NA counts -> not gated.
  out$det_pert  <- out$n_expr_pert / out$n_pert_cells
  out$max_det   <- pmax(out$det_ctrl, out$det_pert, na.rm = TRUE)
  out$expected_detections <- out$thin_arm_cells * out$max_det
  underpowered_gate <- is.finite(out$expected_detections) &
    out$expected_detections < .EFDR_MIN_EXPECTED
  if (on("too_few_cells")) out$empirical_p[underpowered_gate] <- 1
  if (on("too_few_cells")) out$gate_reason <- note_gate(out$gate_reason, underpowered_gate, "too_few_cells")

  # flat_cell_floor: the RETIRED per-experiment cell floor, off by default. Kept so the previous
  # policy is reproducible via `gates` for back-comparison. K is located from the p75 of control
  # detection, which is why it varied 24-58 across experiments and tracked cell-type size
  # composition as much as capture depth.
  min_arm_cells <- NA_integer_; .p_det <- NA_real_
  if (on("flat_cell_floor")) {
    .p_det <- suppressWarnings(stats::quantile(
      out$det_ctrl[is.finite(out$det_ctrl) & out$det_ctrl > 0], .EFDR_DET_PCTL, names = FALSE))
    min_arm_cells <- if (is.finite(.p_det) && .p_det > 0)
      as.integer(ceiling(.EFDR_MIN_EXPECTED / .p_det)) else .EFDR_MIN_ARM_CELLS
    message(sprintf("[efdr] legacy flat cell floor K=%d (p%.0f control detection = %.3f)",
                    min_arm_cells, 100 * .EFDR_DET_PCTL, .p_det))
    flat_floor_gate <- is.finite(out$thin_arm_cells) & out$thin_arm_cells < min_arm_cells
    out$empirical_p[flat_floor_gate] <- 1
    out$gate_reason <- note_gate(out$gate_reason, flat_floor_gate, "flat_cell_floor")
  }

  # empirical_fdr: BH-adjust EACH floor within the cell group, then take the max q -- NOT
  # BH(empirical_p) (the per-gene max compresses the p-distribution so BH crushes real hits).
  # Adjusting each floor in its own distribution lets the tiny ashr q survive while the gene must
  # still clear both q's. A missing tail -> 0 (non-constraining), matching the empirical_p max().
  out <- out %>%
    dplyr::group_by(dplyr::across(all_of(group_col))) %>%
    dplyr::mutate(
      .q_ashr = stats::p.adjust(pmax(.data$p_ashr, 0), "BH"),
      .q_tail = stats::p.adjust(dplyr::coalesce(.data$p_tail, 0), "BH"),
      empirical_fdr = pmax(.data$.q_ashr, .data$.q_tail)) %>%
    dplyr::ungroup()
  out$empirical_fdr[out$empirical_p >= 1] <- 1           # gated (uncallable) genes: fdr = 1 too
  out <- out %>%
    dplyr::select(-".dir", -".a", -".thin_is_ko", -".q_ashr", -".q_tail")  # keep p_ashr/p_tail/thin_arm_cells/K_*
  # Per-experiment efdr diagnostics (attached as an attribute so callers -- e.g. the mcclintock
  # decorator stage -- can persist them; the located min-arm K in particular is otherwise only in
  # the log). Counts are gate MEMBERSHIP (a call may satisfy more than one gate, so they need not sum).
  attr(out, "efdr_stats") <- data.frame(
    min_arm_cells_K       = as.integer(min_arm_cells),   # NA unless flat_cell_floor is on
    p_ctrl_detection      = round(as.numeric(.p_det), 4),
    det_pctl              = .EFDR_DET_PCTL,
    gates_applied         = paste(gates, collapse = ","),
    min_replicate_embryos = as.integer(min_replicate_embryos),
    n_calls               = nrow(out),
    n_ashr_sig            = sum(out$p_ashr < 0.05, na.rm = TRUE),
    n_retained            = sum(out$empirical_p < 0.05, na.rm = TRUE),
    n_retained_down       = sum(out$empirical_p < 0.05 & out$z < 0, na.rm = TRUE),
    n_retained_up         = sum(out$empirical_p < 0.05 & out$z > 0, na.rm = TRUE),
    n_gate_lt2_pseudobulk = sum(gate, na.rm = TRUE),
    n_gate_gene_absent    = sum(absent_gate, na.rm = TRUE),
    n_gate_up_ctrl_absent = sum(up_absent_gate, na.rm = TRUE),
    n_gate_up_empty_gain  = sum(up_empty_gain_gate, na.rm = TRUE),
    n_gate_low_ctrl_det   = sum(low_src_det_gate, na.rm = TRUE),
    n_gate_too_few_cells  = sum(underpowered_gate, na.rm = TRUE),
    n_gate_not_replicated = sum(not_replicated_gate, na.rm = TRUE),
    stringsAsFactors = FALSE)
  out
}
