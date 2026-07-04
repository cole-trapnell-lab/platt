test_that("annotate_empirical_fdr returns calibrated, direction-specific FDR", {
  # Synthetic: two abundance-matched unaffected cell types with a 10% low-expr
  # down-call baseline (the artifact), and one affected cell type at 40% down
  # (30% real signal on top of the 10% artifact baseline). All at low expression.
  mk <- function(cg, n_down, n_total = 200L, M = -4) {
    n_null <- n_total - n_down
    data.frame(
      cell_group = cg,
      id = paste0(cg, "_", seq_len(n_total)),
      log_mean_expression = M,
      perturb_to_ctrl_shrunken_lfc = c(rep(-1, n_down), rep(0.5, n_null)),
      perturb_to_ctrl_p_value = c(rep(0.001, n_down), rep(0.9, n_null)),
      mean_log_sf = log(5),
      stringsAsFactors = FALSE
    )
  }
  deg <- rbind(mk("U1", 20L), mk("U2", 20L), mk("AFF", 80L))
  abund <- c(U1 = 1000, U2 = 1000, AFF = 1000)

  ann <- annotate_empirical_fdr(deg, abundances = abund,
                                unaffected_cell_types = c("U1", "U2"))

  expect_true(all(c("empirical_null_rate", "observed_rate", "empirical_fdr") %in%
                    colnames(ann)))
  expect_equal(nrow(ann), nrow(deg))
  expect_true(all(ann$empirical_fdr >= 0 & ann$empirical_fdr <= 1, na.rm = TRUE))

  # Unaffected cell type: observed down-rate == null down-rate -> FDR ~ 1
  fdr_u <- unique(ann$empirical_fdr[ann$cell_group == "U1" &
                                      ann$perturb_to_ctrl_shrunken_lfc < 0])
  expect_true(all(fdr_u > 0.9))

  # Affected cell type: observed 40% vs null 10% -> FDR ~ 0.25
  fdr_a <- unique(ann$empirical_fdr[ann$cell_group == "AFF" &
                                      ann$perturb_to_ctrl_shrunken_lfc < 0])
  expect_equal(fdr_a, 0.25, tolerance = 0.02)
})

test_that("annotate_empirical_fdr errors on missing required columns", {
  expect_error(annotate_empirical_fdr(data.frame(cell_group = "a")),
               "missing required columns")
})

test_that("annotate_empirical_fdr falls back to lower-envelope self-null", {
  mk <- function(cg, n_down) data.frame(
    cell_group = cg, id = paste0(cg, "_", 1:200), log_mean_expression = -4,
    perturb_to_ctrl_shrunken_lfc = c(rep(-1, n_down), rep(0.5, 200 - n_down)),
    perturb_to_ctrl_p_value = c(rep(0.001, n_down), rep(0.9, 200 - n_down)),
    mean_log_sf = log(5), stringsAsFactors = FALSE)
  deg <- rbind(mk("a", 10L), mk("b", 12L), mk("c", 80L))
  ann <- annotate_empirical_fdr(deg, abundances = c(a = 1000, b = 1000, c = 1000))
  expect_true("empirical_fdr" %in% colnames(ann))
  # the high-rate cell type (c) should be flagged less than the low-rate ones
  fdr_c <- unique(ann$empirical_fdr[ann$cell_group == "c" &
                                      ann$perturb_to_ctrl_shrunken_lfc < 0])
  expect_true(all(fdr_c < 1))
})
