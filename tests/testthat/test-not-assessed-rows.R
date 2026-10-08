# Impact-table rows for cell types the perturbation could not assess are kept and marked, not dropped.

.impact_results <- function() {
  deg <- tibble::tibble(id = "g1", gene_short_name = "a", perturb_to_ctrl_p_value = 0.01)
  pw <- data.table::data.table(pathway = "p", pval = 0.01, padj = 0.02, overlap = 3L, size = 10L,
                               overlapGenes = list(c("a", "b", "c")), ontology = "BP")
  tibble::tibble(
    cell_type = c("ct_present", "ct_below", "ct_ablated"),
    summary = c("full summary", "full summary", "full summary"),
    abundance_code = c("A2 Depletion", NA, "A3 Ablation"),
    present_above_thresh = c(TRUE, FALSE, FALSE),
    has_deg_rows = c(TRUE, TRUE, FALSE),
    degs = list(deg, deg, deg[0, ]),
    pathways = list(pw, pw, pw[0, ]),
    goi_line = c("goi", "goi", "goi"),
    dact_grounding = c("lit", "lit", "lit"),
    has_regulator_goi = c(TRUE, TRUE, FALSE),
    llm_summary = c("LLM text", "No cell-autonomous signal for 'ct_below'", "No cell-autonomous signal for 'ct_ablated'"),
    llm_disrupted_pathways = list(NULL, list(list(name = "x")), NULL),
    llm_other_dysregulated_genes = list("geneA", "geneB", NULL)
  )
}

testthat::test_that("not-assessed cell types keep a row, marked with a reason", {
  out <- mark_not_assessed_cell_types(.impact_results())
  expect_equal(out$cell_type, c("ct_present", "ct_below", "ct_ablated"))
  expect_equal(out$not_assessed_reason, c(NA, "below_presence_threshold", "no_deg_rows"))
  expect_false(any(c("present_above_thresh", "has_deg_rows") %in% names(out)))
  # codes are kept: the ablated cell type's A3 is no longer lost
  expect_equal(out$abundance_code, c("A2 Depletion", NA, "A3 Ablation"))
})

testthat::test_that("not-assessed rows carry no interpretation; assessed rows are untouched", {
  inp <- .impact_results()
  out <- mark_not_assessed_cell_types(inp)
  expect_equal(out[1, setdiff(names(inp), c("present_above_thresh", "has_deg_rows"))],
               inp[1, setdiff(names(inp), c("present_above_thresh", "has_deg_rows"))])
  expect_true(all(is.na(out$llm_summary[2:3])))
  expect_true(all(vapply(out$llm_disrupted_pathways[2:3], is.null, logical(1))))
  expect_true(all(vapply(out$llm_other_dysregulated_genes[2:3], is.null, logical(1))))
  expect_equal(vapply(out$degs, nrow, integer(1)), c(1L, 0L, 0L))
  expect_equal(vapply(out$pathways, nrow, integer(1)), c(1L, 0L, 0L))
  expect_equal(out$goi_line, c("goi", NA, NA))
  expect_equal(out$dact_grounding, c("lit", "", ""))
  expect_equal(out$has_regulator_goi, c(TRUE, FALSE, FALSE))
  expect_match(out$summary[2], "not present above the abundance threshold")
  expect_match(out$summary[3], "returned no results")
})

testthat::test_that("tables without the presence column pass through unchanged", {
  inp <- dplyr::select(.impact_results(), -present_above_thresh)
  expect_identical(mark_not_assessed_cell_types(inp), inp)
})
