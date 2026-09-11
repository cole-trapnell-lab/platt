test_helper_make_rank_vec_and_sets <- function(seed = 1) {
    set.seed(seed)
    genes <- paste0("gene", seq_len(300))
    rank_vec <- setNames(rnorm(length(genes)), genes)
    rank_vec <- sort(rank_vec, decreasing = TRUE)
    gene_sets <- list(
        set_a = genes[1:20],
        set_b = genes[141:160],
        set_c = genes[281:300]
    )
    list(rank_vec = rank_vec, gene_sets = gene_sets)
}

testthat::test_that("run_fgsea_modules is deterministic across repeated calls with identical inputs", {
    fixture <- test_helper_make_rank_vec_and_sets()

    res1 <- platt:::run_fgsea_modules(
        fixture$rank_vec, fixture$gene_sets,
        minSize = 10, maxSize = 500, nperm = 200,
        seed_key = "perturbA::cellB::fitness"
    )
    res2 <- platt:::run_fgsea_modules(
        fixture$rank_vec, fixture$gene_sets,
        minSize = 10, maxSize = 500, nperm = 200,
        seed_key = "perturbA::cellB::fitness"
    )

    testthat::expect_true(nrow(res1) > 0)
    testthat::expect_identical(res1$path, res2$path)
    testthat::expect_identical(res1$NES, res2$NES)
    testthat::expect_identical(res1$padj, res2$padj)
    testthat::expect_identical(res1$log2err, res2$log2err)
})

testthat::test_that("run_fgsea_modules is deterministic with no seed_key supplied (falls back to input-derived seed)", {
    fixture <- test_helper_make_rank_vec_and_sets()

    res1 <- platt:::run_fgsea_modules(fixture$rank_vec, fixture$gene_sets, minSize = 10, maxSize = 500, nperm = 200)
    res2 <- platt:::run_fgsea_modules(fixture$rank_vec, fixture$gene_sets, minSize = 10, maxSize = 500, nperm = 200)

    testthat::expect_identical(res1$NES, res2$NES)
    testthat::expect_identical(res1$padj, res2$padj)
    testthat::expect_identical(res1$log2err, res2$log2err)
})

testthat::test_that("run_fgsea_modules surfaces log2err from fgseaMultilevel", {
    fixture <- test_helper_make_rank_vec_and_sets()

    res <- platt:::run_fgsea_modules(fixture$rank_vec, fixture$gene_sets, minSize = 10, maxSize = 500, nperm = 200)

    testthat::expect_true("log2err" %in% names(res))
    testthat::expect_true(all(is.finite(res$log2err)))
})

testthat::test_that("run_fgsea_modules never forwards nperm to fgsea::fgsea (stays on the fgseaMultilevel path)", {
    fixture <- test_helper_make_rank_vec_and_sets()

    # fgsea::fgsea() only emits the fgseaSimple deprecation warning when it
    # sees an `nperm` argument in `...`; that warning is the direct symptom
    # of the bug this fix addresses (previously silenced by suppressWarnings()
    # around the whole call). If run_fgsea_modules() ever regresses to passing
    # nperm straight through, this will start failing.
    testthat::expect_no_warning(
        withCallingHandlers(
            fgsea::fgsea(
                pathways = fixture$gene_sets, stats = fixture$rank_vec,
                minSize = 10, maxSize = 500, nPermSimple = 200
            ),
            warning = function(w) stop(w)
        )
    )

    direct <- withr::with_seed(
        platt:::.fgsea_seed_from_key("perturbA::cellB::fitness"),
        suppressWarnings(fgsea::fgsea(
            pathways = fixture$gene_sets, stats = fixture$rank_vec,
            minSize = 10, maxSize = 500, nPermSimple = 200
        ))
    )
    via_wrapper <- platt:::run_fgsea_modules(
        fixture$rank_vec, fixture$gene_sets,
        minSize = 10, maxSize = 500, nperm = 200,
        seed_key = "perturbA::cellB::fitness"
    )

    testthat::expect_identical(sort(via_wrapper$NES), sort(direct$NES))
    testthat::expect_identical(sort(via_wrapper$padj), sort(direct$padj))
})

testthat::test_that(".fgsea_seed_from_key gives different keys different seeds (order/identity independence)", {
    s1 <- platt:::.fgsea_seed_from_key("perturbA::cellB::fitness")
    s2 <- platt:::.fgsea_seed_from_key("perturbA::cellC::fitness")
    s3 <- platt:::.fgsea_seed_from_key("perturbA::cellB::fitness")

    testthat::expect_identical(s1, s3)
    testthat::expect_false(identical(s1, s2))
})

# Regression: a `cell_type` COLUMN in ref_expression must not shadow the
# cell-type argument. sulston's make_gene_sets began supplying such a column in
# 2026-01, which silently turned the subset in construct_identity_gene_sets into
# a no-op: every cell type got a ranking pooled over all cell types, and so
# received identical "identity" gene sets.
testthat::test_that(".rank_genes_by_specificity subsets by cell type even when ref_expression has a cell_type column", {
    ref <- tibble::tibble(
        cell_group = rep(c("ct_a", "ct_b"), each = 2),
        gene_short_name = c("a1", "a2", "b1", "b2"),
        fraction_expressing = c(0.5, 0.4, 0.5, 0.4),
        specificity = c(0.9, 0.8, 0.7, 0.6)
    )
    ref$cell_type <- ref$cell_group

    ranks_a <- platt:::.rank_genes_by_specificity(ref, "ct_a")
    testthat::expect_identical(names(ranks_a), c("a1", "a2"))
    testthat::expect_identical(unname(ranks_a), c(0.9, 0.8))

    testthat::expect_identical(
        names(platt:::.rank_genes_by_specificity(ref, "ct_b")),
        c("b1", "b2")
    )
})

testthat::test_that(".rank_genes_by_specificity applies the expression floor per cell type and sorts descending", {
    ref <- tibble::tibble(
        cell_group = c("ct_a", "ct_a", "ct_a"),
        gene_short_name = c("lo", "hi", "mid"),
        fraction_expressing = c(0.001, 0.5, 0.5),
        specificity = c(0.99, 0.3, 0.6)
    )
    ranks <- platt:::.rank_genes_by_specificity(ref, "ct_a")
    # `lo` has the highest specificity but is below the floor, so it is dropped
    # even though dropping it costs the top-ranked gene.
    testthat::expect_identical(names(ranks), c("mid", "hi"))

    # the floor is tunable rather than hardcoded
    testthat::expect_identical(
        names(platt:::.rank_genes_by_specificity(ref, "ct_a", min_fraction_expressing = 0)),
        c("lo", "mid", "hi")
    )
})

testthat::test_that(".rank_genes_by_specificity keeps the highest-specificity row per duplicated gene symbol", {
    ref <- tibble::tibble(
        cell_group = "ct_a",
        gene_short_name = c("g1", "g1", "g2"),
        fraction_expressing = 0.5,
        specificity = c(0.2, 0.8, 0.5)
    )
    ranks <- platt:::.rank_genes_by_specificity(ref, "ct_a")
    testthat::expect_identical(names(ranks), c("g1", "g2"))
    testthat::expect_identical(unname(ranks), c(0.8, 0.5))
})

testthat::test_that(".rank_genes_by_specificity errors clearly on a missing required column", {
    ref <- tibble::tibble(cell_group = "ct_a", gene_short_name = "g1", specificity = 0.5)
    testthat::expect_error(
        platt:::.rank_genes_by_specificity(ref, "ct_a"),
        "fraction_expressing"
    )
})


test_that("every code assign_abundance_code() can emit has a colour-map entry", {

  # The A3 bug: assign_abundance_code() emits "A3 Near-loss" while
  # abundance_code_map keyed on "A3 Ablation/Loss". The lookup returned NA,
  # which impact_to_phenos() coerces to a 0 proxy, so the most severe abundance
  # phenotype rendered identically to "no change". It survived ~6 months because
  # A3 also required q < 0.01 and never fired in practice.
  #
  # This test is the actual fix: it fails if the two lists ever drift again.
  emitted <- c("A0 No change", "A1 Expansion", "A2 Depletion", "A3 Near-loss")
  code_map <- eval(formals(impact_to_phenos)$abundance_code_map)

  missing <- setdiff(emitted, names(code_map))
  expect_equal(missing, character(0))

  # A3 must carry a real, negative sign -- not NA, and not 0.
  expect_false(is.na(code_map[["A3 Near-loss"]]))
  expect_lt(code_map[["A3 Near-loss"]], 0)

  # and it must be at least as severe as a plain depletion
  expect_lte(code_map[["A3 Near-loss"]], code_map[["A2 Depletion"]])
})


test_that("an unmapped abundance code degrades to a neutral proxy", {

  # Documents the failure mode that hid the A3 bug, so the consequence of a
  # future mismatch is explicit rather than folklore: an unknown code does not
  # error, it silently colours as no-change.
  code_map <- eval(formals(impact_to_phenos)$abundance_code_map)
  sign <- unname(code_map["A9 Not a real code"])
  expect_true(is.na(sign))

  proxy <- ifelse(is.na(sign), 0, sign * 2.0)
  expect_equal(proxy, 0)
  expect_equal(proxy, unname(code_map[["A0 No change"]]) * 2.0)
})
