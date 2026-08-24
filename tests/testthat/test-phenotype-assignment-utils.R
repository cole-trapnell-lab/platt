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
