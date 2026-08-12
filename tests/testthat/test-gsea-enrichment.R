test_helper_make_state_gsea_fixture <- function(seed = 1) {
    set.seed(seed)
    genes <- paste0("gene", seq_len(300))
    gene_patterns_within_state_graph <- data.frame(gene_short_name = genes)
    gene_patterns_within_state_graph$pattern_activity_score <- matrix(rnorm(length(genes)), ncol = 1)
    gene_df <- data.frame(
        gene_short_name = c(genes[1:20], genes[141:160], genes[281:300]),
        gs_name = rep(c("set_a", "set_b", "set_c"), each = 20)
    )
    list(gene_patterns_within_state_graph = gene_patterns_within_state_graph, gene_df = gene_df)
}

testthat::test_that("calc_gsea_enrichment_on_state_specific_genes is deterministic across repeated calls", {
    fixture <- test_helper_make_state_gsea_fixture()

    res1 <- platt:::calc_gsea_enrichment_on_state_specific_genes(
        fixture$gene_patterns_within_state_graph, fixture$gene_df, sig_thresh = 1.0
    )
    res2 <- platt:::calc_gsea_enrichment_on_state_specific_genes(
        fixture$gene_patterns_within_state_graph, fixture$gene_df, sig_thresh = 1.0
    )

    testthat::expect_true(nrow(res1) > 0)
    testthat::expect_identical(res1$NES, res2$NES)
    testthat::expect_identical(res1$padj, res2$padj)
    testthat::expect_identical(res1$log2err, res2$log2err)
})

testthat::test_that("calc_gsea_enrichment_on_state_specific_genes surfaces log2err", {
    fixture <- test_helper_make_state_gsea_fixture()

    res <- platt:::calc_gsea_enrichment_on_state_specific_genes(
        fixture$gene_patterns_within_state_graph, fixture$gene_df, sig_thresh = 1.0
    )

    testthat::expect_true("log2err" %in% names(res))
    testthat::expect_true(all(is.finite(res$log2err)))
})

testthat::test_that("calc_gsea_enrichment_on_state_specific_genes respects an explicit seed_key", {
    fixture <- test_helper_make_state_gsea_fixture()

    res1 <- platt:::calc_gsea_enrichment_on_state_specific_genes(
        fixture$gene_patterns_within_state_graph, fixture$gene_df,
        sig_thresh = 1.0, seed_key = "stateA::patternX"
    )
    res2 <- platt:::calc_gsea_enrichment_on_state_specific_genes(
        fixture$gene_patterns_within_state_graph, fixture$gene_df,
        sig_thresh = 1.0, seed_key = "stateA::patternX"
    )

    testthat::expect_identical(res1$NES, res2$NES)
    testthat::expect_identical(res1$padj, res2$padj)
})
