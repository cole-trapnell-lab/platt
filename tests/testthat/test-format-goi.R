# Tests for the impact-table genes-of-interest formatter and the magnitude-arrow
# re-encoding of LLM gene lists (see repos/platt/R/narrative.R). These functions
# are internal; source the file into an isolated env (mirrors test-gene-filtering.R)
# but parent to globalenv so dplyr/tibble/etc. are visible.
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr); library(tidyr)
})

narr_env <- new.env(parent = globalenv())
if (!exists("%||%", envir = narr_env, inherits = TRUE)) {
  assign("%||%", rlang::`%||%`, envir = narr_env)
}
sys.source(testthat::test_path("..", "..", "R", "narrative.R"), envir = narr_env)
magnitude_arrow        <- narr_env$magnitude_arrow
format_goi             <- narr_env$format_goi
build_deg_arrow_lookup <- narr_env$build_deg_arrow_lookup
reencode_gene_arrows   <- narr_env$reencode_gene_arrows
reencode_pathway_arrows <- narr_env$reencode_pathway_arrows

test_that("magnitude_arrow buckets |log2FC| into 1/2/3 arrows with direction", {
  expect_equal(
    magnitude_arrow(
      c(0.5, 1.5, 3, -0.5, -2, NA),
      c("Overexpressed", "Overexpressed", "Overexpressed",
        "Underexpressed", "Underexpressed", NA)
    ),
    c("↑", "↑↑", "↑↑↑",
      "↓", "↓↓↓", "")
  )
})

test_that("magnitude_arrow takes direction from LFC sign when dysreg_type is NA", {
  expect_equal(magnitude_arrow(2.0), "↑↑↑")
  expect_equal(magnitude_arrow(-1.0), "↓↓")
  expect_equal(magnitude_arrow(0), "")        # zero -> no direction
  expect_equal(magnitude_arrow(NA_real_), "") # unknown -> no arrow
})

make_degs <- function() {
  tibble(
    gene_short_name = c("tfA", "tfB", "tfC", "effX", "effY", "effZ", "effW"),
    perturb_to_ctrl_shrunken_lfc = c(0.4, -2.5, 1.2, 3.1, -0.2, 1.8, -4.0),
    dysreg_type = c("Overexpressed", "Underexpressed", "Overexpressed",
                    "Overexpressed", "Underexpressed", "Overexpressed", "Underexpressed")
  )
}

test_that("format_goi puts regulators first, sorts by |LFC| desc, truncates", {
  g <- format_goi(make_degs(), genes_of_interest = c("tfA", "tfB", "tfC"),
                  max_per_compartment = 3)
  # regulators (by |LFC|): tfB 2.5, tfC 1.2, tfA 0.4; then others: effW 4.0,
  # effX 3.1, effZ 1.8, and effY 0.2 truncated -> [+1 more]
  expect_equal(
    g$goi_line,
    paste("tfB↓↓↓", "tfC↑↑", "tfA↑",
          "effW↓↓↓", "effX↑↑↑", "effZ↑↑",
          "[+1 more]", sep = ", ")
  )
})

test_that("format_goi returns 'none' for empty / missing input", {
  expect_equal(format_goi(make_degs()[0, ], c("tfA"))$goi_line, "none")
  expect_equal(format_goi(NULL)$goi_line, "none")
})

test_that("reencode_gene_arrows rewrites to magnitude arrows, sorted, from DEG LFCs", {
  lk <- build_deg_arrow_lookup(make_degs())
  # LLM emitted stale directions; re-encoding trusts the DEG table LFC.
  expect_equal(
    reencode_gene_arrows(c("tfa ↑", "effz ↓", "tfb ↓"), lk),
    c("tfb↓↓↓", "effz↑↑", "tfa↑")
  )
})

test_that("reencode_pathway_arrows re-encodes the dysregulated_genes column", {
  lk <- build_deg_arrow_lookup(make_degs())
  tb <- tibble(name = "p1", description = "d",
               dysregulated_genes = "tfa ↑, effw ↓")
  expect_equal(reencode_pathway_arrows(tb, lk)$dysregulated_genes,
               "effw↓↓↓, tfa↑")
})
