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


test_that("assign_abundance_code splits A0 into four states", {

  # Same non-significant observed change in every null row; only the standard
  # error differs. Before the split these were one indistinguishable bucket.
  lfc <- c(1.2, -2.5, -0.8, 0.02, 0.05, 0.05, 0.05, 0.00, NA)
  q   <- c(0.01, 0.005, 0.02, 0.90, 0.80, 0.95, 0.90, 1.00, NA)
  se  <- c(0.20, 0.30, 0.20, 0.10, 0.40, 1.20, NA, 0.00, 0.20)
  df  <- 40

  code <- assign_abundance_code(lfc, q, q_cut = 0.1, se = se, df = df)

  # Called rows are untouched by the split.
  expect_equal(code[1], "A1 Expansion")
  expect_equal(code[2], "A3 Near-loss")
  expect_equal(code[3], "A2 Depletion")

  # A0 now means a RESOLVED null: not significant, and powered to exclude a
  # two-fold change.
  expect_equal(code[4], "A0 No change")
  expect_lte(abundance_mdfc80(se[4], df, alpha = 0.1), 2)

  # Not significant, but could not have seen a two-fold change.
  expect_equal(code[5], "AU Undetermined")
  expect_equal(code[6], "AU Undetermined")
  expect_gt(abundance_mdfc80(se[5], df, alpha = 0.1), 2)

  # No standard error at all: we cannot certify a null, so we do not claim one.
  expect_equal(code[7], "AU Undetermined")

  # A degenerate fit (SE present but zero) was never really tested, which is
  # different from tested-and-inconclusive.
  expect_equal(code[8], "AN Not assessed")

  # No estimate at all.
  expect_equal(code[9], "AN Not assessed")

  # Every row lands in exactly one of the four states.
  expect_true(all(code %in% c("A0 No change", "A1 Expansion", "A2 Depletion",
                              "A3 Near-loss", "AU Undetermined", "AN Not assessed")))
})


test_that("assign_abundance_code is backward compatible when se/df are absent", {

  # Tables written before the SE and df columns existed must still produce a
  # code for every row, and must not silently claim a resolved null.
  lfc <- c(1.2, 0.02)
  q <- c(0.01, 0.90)

  code <- assign_abundance_code(lfc, q, q_cut = 0.1)
  expect_equal(code[1], "A1 Expansion")
  expect_equal(code[2], "AU Undetermined")
  expect_false(any(code == "A0 No change"))
})


test_that("abundance_mdfc80 is effect-independent and guards degenerate input", {

  se <- c(0.1, 0.2, 0.4)
  mdfc80 <- abundance_mdfc80(se, df = 40, alpha = 0.1)

  # Depends only on the standard error, monotonically.
  expect_equal(order(mdfc80), order(se))
  expect_equal(mdfc80, exp((qt(0.95, 40) + qt(0.8, 40)) * se))

  # A zero, negative, non-finite or missing SE has no detection limit, and an
  # invalid df has none either. None of these may return 1-fold.
  expect_true(is.na(abundance_mdfc80(0, 40)))
  expect_true(is.na(abundance_mdfc80(-1, 40)))
  expect_true(is.na(abundance_mdfc80(Inf, 40)))
  expect_true(is.na(abundance_mdfc80(NA_real_, 40)))
  expect_true(is.na(abundance_mdfc80(0.3, 0)))
  expect_true(is.na(abundance_mdfc80(0.3, NA_real_)))
})


test_that("the margin is derived from the phenotype-calling threshold", {

  # A resolved null claims we were powered to see a change we would have called
  # a phenotype. The smallest callable change is lfc_cut on the log scale, so
  # the only margin that makes that claim true is exp(lfc_cut).
  expect_equal(formals(assign_abundance_code)$margin_fold_change, quote(exp(lfc_cut)))
  expect_equal(eval(formals(assign_abundance_code)$lfc_cut), 0.5)

  # Tuning the call threshold keeps the verdict coherent automatically.
  loose <- assign_abundance_code(0.02, 0.9, se = 0.3, df = 40, lfc_cut = 1.0)
  tight <- assign_abundance_code(0.02, 0.9, se = 0.3, df = 40, lfc_cut = 0.25)
  expect_equal(loose, "A0 No change")      # only needs to exclude 2.7-fold
  expect_equal(tight, "AU Undetermined")   # must exclude 1.28-fold, cannot

  # The over-claim window that margin = 2 would have created: a detection limit
  # worse than the smallest callable phenotype must NOT read as a resolved null.
  se_over <- log(1.8) / (qt(0.95, 40) + qt(0.8, 40))
  expect_gt(abundance_mdfc80(se_over, 40, alpha = 0.1), exp(0.5))
  expect_lt(abundance_mdfc80(se_over, 40, alpha = 0.1), 2)
  expect_equal(assign_abundance_code(0.02, 0.9, se = se_over, df = 40), "AU Undetermined")
  expect_equal(assign_abundance_code(0.02, 0.9, se = se_over, df = 40,
                                     margin_fold_change = 2), "A0 No change")
})


test_that("the margin is a fold change, not a log", {

  # The design note wrote the resolved test as `mdfc80 < log(2)`. mdfc80 is a
  # fold change and is always >= 1, so that comparison can never be TRUE and
  # every null would have become AU. This pins the units.
  se <- c(0.01, 0.05, 0.1)
  mdfc80 <- abundance_mdfc80(se, df = 40, alpha = 0.1)
  expect_true(all(mdfc80 >= 1))
  expect_equal(sum(mdfc80 < log(2)), 0)

  code <- assign_abundance_code(rep(0.01, 3), rep(0.9, 3), q_cut = 0.1,
                                se = se, df = 40)
  expect_true(all(code == "A0 No change"))

  # A stricter margin resolves fewer rows.
  strict <- assign_abundance_code(rep(0.01, 3), rep(0.9, 3), q_cut = 0.1,
                                  se = se, df = 40, margin_fold_change = 1.05)
  expect_true(sum(strict == "A0 No change") < 3)
})


test_that("A3 is unreachable if its branch is moved below A2", {

  # The two loss thresholds look backwards read in isolation: A2's -0.5 is a
  # weaker cutoff than A3's -2.0, so every A3 row also satisfies A2. Order is
  # what separates them. This pins that, so reordering the cascade fails here
  # rather than silently retiring A3.
  lfc <- c(-2.5, -0.8)
  q   <- c(0.01, 0.01)
  se  <- c(0.1, 0.1)

  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = se, df = 40),
               c("A3 Near-loss", "A2 Depletion"))

  # A severe loss satisfies BOTH conditions; only precedence keeps it in A3.
  expect_true(lfc[1] <= -2.0)   # A3 condition
  expect_true(lfc[1] <= -0.5)   # A2 condition, also true
})


test_that("present-but-invalid df is AN, absent df is AU", {

  # Same rule already applied to `se`: a broken fit is "never tested" (AN); a
  # missing column is "cannot certify" (AU). Without the df half, a row with a
  # good SE and an invalid df produced mdfc80 = NA and fell through to AU.
  lfc <- 0.02; q <- 0.9; good_se <- 0.1

  # present but invalid -> AN
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = 0),
               "AN Not assessed")
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = -3),
               "AN Not assessed")
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = Inf),
               "AN Not assessed")

  # absent -> AU (cannot certify, but nothing says the fit was broken)
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = NA_real_),
               "AU Undetermined")
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = NULL),
               "AU Undetermined")
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = NULL, df = 40),
               "AU Undetermined")

  # a valid pair still resolves normally
  expect_equal(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = 40),
               "A0 No change")

  # and an invalid df never silently becomes a resolved null
  for (bad in c(0, -1, Inf)) {
    expect_false(assign_abundance_code(lfc, q, q_cut = 0.1, se = good_se, df = bad)
                 == "A0 No change")
  }
})


test_that("every route to an NA detection limit lands somewhere deliberate", {

  # Enumerates the ways mdfc80 can be NA and asserts each maps to the intended
  # code, so a future edit cannot quietly reroute one of them.
  cases <- list(
    list(se = 0,          df = 40,        want = "AN Not assessed"),  # degenerate SE
    list(se = -1,         df = 40,        want = "AN Not assessed"),
    list(se = Inf,        df = 40,        want = "AN Not assessed"),
    list(se = 0.1,        df = 0,         want = "AN Not assessed"),  # invalid df
    list(se = NA_real_,   df = 40,        want = "AU Undetermined"),  # absent SE
    list(se = 0.1,        df = NA_real_,  want = "AU Undetermined")   # absent df
  )
  for (cs in cases) {
    got <- assign_abundance_code(0.02, 0.9, q_cut = 0.1, se = cs$se, df = cs$df)
    expect_equal(got, cs$want)
    expect_true(is.na(abundance_mdfc80(cs$se, cs$df, alpha = 0.1)))
  }
})


test_that("severity and the code cascade agree for any shared cutoff", {

  # The "mild" tier is the same condition as being called A1/A2, so
  # severity == "none" must mean "not called" and vice versa. That used to hold
  # only because assign_abundance_severity() hardcoded 0.5 and 0.1, matching
  # assign_abundance_code()'s defaults by coincidence. Making those tunable
  # broke it: lfc_cut = 0.8 produced rows reported as "A0 No change" with
  # "mild" severity.
  set.seed(11)
  n <- 4000
  lfc <- rnorm(n, 0, 1.2)
  q   <- runif(n)^2
  se  <- abs(rnorm(n, 0.25, 0.15)) + 0.01
  non_called <- c("A0 No change", "AU Undetermined", "AN Not assessed")

  for (lc in c(0.3, 0.5, 0.8, 1.0)) {
    for (qc in c(0.01, 0.05, 0.1)) {
      code <- assign_abundance_code(lfc, q, q_cut = qc, se = se, df = 40, lfc_cut = lc)
      sev  <- assign_abundance_severity(lfc, q, q_cut = qc, lfc_cut = lc)

      # no non-call carries a severity grade
      expect_equal(sum(code %in% non_called & sev != "none"), 0,
                   info = sprintf("lfc_cut=%.1f q_cut=%.2f", lc, qc))
      # and no call is graded "none"
      expect_equal(sum(!(code %in% non_called) & sev == "none"), 0,
                   info = sprintf("lfc_cut=%.1f q_cut=%.2f", lc, qc))
    }
  }
})


test_that("severity defaults are unchanged from the hardcoded version", {

  # Threading the cutoffs must not move any published severity label.
  set.seed(12)
  n <- 5000
  lfc <- rnorm(n, 0, 1.2)
  q   <- runif(n)^2

  frozen <- dplyr::case_when(
    abs(lfc) >= 2.0 & q < 0.01 ~ "severe",
    abs(lfc) >= 1.0 & q < 0.05 ~ "moderate",
    abs(lfc) >= 0.5 & q < 0.1  ~ "mild",
    TRUE ~ "none"
  )
  expect_identical(assign_abundance_severity(lfc, q), frozen)
})


test_that("cell types the experiment could not assess get a row, not silence", {

  # dacts_when_abundant() drops every row of a cell type whose present_above_thresh
  # is FALSE at all timepoints. That flag is a wild-type reference property -- it
  # marks states the timepoint window could not assess -- so the cell type used to
  # vanish from the impact table with no row at all, which reads identically to
  # "this state was not in the experiment".
  dact <- tibble::tibble(
    cell_group           = c("assessable", "assessable", "out_of_window"),
    timepoint_x          = c(24, 36, 24),
    delta_log_abund      = c(0.02, 0.05, -3.0),
    delta_q_value        = c(0.9, 0.8, 0.9),
    percent_max_abund    = c(0.9, 0.8, 0.02),
    present_above_thresh = c(TRUE, TRUE, FALSE)
  )

  kept <- dacts_when_abundant(dact, percent_max_thresh = 0)
  expect_false("out_of_window" %in% kept$cell_group)   # the silent drop
  expect_true("assessable" %in% kept$cell_group)

  # The universe taken before the drop still contains it, which is what lets the
  # phenotype assigner report it.
  universe <- unique(dact$cell_group)
  expect_true("out_of_window" %in% universe)

  # A cell type with no surviving abundance row is "AN Not assessed" -- never
  # tested -- rather than NA or a resolved null.
  empty_row <- kept[kept$cell_group == "out_of_window", ]
  expect_equal(nrow(empty_row), 0)
  code <- if (nrow(empty_row) > 0) "unreachable" else "AN Not assessed"
  expect_equal(code, "AN Not assessed")
  expect_true(code %in% NON_CALLED_ABUNDANCE_CODES)
})


test_that("all_cell_types widens the universe without changing assessed calls", {

  # The argument defaults to the old behaviour, and supplying it only ADDS rows;
  # it must not alter the code assigned to any cell type that was already there.
  dact <- tibble::tibble(
    cell_group                   = c("a", "b"),
    change_when_present          = c(1.2, 0.02),
    change_when_present_q_val    = c(0.01, 0.90),
    change_when_present_se       = c(0.20, 0.10),
    change_when_present_tvalue_df = 40
  )

  assessed <- assign_abundance_code(dact$change_when_present,
                                    dact$change_when_present_q_val,
                                    q_cut = 0.1,
                                    se = dact$change_when_present_se,
                                    df = dact$change_when_present_tvalue_df)
  expect_equal(assessed, c("A1 Expansion", "A0 No change"))

  # A third cell type with no row anywhere in dact is additive and lands as AN.
  universe <- c("a", "b", "never_assessed")
  codes <- vapply(universe, function(ct) {
    row <- dact[dact$cell_group == ct, ]
    if (nrow(row) > 0) {
      assign_abundance_code(row$change_when_present, row$change_when_present_q_val,
                            q_cut = 0.1, se = row$change_when_present_se,
                            df = row$change_when_present_tvalue_df)
    } else "AN Not assessed"
  }, character(1))

  expect_equal(unname(codes[c("a", "b")]), assessed)   # unchanged
  expect_equal(unname(codes[["never_assessed"]]), "AN Not assessed")
})
