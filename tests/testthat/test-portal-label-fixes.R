# Labels the portal and its narrative prompts read: platt#66, #68, #58(a).

.fgsea <- function(path, NES, padj = 0.001) {
  tibble::tibble(path = path, NES = NES, padj = padj)
}

test_that("F3 severity is graded from its stress modules, not the 'mild' fallback (platt#68)", {
  # The F3 label carries its module names, so an exact sev_map lookup never
  # matched and every F3 row came back "mild".
  sev <- classify_fitness(.fgsea("p53", 2.4))
  expect_match(sev$fitness_label, "^F3 Stress/toxicity \\(p53\\)$")
  expect_equal(unname(sev$severity), "severe")

  mod <- classify_fitness(.fgsea(c("Hypoxia", "UPR_ER"), c(1.7, 1.6)))
  expect_match(mod$fitness_label, "^F3 Stress/toxicity \\(Hypoxia, UPR_ER\\)$")
  expect_equal(unname(mod$severity), "moderate")
})

test_that("labels that match sev_map exactly keep their grading", {
  out <- classify_fitness(.fgsea(c("Proliferation_S", "p53"), c(2.1, 1.6)))
  expect_equal(unname(out$severity[out$fitness_label == "F1 Proliferation increase"]), "severe")
  expect_equal(unname(out$severity[startsWith(out$fitness_label, "F3 ")]), "moderate")
})

test_that("no fitness phenotype has severity 'none', not 'mild'", {
  out <- classify_fitness(.fgsea("p53", 0.4, padj = 0.9))
  expect_equal(out$fitness_label, "F0 No significant phenotype")
  expect_equal(unname(out$severity), "none")
})


# Two timepoints tied at q = 1: one precise and small (certifies A0), one too
# noisy to certify anything (AU).
.tied_dact <- function() {
  tibble::tibble(
    cell_group           = "ct",
    timepoint_x          = c(24, 48),
    present_above_thresh = TRUE,
    percent_max_abund    = c(0.4, 0.9),
    delta_log_abund      = c(0.01, 0.02),
    delta_log_abund_se   = c(0.05, 1.50),
    delta_q_value        = 1,
    df_resid             = 40
  )
}

.code_of <- function(row) {
  assign_abundance_code(row$delta_log_abund, row$delta_q_value,
    se = row$delta_log_abund_se, df = row$df_resid
  )
}

test_that("the fixture's tied rows really carry different codes", {
  codes <- .code_of(.tied_dact())
  expect_setequal(codes, c("A0 No change", "AU Undetermined"))
})

test_that("q ties no longer make the code depend on row order (platt#66)", {
  d <- .tied_dact()
  fwd <- dacts_when_abundant(d, percent_max_thresh = 0)
  rev <- dacts_when_abundant(d[2:1, ], percent_max_thresh = 0)
  expect_equal(nrow(fwd), 1L)
  expect_identical(fwd, rev)
  # A null is certified only when every tied timepoint certifies it.
  expect_equal(.code_of(fwd), "AU Undetermined")
})

test_that("A0 survives a tie when every tied row certifies it", {
  d <- .tied_dact()
  d$delta_log_abund_se <- c(0.05, 0.06)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(.code_of(got), "A0 No change")
  # Remaining ties go to the timepoint of peak wild-type abundance.
  expect_equal(got$timepoint_x, 48)
  expect_identical(got, dacts_when_abundant(d[2:1, ], percent_max_thresh = 0))
})

test_that("a unique minimum q is still the row kept", {
  d <- .tied_dact()
  d$delta_q_value <- c(0.5, 0.9)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 24)
})

test_that("with_ties = TRUE still returns every tied row", {
  expect_equal(nrow(dacts_when_abundant(.tied_dact(), percent_max_thresh = 0, with_ties = TRUE)), 2L)
})

test_that("ties degrade to AU, deterministically, on tables without se/df", {
  d <- .tied_dact()[, setdiff(names(.tied_dact()), c("delta_log_abund_se", "df_resid"))]
  fwd <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_identical(fwd, dacts_when_abundant(d[2:1, ], percent_max_thresh = 0))
  expect_equal(fwd$timepoint_x, 48)
})


test_that("the narrative prompt describes abundance tiers as sizes (platt#58a)", {
  expect_equal(
    abundance_magnitude_words(c("none", "mild", "moderate", "Severe ", NA)),
    c("none", "small", "medium", "large", NA)
  )
  expect_equal(abundance_magnitude_words("something else"), "something else")
})
