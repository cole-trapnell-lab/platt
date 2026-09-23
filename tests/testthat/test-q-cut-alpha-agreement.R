# The abundance cascade and hooke's contrast table must answer "were we powered
# to see a change worth calling here?" at the SAME level.
#
# Two quantities, computed in two packages, are meant to be the same predicate:
#
#   hooke::compare_abundances()  power_status == "powered"
#                                  <=> mdfc80 <= margin_fold_change
#                                  computed at its `alpha` (default 0.05)
#
#   assign_abundance_code()      resolved  <=> abundance_mdfc80 <= margin
#                                  computed at `alpha = q_cut`
#
# They agree only while `q_cut` equals hooke's `alpha`. While they disagreed
# (q_cut 0.1 vs alpha 0.05) the same cell type could read "underpowered" in the
# contrast table and "A0 No change" in the impact table. 0.05 is the stricter of
# the two, so the disagreement ran one way only -- which is exactly why it went
# unnoticed: nothing ever over-claimed, it only under-reported power.
#
# `lfc_cut` cannot drift this way because the margin is DERIVED from it
# (margin_fold_change = exp(lfc_cut)). Alpha is derived nowhere, so it needs a test.

test_that("the cascade's q_cut matches the alpha hooke computes detectability at", {
  skip_if_not_installed("hooke")

  expect_identical(
    formals(assign_abundance_code)$q_cut,
    formals(hooke::compare_abundances)$alpha
  )
  expect_identical(
    formals(assign_abundance_severity)$q_cut,
    formals(hooke::compare_abundances)$alpha
  )
  expect_identical(
    formals(assign_phenotypes)$abundance_q_cut,
    formals(hooke::compare_abundances)$alpha
  )
})

test_that("abundance_mdfc80() reproduces hooke's mdfc80 at the shared alpha", {
  skip_if_not_installed("hooke")

  # Same standard errors and df, computed independently in each package. The
  # cascade sees a single combined SE; hooke combines two. Feeding hooke the
  # per-condition SEs that sum to ours makes the two directly comparable.
  se_combined <- c(0.05, 0.10, 0.1726, 0.20, 0.35, 0.60)
  se_each <- se_combined / sqrt(2)
  df <- 41

  alpha <- formals(hooke::compare_abundances)$alpha
  power <- formals(hooke::compare_abundances)$power

  ours <- abundance_mdfc80(se_combined, df, alpha = alpha, power = power)
  theirs <- hooke:::calculate_mdfc(se_each, se_each, alpha = alpha,
                                   power = power, df = df, base = exp(1))

  expect_equal(ours, theirs, tolerance = 1e-9)
})

test_that("resolved and hooke's power_status pick out the same rows", {
  skip_if_not_installed("hooke")

  se_combined <- c(0.05, 0.10, 0.1726, 0.20, 0.35, 0.60)
  se_each <- se_combined / sqrt(2)
  df <- 41
  lfc_cut <- 0.5
  margin <- exp(lfc_cut)

  alpha <- formals(hooke::compare_abundances)$alpha
  power <- formals(hooke::compare_abundances)$power

  # platt: a non-significant row resolves to A0 exactly when it is powered.
  code <- assign_abundance_code(
    change_when_present = rep(0.02, length(se_combined)),
    change_when_present_q_val = rep(0.9, length(se_combined)),
    se = se_combined, df = df, lfc_cut = lfc_cut
  )
  platt_powered <- code == "A0 No change"

  # hooke: the same question, via the effect-independent power at the margin.
  # The threshold is inlined rather than taken from hooke::calculate_power_status(),
  # which arrived after some installed hooke builds. `power_at_margin >= power`
  # is that function's whole body and is the documented dual of
  # `mdfc80 <= margin_fold_change`.
  pam <- hooke:::calculate_power_at_margin(se_each, se_each, alpha = alpha,
                                           margin_log = lfc_cut, df = df)
  hooke_powered <- !is.na(pam) & pam >= power

  expect_identical(platt_powered, hooke_powered)
  # Guard against a vacuous pass.
  expect_true(any(platt_powered))
  expect_true(any(!platt_powered))
})


test_that("A0 is never claimed for a row whose own estimate clears lfc_cut", {

  # `resolved` tests the DETECTION LIMIT. On its own that certifies a null for
  # any tightly measured row, including one whose point estimate is above the
  # threshold we would have called. A real row from v3.1.0 (WntC59, mesenchymal
  # cell of the meninx): a 2.06-fold expansion that misses q_cut, with an
  # mdfc80 of 1.639 against a 1.649 margin -- resolved by 0.6%.
  lfc <- 0.721; q <- 0.094; se <- 0.162; df <- 12

  expect_lte(abundance_mdfc80(se, df, alpha = 0.05), exp(0.5))   # it IS precise
  expect_gt(abs(lfc), 0.5)                                       # and NOT small

  expect_equal(
    assign_abundance_code(lfc, q, se = se, df = df),
    "AU Undetermined"
  )

  # The guard must not swallow genuine resolved nulls: same precision, small
  # estimate, still A0.
  expect_equal(
    assign_abundance_code(0.02, q, se = se, df = df),
    "A0 No change"
  )

  # Both directions, and right at the boundary.
  expect_equal(assign_abundance_code(-0.721, q, se = se, df = df), "AU Undetermined")
  expect_equal(assign_abundance_code(0.5, q, se = se, df = df), "AU Undetermined")
  expect_equal(assign_abundance_code(0.4999, q, se = se, df = df), "A0 No change")
})
