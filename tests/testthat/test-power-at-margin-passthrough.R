# Regression guard for a silent failure: `plot_phenotypes_glyphs()` keys
# `power_status` on `power_at_margin`, but the augment step that builds its
# input used to join only `power` (observed_power). The impact table carries no
# detectability columns, so nothing else supplied one and every cell type fell
# through to "Underpowered" -- a plot-wide understatement with no error.

test_helper_power_tbl <- function() {
    tibble::tibble(
        cell_group = c("tight", "loose", "degenerate"),
        delta_log_abund = c(0.9, 1.4, 0.2),
        delta_q_value = c(0.001, 0.02, 0.9),
        power = c(0.99, 0.95, 0.05),
        # Only `power_at_margin` reflects precision: "loose" has a big effect
        # measured badly, so observed_power is high and power_at_margin is not.
        power_at_margin = c(0.96, 0.11, NA_real_),
        margin_fold_change = rep(exp(0.5), 3)
    )
}

testthat::test_that(".augment_phenos_with_abundance_summary carries power_at_margin through the join", {
    power_tbl <- test_helper_power_tbl()
    phenos <- tibble::tibble(
        cell_group = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )
    impact_table <- tibble::tibble(
        cell_type = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )

    out <- platt:::.augment_phenos_with_abundance_summary(
        phenos, power_tbl, impact_table, powered_thresh = 0.8
    )

    testthat::expect_true(all(c("power_at_margin", "margin_fold_change", "powered_thresh") %in% names(out)))
    testthat::expect_equal(
        out$power_at_margin[match(c("tight", "loose", "degenerate"), out$cell_group)],
        c(0.96, 0.11, NA_real_)
    )
    testthat::expect_equal(unique(out$margin_fold_change), exp(0.5))
})

testthat::test_that("power_status keys on precision, not on the observed effect", {
    power_tbl <- test_helper_power_tbl()
    phenos <- tibble::tibble(
        cell_group = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )
    impact_table <- tibble::tibble(
        cell_type = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )

    out <- platt:::.augment_phenos_with_abundance_summary(
        phenos, power_tbl, impact_table, powered_thresh = 0.8
    )
    status <- ifelse(
        !is.na(out$power_at_margin) & out$power_at_margin >= out$powered_thresh,
        "Powered", "Underpowered"
    )
    names(status) <- out$cell_group

    testthat::expect_equal(unname(status[["tight"]]), "Powered")
    # Would be "Powered" if this still keyed on `power` (0.95).
    testthat::expect_equal(unname(status[["loose"]]), "Underpowered")
    # NA (degenerate fit, insufficient df, or column absent) is NOT powered.
    testthat::expect_equal(unname(status[["degenerate"]]), "Underpowered")
})

testthat::test_that("a power_tbl predating Hooke 0.0.3 yields NA rather than an error", {
    power_tbl <- test_helper_power_tbl() %>%
        dplyr::select(-power_at_margin, -margin_fold_change)
    phenos <- tibble::tibble(
        cell_group = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )
    impact_table <- tibble::tibble(
        cell_type = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )

    out <- platt:::.augment_phenos_with_abundance_summary(
        phenos, power_tbl, impact_table, powered_thresh = 0.8
    )

    testthat::expect_true("power_at_margin" %in% names(out))
    testthat::expect_true(all(is.na(out$power_at_margin)))
})

testthat::test_that("a stored power_status column survives the join and is preferred", {
    # Hooke >= 0.0.5 writes power_status. The figure should show what the table
    # asserts rather than a second opinion recomputed here, so the stored label
    # wins even where a local recompute would disagree.
    power_tbl <- test_helper_power_tbl()
    power_tbl$power_status <- c("Powered", "Underpowered", "Underpowered")
    # Deliberately inconsistent with power_at_margin (0.96 >= 0.8) to prove
    # which one the pipeline reads.
    power_tbl$power_status[1] <- "Underpowered"

    phenos <- tibble::tibble(
        cell_group = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )
    impact_table <- tibble::tibble(
        cell_type = power_tbl$cell_group,
        abundance_code = "A0 No change"
    )

    out <- platt:::.augment_phenos_with_abundance_summary(
        phenos, power_tbl, impact_table, powered_thresh = 0.8
    )

    testthat::expect_true("power_status" %in% names(out))
    testthat::expect_equal(
        out$power_status[match("tight", out$cell_group)],
        "Underpowered"
    )
})
