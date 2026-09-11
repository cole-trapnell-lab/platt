# The uncallability gate policy of annotate_empirical_fdr_model(). The tail null is stubbed out
# (a surface whose quantiles sit far below any observed |z|, so p_tail is ~0 and drops out of
# max(p_ashr, p_tail)) -- these tests are about which calls the GATES pin to 1, nothing else.
stub_model <- function() {
  taus  <- c(0.5, 0.9, 0.999)
  egrid <- c(-6, -3, 0); pegrid <- c(0, 1, 2, 3)
  Q <- matrix(-1e6, nrow = length(egrid) * length(pegrid), ncol = length(taus))
  list(surf = list(dn = Q, up = Q), taus = taus, egrid = egrid, pegrid = pegrid)
}

# One row per scenario. Each is strongly significant by ashr, so anything that ends at
# empirical_p == 1 was gated and anything below 0.05 was not.
gate_fixture <- function() {
  deg <- data.frame(
    id                   = c("well_det_small_arm", "trace_big_arm", "gain_from_zero",
                             "one_embryo", "absent_both", "complete_depletion"),
    gene_short_name      = c("well_det_small_arm", "trace_big_arm", "gain_from_zero",
                             "one_embryo", "absent_both", "complete_depletion"),
    cell_group           = "CT",
    log_mean_expression  = -3,
    perturb_to_ctrl_shrunken_lfc = c(-2, -2, 2, 2, -2, -2),
    perturb_to_ctrl_p_value      = 1e-8,
    z                    = c(-6, -6, 6, 6, -6, -6),
    stringsAsFactors     = FALSE)
  # arm_cells is per cell type, so every row shares one arm geometry: a 900-cell control arm
  # against a 9-cell perturbation arm. thin_arm_cells = 9 throughout.
  arm_cells <- data.frame(cell_group = "CT", n_ctrl_cells = 900, n_pert_cells = 9)
  counts <- data.frame(
    id = deg$id, gene_short_name = deg$id, cell_group = "CT",
    K_ctrl      = c(2000,  15,   0,    5, 0, 2000),
    K_pert      = c(  20,   0,  40,   40, 0,    0),
    # 38% control detection -> 9 * 0.38 = 3.42 expected in the thin arm: powered.
    # trace_big_arm is 1% -> 0.09 expected, despite its 900-cell control arm.
    # gain_from_zero has nothing in control and 3 of 9 perturbation cells -> 9 * 1/3 = 3.
    # one_embryo is powered the same way (3 of 9) but every detection sits in one embryo.
    n_expr_ctrl = c( 342,   9,   0,    1, 0,  342),
    n_expr_pert = c(   3,   0,   3,    3, 0,    0),
    n_emb_ctrl  = c(  20,   7,   0,    1, 0,   20),
    n_emb_pert  = c(   3,   0,   3,    1, 0,    0),
    stringsAsFactors = FALSE)
  support <- data.frame(cell_group = "CT", n_pert_pb = 9)
  list(deg = deg, arm_cells = arm_cells, counts = counts, support = support)
}

ep <- function(out, id) out$empirical_p[match(id, out$id)]

test_that("the default policy gates on power per call, not on arm size", {
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support)

  # A well-detected gene stays callable in a 9-cell arm: 9 cells x 38% = 3.42 expected >= 3.
  expect_lt(ep(out, "well_det_small_arm"), 0.05)
  # A trace gene is gated even though its control arm is large: 9 x 2% = 0.18 expected.
  expect_equal(ep(out, "trace_big_arm"), 1)
  # A gain from a zero baseline is testable, because the rate comes from the arm that has the
  # gene: 9 cells x (3/9) = 3 expected. The retired gain_from_zero rule killed this outright.
  expect_lt(ep(out, "gain_from_zero"), 0.05)
  # Powered, but every detection sits in a single embryo in either arm -> not_replicated.
  expect_equal(ep(out, "one_embryo"), 1)
  # No counts in either arm: caught by not_replicated (0 embryos), so not_detected_anywhere
  # is not needed to reach it.
  expect_equal(ep(out, "absent_both"), 1)
  # A complete depletion out of a well-replicated control arm survives: the thin arm was
  # powered (3.42 expected) and the gene is replicated in 20 control embryos.
  expect_lt(ep(out, "complete_depletion"), 0.05)

  st <- attr(out, "efdr_stats")
  expect_equal(st$gates_applied, "too_few_embryos,not_replicated,too_few_cells")
  expect_true(is.na(st$min_arm_cells_K))     # K is not located unless flat_cell_floor is on
  expect_equal(st$min_replicate_embryos, 2L)
  expect_equal(st$n_gate_not_replicated, 2L) # one_embryo and absent_both
  expect_equal(st$n_gate_too_few_cells,  2L) # trace_big_arm and absent_both
})

test_that("the previous policy is reproducible via `gates`", {
  f <- gate_fixture()
  legacy <- c("too_few_embryos", "not_detected_anywhere", "gain_from_zero",
              "gain_no_counts", "loss_from_trace", "flat_cell_floor")
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support, gates = legacy)
  # The two policies differ in BOTH directions, which is the point of keeping this path.
  # gain_from_zero killed real ectopic activation on an absolute zero-UMI rule -- this is the
  # shape of the osr1 sox9a call.
  expect_equal(ep(out, "gain_from_zero"), 1)
  # loss_from_trace reached the trace gene too, by a different route (1% control detection).
  expect_equal(ep(out, "trace_big_arm"), 1)
  # ...but the previous policy had no replication requirement, so a call whose every detection
  # sits in one embryo passed.
  expect_lt(ep(out, "one_embryo"), 0.05)

  st <- attr(out, "efdr_stats")
  expect_equal(st$gates_applied, paste(legacy, collapse = ","))
  # K is a single number for the entire table, located from the p75 of control detection --
  # here 38%, so ceil(3/0.38) = 8. It happens to land just under this cell type's 9-cell arm
  # and so fires on nothing; had detection been slightly lower it would have taken out every
  # gene in the cell type, well-expressed ones included. That sensitivity to a distributional
  # summary, rather than to the call being judged, is what the per-call test removes.
  expect_equal(st$min_arm_cells_K, 8L)
  expect_equal(st$p_ctrl_detection, 0.38)
})

test_that("gates are computed even when not applied, and unknown names error", {
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support,
                                      gates = character())
  # Nothing applied -> empirical_p is max(p_ashr, p_tail) for every row...
  expect_true(all(out$empirical_p < 0.05))
  # ...but membership is still counted, which is what makes the diagnostics comparable
  # across policies.
  expect_gt(attr(out, "efdr_stats")$n_gate_not_replicated, 0)

  expect_error(annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                            counts = f$counts, support = f$support,
                                            gates = c("too_few_cells", "no_such_gate")),
               "unknown gate")
})

test_that("missing counts disable the count-driven gates rather than gating everything", {
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      support = f$support)
  expect_true(all(out$empirical_p < 0.05))
})


# gate_reason: WHY a call was pinned, published alongside WHETHER it was.

gr <- function(out, id) out$gate_reason[match(id, out$id)]

test_that("gate_reason names the gate that pinned each call", {
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support)

  # Ungated calls carry no reason.
  expect_true(is.na(gr(out, "well_det_small_arm")))
  expect_true(is.na(gr(out, "complete_depletion")))

  # Gated calls name their gate.
  expect_equal(gr(out, "trace_big_arm"), "too_few_cells")
  expect_equal(gr(out, "one_embryo"), "not_replicated")

  # The reason is present exactly where the call was pinned, and nowhere else.
  expect_equal(is.na(out$gate_reason), out$empirical_p < 1)
})


test_that("gate_reason explains pins that expected_detections cannot", {
  # This is the case the column exists for. `one_embryo` clears the power gate
  # -- 3 expected detections, at the threshold -- so every published count
  # column says it was testable, yet it is pinned to 1 by not_replicated. Before
  # gate_reason it was indistinguishable from a gene that was tested and found
  # null, which is the DEG analogue of the A0 problem.
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support)

  row <- out[match("one_embryo", out$id), ]
  expect_gte(row$expected_detections, 3)      # passes the power gate
  expect_equal(row$empirical_p, 1)            # pinned anyway
  expect_equal(row$gate_reason, "not_replicated")
})


test_that("a gate that is off records no reason", {
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support,
                                      gates = character(0))

  expect_true(all(is.na(out$gate_reason)))
  expect_true(all(out$empirical_p < 1))
})


test_that("gate_reason accumulates when several gates fire on one call", {
  f <- gate_fixture()
  all_gates <- c("too_few_embryos", "not_replicated", "too_few_cells",
                 "not_detected_anywhere", "gain_from_zero", "gain_no_counts",
                 "loss_from_trace", "flat_cell_floor")
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support,
                                      gates = all_gates)

  # absent_both trips four of them; all are recorded, in evaluation order.
  reasons <- strsplit(gr(out, "absent_both"), ";", fixed = TRUE)[[1]]
  expect_true(all(c("not_replicated", "not_detected_anywhere",
                    "loss_from_trace", "too_few_cells") %in% reasons))
  expect_false(anyDuplicated(reasons) > 0)

  # and every named reason is a real gate name
  named <- unique(unlist(strsplit(na.omit(out$gate_reason), ";", fixed = TRUE)))
  expect_true(all(named %in% all_gates))
})


test_that("recording the reason does not change the verdict", {
  # The column is purely additive: it must not move empirical_p or empirical_fdr
  # under any gate configuration, including the retired full policy.
  f <- gate_fixture()
  all_gates <- c("too_few_embryos", "not_replicated", "too_few_cells",
                 "not_detected_anywhere", "gain_from_zero", "gain_no_counts",
                 "loss_from_trace", "flat_cell_floor")

  for (gset in list(.EFDR_DEFAULT_GATES, all_gates, character(0))) {
    out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                        counts = f$counts, support = f$support,
                                        gates = gset)
    # Pinned exactly when a reason was recorded, and pinned to exactly 1.
    expect_equal(out$empirical_p >= 1, !is.na(out$gate_reason))
    expect_true(all(out$empirical_p[!is.na(out$gate_reason)] == 1))
    expect_true(all(out$empirical_fdr[!is.na(out$gate_reason)] == 1))
  }
})


test_that("expected_detections is published, not just used internally", {
  # The design note assumed this was computed and discarded. It is not: it
  # survives on the returned table, and build_efdr_model.R writes that table
  # whole, so it already reaches the published DEG file.
  f <- gate_fixture()
  out <- annotate_empirical_fdr_model(stub_model(), f$deg, arm_cells = f$arm_cells,
                                      counts = f$counts, support = f$support)

  expect_true("expected_detections" %in% names(out))
  expect_equal(out$expected_detections, out$thin_arm_cells * out$max_det)
  # effect-independent: a function of arm size and detection rate only
  expect_false(any(is.na(out$expected_detections[is.finite(out$max_det)])))
})
