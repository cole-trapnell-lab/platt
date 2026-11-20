test_that("one state is returned unchanged", {
  est <- matrix(c(-2, 0.5), ncol = 1, dimnames = list(c("a", "b"), "s1"))
  se <- matrix(c(0.3, 0.1), ncol = 1, dimnames = list(c("a", "b"), "s1"))
  out <- platt:::combine_log_states(est, se)
  expect_equal(unname(out$effect), c(-2, 0.5))
  expect_equal(unname(out$se), c(0.3, 0.1))
})

test_that("several states combine to the log of their total", {
  est <- log(matrix(c(1e-3, 2e-3, 5, 5), ncol = 2, byrow = TRUE, dimnames = list(c("a", "b"), c("p1", "p2"))))
  se <- matrix(c(0.2, 0.4, 0.1, 0.1), ncol = 2, byrow = TRUE)
  out <- platt:::combine_log_states(est, se)
  expect_equal(unname(out$effect), log(c(3e-3, 10)))
  w <- c(1, 2) / 3
  expect_equal(unname(out$se)[1], sqrt(sum(w^2 * c(0.2, 0.4)^2)))
  expect_equal(unname(out$se)[2], sqrt(2 * 0.5^2 * 0.1^2))
})

test_that("stable for very small and very large log levels", {
  est <- matrix(c(-800, -801, 800, 799), ncol = 2, byrow = TRUE)
  se <- matrix(0.1, 2, 2)
  out <- platt:::combine_log_states(est, se)
  expect_true(all(is.finite(out$effect)))
  expect_equal(out$effect[1], -800 + log(1 + exp(-1)), ignore_attr = TRUE)
})

test_that("infinite standard errors propagate only through states that carry weight", {
  est <- matrix(c(0, -Inf, 0, 0), ncol = 2, byrow = TRUE)
  se <- matrix(c(0.1, Inf, 0.1, Inf), ncol = 2, byrow = TRUE)
  out <- platt:::combine_log_states(est, se)
  expect_equal(unname(out$se)[1], 0.1)
  expect_true(is.infinite(out$se[2]))
})

test_that("a state unchanged from two equal parents reads as the total, not the product", {
  # rate r in the state and in both parents: the contrast is log(r / 2r) = -log 2 at every level,
  # where summing log estimates gave log(r / r^2), which depends on r
  r <- c(1e-3, 1e-2, 1e-1)
  pem <- log(cbind(s = r, p1 = r, p2 = r)); rownames(pem) <- paste0("g", 1:3)
  psem <- matrix(0.1, 3, 3, dimnames = dimnames(pem))
  res <- platt:::contrast_helper("s", c("p1", "p2"), PEM = pem, PSEM = psem)
  expect_equal(unname(res$raw_lfc), rep(-log(2), 3))
})

test_that("the total matches the direct computation", {
  set.seed(1)
  est <- matrix(rnorm(60, mean = -6, sd = 2), ncol = 3, dimnames = list(paste0("g", 1:20), c("p1", "p2", "p3")))
  se <- matrix(runif(60, 0.05, 0.5), ncol = 3, dimnames = dimnames(est))
  out <- platt:::combine_log_states(est, se)
  expect_equal(out$effect, log(rowSums(exp(est))))
  w <- exp(est) / rowSums(exp(est))
  expect_equal(out$se, sqrt(rowSums(w^2 * se^2)))
  expect_identical(names(out$se), rownames(est))
})

test_that("standard errors are matched to estimates by row name", {
  est <- log(matrix(c(1e-3, 2e-3, 5, 1), ncol = 2, byrow = TRUE, dimnames = list(c("a", "b"), c("p1", "p2"))))
  se <- matrix(c(0.2, 0.4, 0.1, 0.3), ncol = 2, byrow = TRUE, dimnames = dimnames(est))
  expected <- platt:::combine_log_states(est, se)
  expect_equal(platt:::combine_log_states(est, se[c("b", "a"), ]), expected)
  # one state: the standard error must follow the row name too
  expect_equal(platt:::combine_log_states(est[, 1, drop = FALSE], se[c("b", "a"), 1, drop = FALSE])$se, c(a = 0.2, b = 0.1))
  se_other <- se; rownames(se_other) <- c("a", "c")
  expect_error(platt:::combine_log_states(est, se_other), "same row names")
  expect_error(platt:::combine_log_states(est, se[, 1, drop = FALSE]), "same dimensions")
})

test_that("a row with no finite total still has a defined standard error", {
  est <- matrix(c(-Inf, -Inf, 0, -Inf), ncol = 2, byrow = TRUE)
  se <- matrix(c(0.3, 0.4, 0.1, 0.2), ncol = 2, byrow = TRUE)
  out <- platt:::combine_log_states(est, se)
  expect_equal(unname(out$effect), c(-Inf, 0))
  expect_equal(unname(out$se), c(sqrt(0.5^2 * 0.3^2 + 0.5^2 * 0.4^2), 0.1))
})

test_that("with recentring, a gene unchanged from its parents reads zero for any number of parents", {
  # the configuration classify_gene_patterns uses: cv_threshold, n and an effect mode interval
  set.seed(2)
  n_genes <- 2000
  n_up <- 100
  for (k in c(2, 3, 5)) {
    parents <- paste0("p", seq_len(k))
    level <- runif(n_genes, log(1e-4), log(1e-1))
    truth <- matrix(level, n_genes, k + 1, dimnames = list(paste0("g", seq_len(n_genes)), c("s", parents)))
    truth[seq_len(n_up), "s"] <- truth[seq_len(n_up), "s"] + 3 # genes that really are higher in the state
    psem <- matrix(0.05, n_genes, k + 1, dimnames = dimnames(truth))
    pem <- truth + rnorm(length(truth), sd = 0.05)
    res <- platt:::contrast_helper("s", parents,
      PEM = pem, PSEM = psem, n = 500,
      ash.control = list(mode = c(-10, 10)), cv_threshold = 100
    )
    unchanged <- res[-seq_len(n_up), ]
    expect_equal(res$coefficient_mode[1], -log(k), tolerance = 0.05)
    expect_lt(abs(median(unchanged$raw_lfc)), 0.02)
    expect_true(all(abs(unchanged$shrunken_lfc) < 1))
    expect_true(all(res$shrunken_lfc[seq_len(n_up)] > 1))
    # without recentring the same genes sit at -log(k), past the default threshold of 1 from 3 parents up
    raw <- platt:::contrast_helper("s", parents, PEM = pem, PSEM = psem, n = 500, cv_threshold = 100)
    expect_equal(median(raw$raw_lfc[-seq_len(n_up)]), -log(k), tolerance = 0.05)
  }
})

test_that("a state with no counts (negligible weight, infinite SE) adds nothing to the total's SE", {
  # the closed-form fit gives a state with no counts log(abs_expr_thresh) - 20 and an infinite SE
  sentinel <- log(1e-3) - 20
  est <- rbind(g1 = c(log(0.05), sentinel), g2 = c(log(0.05), log(0.02)), g3 = c(sentinel, sentinel))
  se <- rbind(g1 = c(0.2, Inf), g2 = c(0.2, Inf), g3 = c(Inf, Inf))
  out <- platt:::combine_log_states(est, se)
  expect_equal(unname(out$se["g1"]), 0.2, tolerance = 1e-8)       # only the counted state contributes
  expect_equal(unname(out$effect["g1"]), log(0.05), tolerance = 1e-6)
  expect_true(is.infinite(out$se["g2"]))                            # unknown level with real weight: still infinite
  expect_true(is.infinite(out$se["g3"]))                            # no counts anywhere: nothing is known
  # the tolerance is an argument; 0 gives the old rule (any non-zero weight propagates)
  expect_true(is.infinite(platt:::combine_log_states(est, se, negligible_weight = 0)$se["g1"]))
})
