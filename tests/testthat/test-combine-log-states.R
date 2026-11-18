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
