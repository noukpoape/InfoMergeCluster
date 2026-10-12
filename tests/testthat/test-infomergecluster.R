dat <- simulate_sparse_longitudinal(N = 120, K = 2, seed = 42)
long <- dat[, c("id", "time", "y")]

test_that("imposed K returns exactly K clusters", {
  fit <- infomergecluster(long, K = 2, verbose = FALSE, seed = 1)
  expect_s3_class(fit, "infomergecluster")
  expect_equal(fit$K, 2L)
  expect_true(fit$K_fixed)
  expect_length(fit$cluster, 120)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-8)
  expect_null(fit$outer)
})

test_that("automatic selection returns a valid fit", {
  fit <- infomergecluster(long, K_max = 3, verbose = FALSE, seed = 1)
  expect_false(fit$K_fixed)
  expect_true(fit$K >= 1 && fit$K <= 3)
  expect_true(all(c("K", "BIC", "ICL") %in% names(fit$table)))
})

test_that("invalid K is rejected", {
  expect_error(infomergecluster(long, K = 0, verbose = FALSE))
  expect_error(infomergecluster(long, K = 2.5, verbose = FALSE))
  expect_error(infomergecluster(long, K = 1000, verbose = FALSE))
})

test_that("predict and fitted_curves have the expected shape", {
  fit <- infomergecluster(long, K = 2, verbose = FALSE, seed = 1)
  pr <- predict(fit, times = seq(0, 1, length.out = 11), interval = "confidence")
  expect_equal(nrow(pr), 11)
  expect_true(all(c("C1", "C1_lwr", "C1_upr") %in% names(pr)))
  fc <- fitted_curves(fit, times = seq(0, 1, length.out = 5))
  expect_equal(dim(fc), c(120L, 5L))
})

test_that("list input is accepted", {
  lst <- split(long[, c("time", "y")], long$id)
  fit <- infomergecluster(lst, K = 2, verbose = FALSE, seed = 3)
  expect_equal(fit$K, 2L)
  expect_setequal(names(fit$cluster), names(lst))
})

test_that("adjusted_rand is 1 for identical partitions", {
  a <- rep(1:3, each = 5)
  expect_equal(adjusted_rand(a, a), 1)
  expect_equal(adjusted_rand(a, c(3, 1, 2)[a]), 1)
})
