set.seed(42)
gen <- f.data.generation(senario = 1, n.k = c(60, 60))
long <- gen$df

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
  lst <- split(long[, c("x", "y")], long$id)
  fit <- infomergecluster(lst, K = 2, verbose = FALSE, seed = 3)
  expect_equal(fit$K, 2L)
  expect_setequal(names(fit$cluster), names(lst))
})

test_that("adjusted_rand is 1 for identical partitions", {
  a <- rep(1:3, each = 5)
  expect_equal(adjusted_rand(a, a), 1)
  expect_equal(adjusted_rand(a, c(3, 1, 2)[a]), 1)
})

test_that("scenario 1 with its reference settings runs with both data formats", {
  set.seed(1)
  g1 <- f.data.generation(senario = 1, n.k = c(43, 57), n.cts = c(1.8, 2.4), n.min = 2)
  f_df <- infomergecluster(g1$df, K = 2, re = "intercept", verbose = FALSE, seed = 1)
  f_ls <- infomergecluster(g1$liste, K = 2, re = "intercept", verbose = FALSE, seed = 1)
  expect_equal(f_df$K, 2L)
  expect_equal(f_ls$K, 2L)
  expect_length(f_df$cluster, 100)
})

test_that("f.data.generation reproduces the reference group sizes", {
  ref <- list(
    list(n.k = c(43, 57),             n.cts = c(1.8, 2.4),           n.min = 2),
    list(n.k = c(35, 43, 72),         n.cts = c(2, 2.3, 2.7),        n.min = 3),
    list(n.k = c(34, 67, 71, 28),     n.cts = c(0.9, 1.5, 2.2, 2.4), n.min = 3),
    list(n.k = c(60, 85, 45, 70, 50), n.cts = rep(0, 5),             n.min = 3)
  )
  for (s in 1:4) {
    g <- do.call(f.data.generation, c(list(senario = s), ref[[s]]))
    expect_equal(as.vector(table(attr(g$liste, "groupe"))), ref[[s]]$n.k)
    expect_false(anyNA(g$df$y))
    expect_gte(min(sapply(g$liste, nrow)), ref[[s]]$n.min)
  }
})

test_that("covariates and covariate-dependent weights", {
  set.seed(2)
  dat <- long[, c("id", "x", "y")]
  n_i <- as.vector(table(factor(dat$id, levels = unique(dat$id))))
  dat$sex <- rep(sample(c("F", "M"), length(n_i), replace = TRUE), n_i)
  dat$age <- rep(round(runif(length(n_i), 20, 60)), n_i)
  fit <- infomergecluster(dat, K = 2, n_cov = 2, verbose = FALSE, seed = 1)
  expect_equal(dim(fit$gamma), c(2L, 2L))
  expect_null(fit$alpha)
  pr <- predict(fit, times = c(0, 0.5, 1), covariates = c(1, 40))
  expect_equal(nrow(pr), 3)
  fit_w <- infomergecluster(dat, K = 2, n_cov = 2, weights = "covariates",
                            verbose = FALSE, seed = 1)
  expect_equal(dim(fit_w$alpha), c(3L, 2L))
  expect_error(infomergecluster(long, K = 2, weights = "covariates", verbose = FALSE))
})

test_that("infomergecluster_citation returns the BibTeX entry", {
  ref <- infomergecluster_citation()
  expect_match(ref, "^@manual\\{noukpoape2026infomergecluster,")
  expect_match(ref, "author=\\{Noukpoape, Komlan Midodzi and Preda, Cristian\\}")
  expect_match(infomergecluster_citation("text"), "^Noukpoape, K. M. and Preda, C. \\(2026\\)")
  bib <- tempfile(fileext = ".bib")
  infomergecluster_citation(file = bib)
  expect_true(file.exists(bib))
})
