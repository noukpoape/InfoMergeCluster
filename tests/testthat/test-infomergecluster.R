test_that("infomergecluster recovers two well-separated groups", {
  set.seed(1)
  g <- f.data.generation(senario = 1, n.k = c(40, 50), n.cts = c(0, 2.5), n.min = 3)
  fit <- infomergecluster(g$liste, K_max = 4, seed = 1, verbose = FALSE)
  expect_s3_class(fit, "infomergecluster")
  expect_false(fit$K_fixed)
  expect_equal(fit$K, 2L)
  expect_gt(adjusted_rand(fit$cluster, attr(g$liste, "groupe")), 0.8)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-8)
  expect_true(fit$sigma2 > 0)

  pr <- predict(fit, times = c(0.1, 0.5, 0.9), interval = "prediction")
  expect_named(pr, c("time", "C1", "C2", "C1_lwr", "C1_upr", "C2_lwr", "C2_upr"),
               ignore.order = TRUE)
  expect_true(all(pr$C1_lwr < pr$C1 & pr$C1 < pr$C1_upr))

  fc <- fitted_curves(fit, times = c(0.2, 0.8))
  expect_equal(dim(fc), c(90L, 2L))

  expect_output(print(fit), "Selected K")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fit, envelope = TRUE))
})

test_that("a user-supplied K is used as is, without selection", {
  set.seed(3)
  g <- f.data.generation(senario = 2, n.k = c(30, 30, 30), n.cts = c(0, 1, 2), n.min = 3)
  fit <- infomergecluster(g$df, K = 2, seed = 1, verbose = FALSE)
  expect_equal(fit$K, 2L)
  expect_true(fit$K_fixed)
  expect_equal(ncol(fit$beta), 2L)
  expect_null(fit$outer)
  expect_output(print(fit), "K imposed")

  fit1 <- infomergecluster(g$df, K = 1, verbose = FALSE)
  expect_equal(fit1$K, 1L)
  expect_true(all(fit1$cluster == 1L))

  expect_error(infomergecluster(g$df, K = 0, verbose = FALSE), "integer")
  expect_error(infomergecluster(g$df, K = 60, verbose = FALSE), "too large")
})

test_that("covariates and covariate-dependent weights", {
  set.seed(4)
  g <- f.data.generation(senario = 1, n.k = c(40, 40), n.cts = c(0, 2.5), n.min = 3)
  d <- g$df[, c("id", "x", "y")]
  grp <- attr(g$liste, "groupe")
  age <- stats::rnorm(80, 50, 10) + 10 * (grp == 2)
  d$age <- age[d$id]
  d$y <- d$y + 0.05 * (d$age - 50)

  fit <- infomergecluster(d, K = 2, n_cov = 1, seed = 1, verbose = FALSE)
  expect_equal(dim(fit$gamma), c(1L, 2L))
  expect_true(all(abs(fit$gamma - 0.05) < 0.05))
  expect_null(fit$alpha)
  p1 <- predict(fit, times = 0.5, covariates = 60)
  p0 <- predict(fit, times = 0.5)
  expect_equal(unname(unlist(p1[, -1] - p0[, -1])),
               unname((60 - fit$covariates$center) * fit$gamma[1, ]), tolerance = 1e-10)
  expect_error(predict(fit, covariates = c(1, 2)), "value")

  fitw <- infomergecluster(d, K = 2, n_cov = 1, weights = "covariates", seed = 1, verbose = FALSE)
  expect_equal(dim(fitw$alpha), c(2L, 2L))
  expect_equal(unname(fitw$alpha[, 1]), c(0, 0))
  expect_output(print(fitw), "multinomial logistic")

  expect_error(infomergecluster(d, weights = "covariates", verbose = FALSE), "n_cov")
})
