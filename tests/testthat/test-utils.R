test_that("adjusted_rand is 1 for identical partitions", {
  expect_equal(adjusted_rand(c(1, 1, 2, 2, 3), c("a", "a", "b", "b", "c")), 1)
})

test_that("list and data-frame inputs give the same long format", {
  set.seed(1)
  g <- f.data.generation(senario = 1, n.k = c(5, 6))
  a <- InfoMergeCluster:::.im_as_long(g$df)
  b <- InfoMergeCluster:::.im_as_long(g$liste)
  expect_equal(nrow(a), nrow(g$df))
  expect_equal(a$t, b$t)
  expect_equal(a$y, b$y)
  expect_equal(length(unique(b$id)), 11)
})

test_that("covariate columns are read by position", {
  d <- data.frame(id = rep(1:3, each = 2), time = rep(c(0, 1), 3),
                  y = 1:6, age = rep(c(30, 40, 50), each = 2), sex = rep(c("F", "M", "F"), each = 2))
  a <- InfoMergeCluster:::.im_as_long(d, n_cov = 2)
  expect_identical(names(a), c("id", "t", "y", "z1", "z2"))
  expect_equal(a$y, 1:6)
  X <- InfoMergeCluster:::.im_covariate_design(a[c("z1", "z2")], rep(TRUE, 6))
  expect_equal(ncol(X$X), 2L)
  expect_equal(unname(colMeans(X$X)), c(0, 0), tolerance = 1e-12)
})

test_that("the penalty matrix vanishes on polynomials of degree < m", {
  basis <- InfoMergeCluster:::.im_make_basis(seq(0, 10, length.out = 50), df = 10, degree = 3)
  P <- InfoMergeCluster:::.im_penalty_matrix(basis, 2)
  # B-spline coefficients of a straight line: Greville abscissae
  k <- basis$knots
  grev <- vapply(seq_len(basis$p), function(j) mean(k[(j + 1):(j + 3)]), numeric(1))
  beta <- 2 + 0.5 * grev
  expect_lt(abs(sum(beta * (P %*% beta))), 1e-8)
  expect_equal(qr(P)$rank, basis$p - 2)
})

test_that("f.data.generation returns the same curves in both formats", {
  set.seed(5)
  g <- f.data.generation(senario = 3, n.k = c(4, 5, 6), n.cts = c(0, 1, 2))
  expect_equal(length(g$liste), 15L)
  expect_equal(attr(g$liste, "groupe"), rep(1:3, c(4, 5, 6)))
  expect_equal(unname(do.call(rbind, g$liste)[, "y"]), g$df$y)
})
