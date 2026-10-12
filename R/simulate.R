# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi NOUKPOAPE <noukpoape@gmail.com>

#' Simulate sparse longitudinal data
#'
#' Simulates subjects from K distinct mean trajectories on [0, 1], with a
#' random intercept, a random slope and measurement noise.
#'
#' @param N number of subjects.
#' @param K number of clusters (at most 6).
#' @param n_range range of the number of measurements per subject.
#' @param sd_int,sd_slope standard deviations of the random intercept and slope.
#' @param sd_eps standard deviation of the measurement noise.
#' @param props cluster proportions (default: equal).
#' @param seed optional random seed.
#' @return A long data.frame with columns \code{id}, \code{time}, \code{y}
#'   and \code{true_class}.
#' @export
simulate_sparse_longitudinal <- function(N = 1000, K = 3, n_range = c(2, 5),
                                         sd_int = 0.5, sd_slope = 0.5, sd_eps = 0.3,
                                         props = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  curves <- list(
    function(t) 2 * sin(2 * pi * t),
    function(t) 3 * t - 1.5,
    function(t) -2 * t^2 + 1,
    function(t) 1.5 * cos(3 * pi * t),
    function(t) 2 * exp(-20 * (t - 0.5)^2) - 1,
    function(t) -3 * t + 1.5
  )
  if (K > length(curves)) stop("K <= 6 for this simulator.")
  if (is.null(props)) props <- rep(1 / K, K)
  z <- sample.int(K, N, replace = TRUE, prob = props)
  n <- sample(n_range[1]:n_range[2], N, replace = TRUE)
  id <- rep(seq_len(N), n)
  t <- stats::runif(sum(n))
  zz <- z[id]
  y <- numeric(length(t))
  for (k in seq_len(K)) y[zz == k] <- curves[[k]](t[zz == k])
  y <- y + stats::rnorm(N, 0, sd_int)[id] +
    stats::rnorm(N, 0, sd_slope)[id] * (t - 0.5) + stats::rnorm(length(t), 0, sd_eps)
  data.frame(id = id, time = t, y = y, true_class = zz)
}

#' Adjusted Rand index
#'
#' Adjusted Rand index between two partitions (no external dependency).
#'
#' @param a,b two vectors of cluster labels of the same length.
#' @return The adjusted Rand index (1 = identical partitions).
#' @export
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  cmb <- function(x) sum(choose(x, 2))
  s_ij <- cmb(tab)
  s_a <- cmb(rowSums(tab))
  s_b <- cmb(colSums(tab))
  e <- s_a * s_b / choose(sum(tab), 2)
  (s_ij - e) / ((s_a + s_b) / 2 - e)
}
