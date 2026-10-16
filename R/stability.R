# reruns infomergecluster() with several seeds; returns the fit with the best
# criterion, supplemented with a stability diagnostic (frequency of K, ARI between runs).
# The arguments of infomergecluster() go through '...' (do not include 'seed').

#' Stability of the partition across random seeds
#'
#' Runs \code{\link{infomergecluster}} with several seeds and measures the
#' agreement between the resulting partitions (adjusted Rand index).
#'
#' @param data data (see \code{\link{infomergecluster}}).
#' @param ... further arguments of \code{\link{infomergecluster}} (do not
#'   include \code{seed}).
#' @param seeds seeds used.
#' @param verbose print progress messages.
#' @return The fit with the best criterion, with an extra element
#'   \code{stability}: seeds, K obtained and their frequencies, criteria,
#'   matrix of pairwise ARIs between runs and mean ARI.
#' @examples
#' \donttest{
#' set.seed(1)
#' g <- f.data.generation(senario = 1, n.k = c(30, 40), n.cts = c(0, 2))
#' st <- infomergecluster_stability(g$liste, K_max = 3, seeds = 1:3)
#' st$stability$K_freq
#' }
#' @export

infomergecluster_stability <- function(data, ..., seeds = 1:10, verbose = FALSE) {
  fits <- lapply(seeds, function(s) {
    if (verbose) message(sprintf("seed %s", s))
    infomergecluster(data, ..., seed = s, verbose = FALSE)
  })
  Ks <- vapply(fits, function(f) as.integer(f$K), integer(1))
  crit <- vapply(fits, function(f) {
    row <- f$table[f$table$K == f$K, ]
    if (f$criterion == "BIC") row$BIC else row$ICL
  }, numeric(1))
  n <- length(fits)
  ari <- diag(1, n)
  if (n > 1) {
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      a <- fits[[i]]$cluster
      ari[i, j] <- ari[j, i] <- adjusted_rand(a, fits[[j]]$cluster[names(a)])
    }
  }
  best <- fits[[which.min(crit)]]
  best$stability <- list(seeds = seeds, K = Ks, K_freq = table(Ks), criterion = crit,
                         ari = ari, mean_ari = if (n > 1) mean(ari[upper.tri(ari)]) else 1)
  best
}
