# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi NOUKPOAPE <noukpoape@gmail.com>

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
