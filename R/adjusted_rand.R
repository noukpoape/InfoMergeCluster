# adjusted Rand index (no external dependency)

#' Adjusted Rand index
#'
#' Computes the adjusted Rand index (Hubert and Arabie, 1985) between two
#' partitions, without external dependencies.
#'
#' @param a,b two label vectors of the same length.
#' @return A number: 1 for identical partitions, about 0 for independent
#'   partitions.
#' @examples
#' adjusted_rand(c(1, 1, 2, 2), c("a", "a", "b", "b"))
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
