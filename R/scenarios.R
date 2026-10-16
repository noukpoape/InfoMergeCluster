#' Generate test data (scenarios 1 to 4)
#'
#' \code{f.data.generation} simulates curves split into groups following one
#' of four scenarios; \code{fsenario} evaluates one curve of scenario
#' \code{senario} for group \code{k}.
#'
#' A curve of group \code{k} is
#' \eqn{U(-a, a) + cst + shape(x, w_k) + N(2, sd)}, with one draw of each
#' random shift per curve; \eqn{a}, \eqn{sd}, the shape and the group
#' frequencies \eqn{w_k} depend on the scenario. Each curve has
#' \code{n.min + Binomial(8, 0.5)} measurements, at sorted uniform times on
#' [0, 1].
#'
#' @param x observation times.
#' @param senario scenario number (1 to 4; any other value gives
#'   scenario 4).
#' @param k group number.
#' @param cst constant added to the group curve.
#' @param n.k group sizes.
#' @param n.cts group constants (same length as \code{n.k}).
#' @param n.min minimum number of measurements per curve.
#' @return \code{fsenario}: vector of curve values.
#'   \code{f.data.generation}: list of two elements holding the same data,
#'   \item{liste}{list of \code{(x, y)} matrices, one per curve, with
#'     attribute \code{"groupe"};}
#'   \item{df}{long data frame with columns \code{id}, \code{groupe},
#'     \code{x}, \code{y}.}
#' @examples
#' set.seed(2)
#' g <- f.data.generation(senario = 2, n.k = c(20, 30, 25), n.cts = c(0, 1, 2))
#' str(g$df)
#' plot_data(g$liste, legend = TRUE)
#' @name f.data.generation
NULL

# =========================================================================== #
#  Test data generator (scenarios 1 to 4)                                     #
# =========================================================================== #

# parameters by scenario: uniform shift U(-a, a), standard deviation sd of the
# normal shift N(2, sd), frequencies w by group (the last one applies to higher k)
.fsen_par <- list(
  list(a = 1/3,  sd = 0.4, w = c(1.8, 2.2)),
  list(a = 0.5,  sd = 0.4, w = c(1.5, 1.7, 1.9)),
  list(a = 1/3,  sd = 0.4, w = c(1.5, 1.7, 1.9, 1.6)),
  list(a = 0.25, sd = 0.3, w = c(1, 1.2, 1.4, 1.6, 1.8))
)

# mean shape of the curves, by scenario (w: frequency of the group)
.fsen_shape <- list(
  function(x, w) -cos(1.8*pi*x) - sin(w*x),
  function(x, w) -cos(w*pi*x^3),
  function(x, w)  sin(w*pi*x) + sin(pi*x^2),
  function(x, w)  cos(w*pi*x) - x^2
)

# curve of group k: U(-a, a) + cst + shape + N(2, sd), one draw of each shift
# per curve; scenario 4 for any value of 'senario' outside 1:3.
# Same order of random draws (runif then rnorm) as the version in
# CodeLatex_InfoMerge.txt: identical data for a given seed

#' @rdname f.data.generation
#' @export

fsenario <- function(x, senario, k, cst = 0) {
  s <- if (senario %in% 1:3) senario else 4L
  p <- .fsen_par[[s]]
  w <- p$w[min(k, length(p$w))]
  stats::runif(1, -p$a, p$a) + cst + .fsen_shape[[s]](x, w) + stats::rnorm(1, 2, p$sd)
}

#' @rdname f.data.generation
#' @export

f.data.generation <- function(senario = 1, n.k = c(20, 30), n.cts = c(1.8, 2.4), n.min = 2) {
  
  # Output: list(liste = ..., df = ...)
  #   liste : list of matrices (x, y), one per curve ("groupe" attribute)
  #   df    : data.frame (id, groupe, x, y)
  
  K <- length(n.k)
  
  courbes.liste <- vector("list", sum(n.k))  # list format
  courbes.df    <- vector("list", sum(n.k))  # data.frame format
  id <- 0
  
  for (k in seq_len(K)) {
    for (i in seq_len(n.k[k])) {
      id <- id + 1
      
      # A single simulation per curve, used for both formats:
      # both outputs therefore contain exactly the same data
      x <- sort(stats::runif(n.min + stats::rbinom(1, 8, 0.5)))
      y <- fsenario(x, senario, k, n.cts[k])
      
      courbes.liste[[id]] <- cbind(x = x, y = y)
      courbes.df[[id]]    <- data.frame(id = id, groupe = k, x = x, y = y)
    }
  }
  
  # List format, with the group of each curve as an attribute
  df.liste <- courbes.liste
  attr(df.liste, "groupe") <- rep(seq_len(K), times = n.k)
  
  # data.frame format
  df.frame <- do.call(rbind, courbes.df)
  
  return(list(liste = df.liste, df = df.frame))
}
