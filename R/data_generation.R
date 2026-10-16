# Part of the InfoMergeCluster package.
# Authors: Komlan Midodzi Noukpoape <noukpoape@gmail.com>, Cristian Preda

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
# Random draws in the order runif then rnorm, as in the reference
# implementation of the scenarios.
fsenario <- function(x, senario, k, cst = 0) {
  s <- if (senario %in% 1:3) senario else 4L
  p <- .fsen_par[[s]]
  w <- p$w[min(k, length(p$w))]
  runif(1, -p$a, p$a) + cst + .fsen_shape[[s]](x, w) + rnorm(1, 2, p$sd)
}

#' Generate test data (scenarios 1 to 4)
#'
#' Generates curves from several groups following one of the four scenarios
#' of Zhou, Zhang and Tu (2023), and returns them both as a list of matrices
#' and as a long data.frame (the two formats contain exactly the same data).
#'
#' Each curve is observed at \code{n.min + Binomial(8, 1/2)} sorted uniform
#' times on [0, 1]. A curve of group \eqn{k} is
#' \eqn{y_{ij} = f_{0k}(t_{ij}) + u_i}, where
#' \eqn{u_i = U_i + E_i}, \eqn{U_i \sim U(-a, a)} and
#' \eqn{E_i \sim N(0, s^2)}, is a random intercept (one draw of each per
#' curve) and no measurement error is added. The mean curves and the
#' parameters \eqn{(a, s)} are
#' \tabular{ll}{
#'   S1: \tab \eqn{f_{0k}(t) = 2 + c_k - \cos(1.8\pi t) - \sin(\omega_k t)}, \eqn{a = 1/3}, \eqn{s = 0.4}, \cr
#'   S2: \tab \eqn{f_{0k}(t) = 2 + c_k - \cos(\omega_k \pi t^3)}, \eqn{a = 1/2}, \eqn{s = 0.4}, \cr
#'   S3: \tab \eqn{f_{0k}(t) = 2 + c_k + \sin(\omega_k \pi t) + \sin(\pi t^2)}, \eqn{a = 1/3}, \eqn{s = 0.4}, \cr
#'   S4: \tab \eqn{f_{0k}(t) = 2 + c_k + \cos(\omega_k \pi t) - t^2}, \eqn{a = 1/4}, \eqn{s = 0.3}, \cr
#' }
#' with frequencies \eqn{\omega = (1.8, 2.2)} in S1, \eqn{(1.5, 1.7, 1.9)} in
#' S2, \eqn{(1.5, 1.7, 1.9, 1.6)} in S3 and \eqn{(1, 1.2, 1.4, 1.6, 1.8)} in
#' S4 (the last value is used for any further group), and level shifts
#' \eqn{c_k} given by \code{n.cts}.
#'
#' \strong{Reference settings of each scenario.} The reference data sets of
#' the four scenarios (study A of the simulation study of InfoMergeCluster)
#' use the following group sizes \code{n.k}, shifts \code{n.cts} and
#' \code{n.min}:
#' \tabular{cccccc}{
#'   \strong{Scenario} \tab \strong{K} \tab \strong{N} \tab \strong{Group sizes} \code{n.k} \tab \strong{Shifts} \code{n.cts} \tab \code{n.min} \cr
#'   S1 \tab 2 \tab 100 \tab \strong{\code{c(43, 57)}} \tab \code{c(1.8, 2.4)} \tab 2 \cr
#'   S2 \tab 3 \tab 150 \tab \strong{\code{c(35, 43, 72)}} \tab \code{c(2, 2.3, 2.7)} \tab 3 \cr
#'   S3 \tab 4 \tab 200 \tab \strong{\code{c(34, 67, 71, 28)}} \tab \code{c(0.9, 1.5, 2.2, 2.4)} \tab 3 \cr
#'   S4 \tab 5 \tab 310 \tab \strong{\code{c(60, 85, 45, 70, 50)}} \tab \code{rep(0, 5)} \tab 3 \cr
#' }
#' The default values of the arguments (\code{n.k = c(20, 30)},
#' \code{n.cts = c(1.8, 2.4)}, \code{n.min = 2}) give a small two-group data
#' set; give the values above to reproduce the reference settings.
#'
#' @references
#' Zhou, J., Zhang, Y. and Tu, W. (2023). clusterMLD: An efficient hierarchical clustering method for multivariate longitudinal data. Journal of Computational and Graphical Statistics, 32(3), 1131--1144. \doi{10.1080/10618600.2022.2149540}
#'
#' @param senario scenario (1 to 4).
#' @param n.k number of curves in each group (see the reference settings in Details).
#' @param n.cts constant level shift of each group, one value per group.
#' @param n.min minimum number of measurements per curve.
#' @return A list with \code{liste} (list of matrices (x, y), one per curve,
#'   with a \code{"groupe"} attribute) and \code{df} (data.frame with columns
#'   \code{id}, \code{groupe}, \code{x}, \code{y}).
#' @examples
#' # scenario 1 with its reference settings: 2 groups of 43 and 57 curves
#' set.seed(1)
#' gen <- f.data.generation(senario = 1, n.k = c(43, 57), n.cts = c(1.8, 2.4), n.min = 2)
#' plot.data(gen$liste, legend = TRUE)
#' plot.data(gen$df, data.forma = 2)
#' table(attr(gen$liste, "groupe"))
#'
#' # reference settings of the other scenarios
#' s2 <- f.data.generation(2, n.k = c(35, 43, 72), n.cts = c(2, 2.3, 2.7), n.min = 3)
#' s3 <- f.data.generation(3, n.k = c(34, 67, 71, 28), n.cts = c(0.9, 1.5, 2.2, 2.4), n.min = 3)
#' s4 <- f.data.generation(4, n.k = c(60, 85, 45, 70, 50), n.cts = rep(0, 5), n.min = 3)
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
      x <- sort(runif(n.min + rbinom(1, 8, 0.5)))
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

#' Plot the trajectories with their true groups
#'
#' Draws the curves generated by \code{\link{f.data.generation}}, one
#' color per true group.
#'
#' @param data output of \code{\link{f.data.generation}}: either the list of
#'   matrices (\code{data.forma = 1}) or the data.frame (\code{data.forma = 2}).
#' @param n.k group sizes (only needed for a list without a \code{"groupe"} attribute).
#' @param data.forma 1 for a list of matrices, 2 for a data.frame.
#' @param xlab,ylab,main axis labels and title.
#' @param col colors, one per group.
#' @param legend \code{TRUE} to draw a legend of the groups (default \code{FALSE}).
#' @param ... further arguments passed to \code{plot}.
#' @return \code{NULL}, invisibly.
#' @export
plot.data <- function(data, n.k = NULL, data.forma = 1,
                      xlab = "Time (days)", ylab = "Weight",
                      main = "Trajectories (true labels)",
                      col = NULL, legend = FALSE, ...) {
  
  # ---- 1. Extraction of the coordinates and groups -------------------------
  if (data.forma == 1) {
    # data: list of matrices (x, y), one per curve
    groupe <- attr(data, "groupe")
    if (is.null(groupe)) {
      if (is.null(n.k)) {
        stop("Supply 'n.k' (group sizes) or a 'groupe' attribute on 'data'.")
      }
      if (sum(n.k) != length(data)) {
        stop("sum(n.k) must equal the number of curves in 'data'.")
      }
      groupe <- rep(seq_along(n.k), times = n.k)
    }
    t_all <- unlist(lapply(data, function(d) d[, 1]))
    y_all <- unlist(lapply(data, function(d) d[, 2]))
    
  } else {
    # data: data.frame with columns id, groupe, x, y
    if (!"groupe" %in% names(data)) {
      stop("The data.frame must contain a 'groupe' column.")
    }
    t_all <- data$x
    y_all <- data$y
    ids   <- unique(data$id)
    groupe <- data$groupe[match(ids, data$id)]  # one group per curve
  }
  
  # ---- 2. Colors: one per group, whatever the number of groups --------------
  niveaux <- sort(unique(groupe))
  K <- length(niveaux)
  if (is.null(col)) {
    col <- if (K <= 8) seq_len(K) else hcl.colors(K, "Dark 3")
  }
  col_courbe <- col[match(groupe, niveaux)]
  
  # ---- 3. Plot ---------------------------------------------------------------
  plot(NA, xlim = range(t_all), ylim = range(y_all),
       xlab = xlab, ylab = ylab, main = main, ...)
  
  if (data.forma == 1) {
    for (i in seq_along(data)) {
      lines(data[[i]], col = col_courbe[i], lty = 1, type = "b")
    }
  } else {
    rows <- split(seq_len(nrow(data)), factor(data$id, levels = ids))
    for (i in seq_along(ids)) {
      r <- rows[[i]]
      lines(data$x[r], data$y[r], col = col_courbe[i], lty = 1, type = "b")
    }
  }
  
  if (legend) {
    legend("topright", legend = paste("Group", niveaux),
           col = col, lty = 1, pch = 1, bty = "n")
  }
  
  
  invisible(NULL)
}
