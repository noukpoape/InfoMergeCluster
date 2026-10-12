# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi Noukpoape <noukpoape@gmail.com>

# =========================================================================== #
#  Test data generator (scenarios 1 to 4)                                     #
# =========================================================================== #

# frequencies by scenario and group (the last listed group applies to higher k)
.fsen_freq <- list(c(1.8, 2.2), c(1.5, 1.7, 1.9), c(1.5, 1.7, 1.9, 1.6),
                   c(1, 1.2, 1.4, 1.6, 1.8))

# same order of operations and random draws as the previous version:
# the generated data are identical for a given seed
fsenario <- function(x, senario, k, cst = 0) {
  s <- if (senario %in% 1:3) senario else 4L
  w <- .fsen_freq[[s]]
  w <- w[min(k, length(w))]
  switch(s,
         rnorm(1, -1/3, 1/3) + cst - cos(1.8*pi*x) - sin(w*x) + rnorm(1, 2, 0.4),
         rnorm(1, -0.5, 0.5) + cst - cos(w*pi*x^3) + rnorm(1, 2, 0.4),
         rnorm(1, -1/3, 1/3) + cst + sin(w*pi*x) + sin(pi*x^2) + rnorm(1, 2, 0.4),
         rnorm(1, -0.25, 0.25) + cst + cos(w*pi*x) - x^2 + rnorm(1, 2, 0.3))
}

# reference group sizes n.k of each scenario (used when n.k = NULL);
# the reference shifts are n.cts = 1.8 + 0.6 * (k - 1) for group k
.fsen_nk <- list(c(20, 30), c(20, 30, 25), c(20, 30, 25, 35), c(20, 30, 25, 35, 40))

#' Generate test data (scenarios 1 to 4)
#'
#' Generates curves from several groups following one of four scenarios,
#' and returns them both as a list of matrices and as a long data.frame
#' (the two formats contain exactly the same data).
#'
#' Each curve is observed at \code{n.min + Binomial(8, 0.5)} sorted uniform
#' times on [0, 1]; its group sets the mean shape and the constant shift.
#'
#' \strong{Reference settings of each scenario.} When \code{n.k} and
#' \code{n.cts} are left to \code{NULL}, the following default values are used:
#' \tabular{clll}{
#'   \strong{Scenario} \tab \strong{Groups} \tab \strong{Group sizes} \code{n.k} \tab \strong{Shifts} \code{n.cts} \cr
#'   1 \tab 2 \tab \code{c(20, 30)} (50 curves) \tab \code{c(1.8, 2.4)} \cr
#'   2 \tab 3 \tab \code{c(20, 30, 25)} (75 curves) \tab \code{c(1.8, 2.4, 3.0)} \cr
#'   3 \tab 4 \tab \code{c(20, 30, 25, 35)} (110 curves) \tab \code{c(1.8, 2.4, 3.0, 3.6)} \cr
#'   4 \tab 5 \tab \code{c(20, 30, 25, 35, 40)} (150 curves) \tab \code{c(1.8, 2.4, 3.0, 3.6, 4.2)} \cr
#' }
#' Other group sizes can be given in \code{n.k}; the default shifts are then
#' \code{1.8 + 0.6 * (k - 1)} for group \code{k}.
#'
#' @param senario scenario (1 to 4).
#' @param n.k number of curves in each group (\code{NULL} = reference sizes
#'   of the scenario, see Details).
#' @param n.cts constant shift of each group, one value per group
#'   (\code{NULL} = reference shifts, see Details).
#' @param n.min minimum number of measurements per curve.
#' @return A list with \code{liste} (list of matrices (x, y), one per curve,
#'   with a \code{"groupe"} attribute) and \code{df} (data.frame with columns
#'   \code{id}, \code{groupe}, \code{x}, \code{y}).
#' @examples
#' # scenario 1 with its reference settings: 2 groups of 20 and 30 curves
#' set.seed(1)
#' gen <- f.data.generation(senario = 1)
#' plot.data(gen$liste, legend = TRUE)
#' plot.data(gen$df, data.forma = 2)
#'
#' # reference settings of the other scenarios
#' table(attr(f.data.generation(senario = 2)$liste, "groupe"))
#' table(attr(f.data.generation(senario = 4)$liste, "groupe"))
#' @export
f.data.generation <- function(senario = 1, n.k = NULL, n.cts = NULL, n.min = 2) {
  
  # Output: list(liste = ..., df = ...)
  #   liste : list of matrices (x, y), one per curve ("groupe" attribute)
  #   df    : data.frame (id, groupe, x, y)
  
  # reference settings of the scenario when n.k / n.cts are not supplied
  s <- if (senario %in% 1:3) senario else 4L
  if (is.null(n.k)) n.k <- .fsen_nk[[s]]
  if (is.null(n.cts)) n.cts <- 1.8 + 0.6 * (seq_along(n.k) - 1)
  if (length(n.cts) < length(n.k))
    stop("'n.cts' must contain one value per group (length(n.k) = ", length(n.k), ").")
  
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
