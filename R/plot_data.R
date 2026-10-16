#' Plot trajectories with their true groups
#'
#' Plots the curves produced by \code{\link{f.data.generation}}, coloured
#' by group.
#'
#' @param data list of \code{(x, y)} matrices (\code{data.forma = 1}) or
#'   data frame with columns \code{id}, \code{groupe}, \code{x}, \code{y}
#'   (\code{data.forma = 2}).
#' @param n.k group sizes; only needed when the list has no
#'   \code{"groupe"} attribute.
#' @param data.forma 1 (list) or 2 (data frame).
#' @param xlab,ylab,main labels and title.
#' @param col group colours (\code{NULL} = automatic).
#' @param legend \code{TRUE} to draw a legend of the groups.
#' @param ... further arguments passed to
#'   \code{\link[graphics]{plot.default}}.
#' @return \code{NULL}, invisibly.
#' @examples
#' set.seed(2)
#' g <- f.data.generation(senario = 1)
#' plot_data(g$df, data.forma = 2)
#' @export

plot_data <- function(data, n.k = NULL, data.forma = 1,
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
    col <- if (K <= 8) seq_len(K) else grDevices::hcl.colors(K, "Dark 3")
  }
  col_courbe <- col[match(groupe, niveaux)]
  
  # ---- 3. Plot ---------------------------------------------------------------
  graphics::plot(NA, xlim = range(t_all), ylim = range(y_all),
       xlab = xlab, ylab = ylab, main = main, ...)
  
  if (data.forma == 1) {
    for (i in seq_along(data)) {
      graphics::lines(data[[i]], col = col_courbe[i], lty = 1, type = "b")
    }
  } else {
    rows <- split(seq_len(nrow(data)), factor(data$id, levels = ids))
    for (i in seq_along(ids)) {
      r <- rows[[i]]
      graphics::lines(data$x[r], data$y[r], col = col_courbe[i], lty = 1, type = "b")
    }
  }
  
  if (legend) {
    graphics::legend("topright", legend = paste("Group", niveaux),
           col = col, lty = 1, pch = 1, bty = "n")
  }
  
  
  invisible(NULL)
}
