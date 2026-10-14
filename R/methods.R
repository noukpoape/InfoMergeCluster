# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi Noukpoape <noukpoape@gmail.com>

# =========================================================================== #
#  S3 methods and helper functions                                            #
# =========================================================================== #

#' Print an InfoMergeCluster fit
#'
#' Prints the number of clusters, the proportions, the variance parameters,
#' the priors and the table of criteria of a fit returned by
#' \code{\link{infomergecluster}}.
#'
#' @param x an object of class \code{"infomergecluster"}.
#' @param ... unused.
#' @return \code{x}, invisibly.
#' @export
print.infomergecluster <- function(x, ...) {
  cat("InfoMergeCluster: clustering of sparse longitudinal data\n")
  cat(sprintf("  %d subjects, %d measurements, B-spline basis p = %d, random effect: %s\n",
              x$N, x$n_obs, x$basis$p, x$re))
  if (isTRUE(x$K_fixed)) {
    cat(sprintf("  K imposed by the user: %d (not estimated)\n", x$K))
  } else {
    cat(sprintf("  Selected K (%s): %d\n", x$criterion, x$K))
  }
  cat(sprintf("  Proportions: %s\n", paste(sprintf("%.3f", x$pi), collapse = ", ")))
  cat(sprintf("  sigma2 = %.4g; trace(G) = %.4g; lambda = %.3g; effective df (total) = %.1f\n",
              x$sigma2, sum(diag(x$G)), x$lambda, x$edf))
  if (!is.null(x$prior)) {
    cat(sprintf("  Priors: G ~ IW(nu0 = %.3g, Psi0 [trace %.3g]); sigma2 ~ IG(%.3g, %.3g)\n",
                x$prior$nu0, sum(diag(x$prior$Psi0)), x$prior$a_sigma, x$prior$b_sigma))
    cat(sprintf("    nugget = %.4g (%s), sill = %.4g; prior modes: sigma2 = %.4g, G = %.4g I\n",
                x$prior$nugget, x$prior$nugget_source, x$variogram$sill,
                x$prior$b_sigma / (x$prior$a_sigma + 1),
                mean(diag(x$prior$Psi0)) / (x$prior$nu0 + x$prior$q + 1)))
  }
  if (!is.null(x$gamma)) {
    cat(sprintf("  Covariates (centered, cluster-specific effects gamma_k): %s\n",
                paste(rownames(x$gamma), collapse = ", ")))
    print(round(x$gamma, 4))
  }
  if (identical(x$weights, "covariates") && !is.null(x$alpha)) {
    cat("  Weights depending on the covariates (multinomial logistic, reference C1):\n")
    print(round(x$alpha, 4))
  }
  cat(sprintf("  computing time = %.2f s\n\n", x$time_sec))
  if (!is.null(x$var_y) && x$sigma2 < 1e-3 * x$var_y) {
    cat("  WARNING: near-zero residual variance (sigma2 < 0.1% of var(y)).\n",
        " In this regime (data without measurement error), several partitions have\n",
        " nearly equal criteria: the selected K may vary with the seed or a_dir.\n",
        " See infomergecluster_stability() or increase df.\n\n")
  }
  if (!is.null(x$stability)) {
    st <- x$stability
    cat(sprintf("  Stability over %d seeds: mean ARI between partitions = %.2f\n",
                length(st$seeds), st$mean_ari))
    cat("  Frequency of the K obtained:",
        paste(sprintf("K=%s: %d", names(st$K_freq), as.integer(st$K_freq)), collapse = ", "),
        "\n\n")
  }
  tab <- x$table
  num <- c("loglik", "edf", "n_par", "BIC", "ICL", "entropy")
  tab[num] <- lapply(tab[num], round, 2)
  tab$min_pi <- round(tab$min_pi, 3)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Mean trajectories of the clusters
#'
#' Evaluates the mean curve of each cluster on a time grid, with an optional
#' pointwise envelope.
#'
#' With \code{interval = "confidence"}, the band reflects the uncertainty on
#' the mean curve, \eqn{f_k(t) \pm z \sqrt{B(t)' H_k^{-1} B(t)}}. With
#' \code{interval = "prediction"}, it contains the measurements of a subject
#' of the cluster, with variance \eqn{B' H_k^{-1} B + Z G Z' + \sigma^2}.
#'
#' @param object an object of class \code{"infomergecluster"}.
#' @param times time grid (default: 200 points over the observed range).
#' @param interval \code{"none"}, \code{"confidence"} or \code{"prediction"}.
#' @param level level of the envelope.
#' @param covariates \code{NULL} (curves at the mean covariate values) or a
#'   numeric vector of covariate values on the original scale, in the order
#'   of \code{rownames(object$gamma)} (dummy-coded for factors).
#' @param ... unused.
#' @return A data.frame with a \code{time} column, one column per cluster
#'   (\code{C1}, \code{C2}, ...) and, if an envelope is requested, the bounds
#'   \code{<Ck>_lwr} and \code{<Ck>_upr}.
#' @export
predict.infomergecluster <- function(object, times = NULL,
                                     interval = c("none", "confidence", "prediction"),
                                     level = 0.95, covariates = NULL, ...) {
  interval <- match.arg(interval)
  if (is.null(times)) times <- seq(object$basis$lo, object$basis$hi, length.out = 200)
  Bt <- .im_basis_eval(object$basis, times)
  r <- if (is.null(object$gamma)) 0 else nrow(object$gamma)
  xc <- numeric(r)
  if (r > 0 && !is.null(covariates)) {
    if (length(covariates) != r)
      stop(sprintf("'covariates' must have %d value(s), in the order: %s", r,
                   paste(rownames(object$gamma), collapse = ", ")))
    xc <- as.numeric(covariates) - as.numeric(object$covariates$center)
  }
  At <- if (r > 0) cbind(Bt, matrix(xc, nrow(Bt), r, byrow = TRUE)) else Bt
  coef <- if (r > 0) rbind(object$beta, object$gamma) else object$beta
  fit <- At %*% coef
  labs <- colnames(object$beta)
  out <- data.frame(time = times, fit)
  names(out)[-1] <- labs
  if (interval == "none") return(out)
  if (is.null(object$Vbeta))
    stop("Object fitted with an earlier version: rerun infomergecluster() to obtain the envelopes.")
  if (level <= 0 || level >= 1) stop("'level' must be in ]0, 1[.")
  z <- stats::qnorm(1 - (1 - level) / 2)
  var_extra <- 0
  if (interval == "prediction") {
    Zt <- .im_Z(object$basis, object$re, times, Bt)
    var_extra <- rowSums((Zt %*% object$G) * Zt) + object$sigma2
  }
  for (k in seq_along(labs)) {
    se <- sqrt(pmax(rowSums((At %*% object$Vbeta[[k]]) * At), 0) + var_extra)
    out[[paste0(labs[k], "_lwr")]] <- fit[, k] - z * se
    out[[paste0(labs[k], "_upr")]] <- fit[, k] + z * se
  }
  out
}

#' Predicted individual trajectories
#'
#' Computes the predicted trajectory of each subject (BLUP of the random
#' effect, conditional on its MAP cluster) on a time grid. With covariates,
#' they are held at the subject's mean values.
#'
#' @param object an object of class \code{"infomergecluster"} fitted with
#'   \code{keep_data = TRUE}.
#' @param times time grid (default: 100 points over the observed range).
#' @return A matrix with one row per subject and one column per time.
#' @export
fitted_curves <- function(object, times = NULL) {
  if (is.null(object$data)) stop("Rerun infomergecluster() with keep_data = TRUE.")
  if (is.null(times)) times <- seq(object$basis$lo, object$basis$hi, length.out = 100)
  d <- object$data
  ids <- names(object$cluster)
  Bt <- .im_basis_eval(object$basis, times)
  Zt <- .im_Z(object$basis, object$re, times, Bt)
  G <- object$G
  s2 <- object$sigma2
  rows <- split(seq_len(nrow(d)), factor(d$id, levels = unique(d$id)))
  res <- matrix(NA_real_, length(ids), length(times), dimnames = list(ids, NULL))
  for (i in seq_along(ids)) {
    r_i <- rows[[ids[i]]]
    Bi <- .im_basis_eval(object$basis, d$t[r_i])
    Zi <- .im_Z(object$basis, object$re, d$t[r_i], Bi)
    Vi <- Zi %*% G %*% t(Zi) + diag(s2, nrow(Zi))
    bk <- object$beta[, object$cluster[i]]
    shift_obs <- 0
    shift_grid <- 0
    if (!is.null(object$gamma)) {
      Xi <- as.matrix(d[r_i, object$covariates$names, drop = FALSE])
      gk <- object$gamma[, object$cluster[i]]
      shift_obs <- as.vector(Xi %*% gk)
      shift_grid <- sum(colMeans(Xi) * gk)   # covariates held at the subject mean
    }
    u <- G %*% t(Zi) %*% solve(Vi, d$y[r_i] - Bi %*% bk - shift_obs)
    res[i, ] <- Bt %*% bk + shift_grid + Zt %*% u
  }
  res
}

#' Plot an InfoMergeCluster fit
#'
#' Draws the individual trajectories colored by cluster and the mean curve
#' of each cluster, with an optional envelope.
#'
#' @param x an object of class \code{"infomergecluster"}.
#' @param max_subjects maximum number of subjects drawn (random subsample beyond).
#' @param envelope \code{TRUE} to draw an envelope around the mean curves.
#' @param envelope_type \code{"confidence"} (uncertainty on the mean) or
#'   \code{"prediction"} (spread of the subjects in the cluster).
#' @param xlab,ylab axis labels.
#' @param main title; a format with \code{\%d} for K and \code{\%s} for the
#'   criterion (or "imposed").
#' @param level envelope level.
#' @param show_data whether to draw the individual trajectories.
#' @param legend \code{TRUE} to draw the legend.
#' @param legend_pos legend position (\code{"topleft"}, \code{"topright"}, ...).
#' @param ... further arguments passed to \code{plot}.
#' @return \code{x}, invisibly.
#' @export
plot.infomergecluster <- function(x, max_subjects = 300, envelope = FALSE,
                                  envelope_type = c("confidence", "prediction"),
                                  xlab = "Time (days)", ylab = "Weight",
                                  main = "InfoMergeCluster: K = %d (%s)",
                                  level = 0.95, show_data = TRUE,
                                  legend = TRUE, legend_pos = "topleft", ...) {
  envelope_type <- match.arg(envelope_type)
  cols <- grDevices::hcl.colors(max(x$K, 2), "Dark 3")[seq_len(x$K)]
  labs <- colnames(x$beta)
  mc <- predict(x, interval = if (envelope) envelope_type else "none", level = level)
  yl <- range(mc[, -1])
  if (show_data && !is.null(x$data)) yl <- range(yl, x$data$y)
  graphics::plot(NA, xlim = c(x$basis$lo, x$basis$hi), ylim = yl,
                 xlab = xlab, ylab = ylab,
                 main = sprintf(main, x$K, if (isTRUE(x$K_fixed)) "imposed" else x$criterion),
                 ...)
  if (show_data && !is.null(x$data)) {
    ids <- names(x$cluster)
    show <- if (length(ids) > max_subjects) sample(ids, max_subjects) else ids
    d <- x$data[as.character(x$data$id) %in% show, ]
    for (s in split(d, as.character(d$id))) {
      k <- x$cluster[as.character(s$id[1])]
      graphics::lines(s$t, s$y, col = grDevices::adjustcolor(cols[k], 0.25))
      graphics::points(s$t, s$y, col = grDevices::adjustcolor(cols[k], 0.4),
                       pch = 16, cex = 0.4)
    }
  }
  if (envelope) {
    for (k in seq_len(x$K)) {
      lo <- mc[[paste0(labs[k], "_lwr")]]
      up <- mc[[paste0(labs[k], "_upr")]]
      graphics::polygon(c(mc$time, rev(mc$time)), c(lo, rev(up)),
                        col = grDevices::adjustcolor(cols[k], 0.25), border = NA)
      graphics::lines(mc$time, lo, col = cols[k], lty = 2, lwd = 1)
      graphics::lines(mc$time, up, col = cols[k], lty = 2, lwd = 1)
    }
  }
  for (k in seq_len(x$K)) graphics::lines(mc$time, mc[[labs[k]]], col = cols[k], lwd = 3)
  if (isTRUE(legend)) {
    leg_title <- if (envelope) {
      sprintf("%g%% %s envelope", 100 * level,
              if (envelope_type == "confidence") "confidence" else "prediction")
    } else NULL
    graphics::legend(legend_pos, legend = sprintf("%s (%.0f%%)", labs, 100 * x$pi),
                     col = cols, lwd = 3, bty = "n", title = leg_title)
  }
  invisible(x)
}

#' Stability of the clustering across random seeds
#'
#' Reruns \code{\link{infomergecluster}} with several seeds and returns the
#' fit with the best criterion, supplemented with a stability diagnostic.
#'
#' @param data data, as in \code{\link{infomergecluster}}.
#' @param ... arguments passed to \code{\link{infomergecluster}} (except \code{seed}).
#' @param seeds vector of seeds.
#' @param verbose print the current seed.
#' @return The best fit (class \code{"infomergecluster"}) with an extra
#'   element \code{stability}: the K of each run, their frequencies, the
#'   criteria, the matrix of adjusted Rand indices between runs and its mean.
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
