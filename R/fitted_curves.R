# predicted individual trajectories (BLUP conditional on the MAP cluster)

#' Predicted individual trajectories
#'
#' Predicts the trajectory of each subject (BLUP of the random effect,
#' conditional on its MAP cluster). With covariates, they are held at the
#' subject's mean value.
#'
#' @param object an object of class \code{"infomergecluster"} fitted with
#'   \code{keep_data = TRUE}.
#' @param times evaluation times (\code{NULL} = 100 points spanning the
#'   observed range).
#' @return A subjects x times matrix.
#' @examples
#' set.seed(1)
#' g <- f.data.generation(senario = 1, n.k = c(30, 40), n.cts = c(0, 2))
#' fit <- infomergecluster(g$liste, K = 2, verbose = FALSE)
#' fc <- fitted_curves(fit)
#' dim(fc)
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
