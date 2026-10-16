# =========================================================================== #
#  Criteria and utilities                                                     #
# =========================================================================== #

.im_criteria_table <- function(fits, q, N, penalty) {
  rows <- lapply(fits, function(f) {
    if (is.null(f)) return(NULL)
    K <- f$K
    n_alpha <- if (is.null(f$n_alpha)) 0 else f$n_alpha
    dK <- (K - 1) + n_alpha + f$edf + q * (q + 1) / 2 + 1 + as.numeric(penalty)
    tt <- f$tau[f$tau > 0]
    ent <- -sum(tt * log(tt))
    bic <- -2 * f$ll + dK * log(N)
    data.frame(K = K, loglik = f$ll, edf = f$edf, n_par = dK, BIC = bic,
               ICL = bic + 2 * ent, entropy = ent, min_pi = min(f$pi))
  })
  do.call(rbind, rows)
}

.im_onehot <- function(lab, K, eps = 0.01) {
  tau <- matrix(eps / K, length(lab), K)
  tau[cbind(seq_along(lab), lab)] <- 1 - eps + eps / K
  tau
}

.im_inv <- function(M) chol2inv(chol(M))

.im_logdet <- function(M) 2 * sum(log(diag(chol(M))))

.im_psd_floor <- function(M, floor_abs) {
  M <- (M + t(M)) / 2
  e <- eigen(M, symmetric = TRUE)
  v <- pmax(e$values, floor_abs)
  e$vectors %*% (v * t(e$vectors))
}
