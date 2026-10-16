# =========================================================================== #
#  Priors on the variance structure (definition "priors")                     #
# =========================================================================== #

# empirical within-subject variogram of the residuals (Diggle, 1988)
#   pairs (j < k) of the same subject: u = |t_j - t_k|, v = (r_j - r_k)^2 / 2
#   nugget: intercept of the regression of v on u, restricted to the
#           lags u <= quantile lag_prob (at least min_pairs pairs);
#           mean of v over these pairs if the intercept is not positive.
#   sill:   var(r).
.im_variogram <- function(idx, t, r, lag_prob = 0.25, min_pairs = 10) {
  sill <- stats::var(r)
  grp <- split(seq_along(r), idx)
  grp <- grp[lengths(grp) >= 2]
  if (length(grp) == 0)
    return(list(nugget = NA_real_, sill = sill, n_pairs = 0L, lag_max = NA_real_))
  pairs <- do.call(rbind, lapply(grp, function(ix) {
    cm <- utils::combn(length(ix), 2)
    a <- ix[cm[1, ]]
    b <- ix[cm[2, ]]
    cbind(abs(t[a] - t[b]), 0.5 * (r[a] - r[b])^2)
  }))
  u <- pairs[, 1]
  v <- pairs[, 2]
  n_pairs <- length(u)
  if (n_pairs < min_pairs)
    return(list(nugget = NA_real_, sill = sill, n_pairs = n_pairs, lag_max = NA_real_))
  u_max <- stats::quantile(u, lag_prob, names = FALSE)
  sel <- which(u <= u_max)
  if (length(sel) < min_pairs) sel <- order(u)[seq_len(min_pairs)]
  nug <- NA_real_
  if (length(unique(u[sel])) >= 2) {
    cf <- stats::lm.fit(cbind(1, u[sel]), v[sel])$coefficients
    nug <- cf[1]
  }
  if (!is.finite(nug) || nug <= 0) nug <- mean(v[sel])
  list(nugget = unname(min(nug, sill)), sill = sill, n_pairs = n_pairs,
       lag_max = max(u[sel]))
}

# prior hyperparameters, calibrated on the variogram (not user-settable)
#   G ~ IW_q(nu0 = q, Psi0),  s2 ~ IG(a_sigma = 1, b_sigma)
#   prior mode of s2 = nugget; prior mode of G = g I_q, with g such that the
#   mean variance of Z(t)'u at the observed times equals
#   g_share (sill - nugget), g_share = 0.5
#   vario: output of .im_variogram; Zm: random-effect design
.im_prior <- function(q, vy, vario, Zm) {
  g_share <- 0.5
  nu0 <- q
  a <- 1
  # estimated nugget; fallback if the variogram is unusable
  if (is.finite(vario$nugget)) {
    nugget <- vario$nugget
    src <- "variogram"
  } else {
    nugget <- 0.05 * vario$sill
    src <- "fallback 5% of the sill (too few within-subject pairs)"
  }
  nugget <- max(nugget, 1e-8 * vy)  # safeguard on the hyperparameter (b_sigma > 0)
  proc <- max(vario$sill - nugget, 1e-8 * vy)
  g_mode <- g_share * proc / mean(rowSums(Zm^2))
  Psi0 <- diag((nu0 + q + 1) * g_mode, q)
  b <- (a + 1) * nugget
  list(nu0 = nu0, Psi0 = Psi0, a_sigma = a, b_sigma = b, q = q,
       nugget = nugget, nugget_source = src, g_mode = g_mode, g_share = g_share)
}

# log p0(G, s2), up to a constant
#   -(nu0 + q + 1)/2 log|G| - tr(Psi0 G^-1)/2 - (a + 1) log s2 - b / s2
.im_log_prior <- function(G, s2, pr) {
  q <- nrow(G)
  -0.5 * (pr$nu0 + q + 1) * .im_logdet(G) - 0.5 * sum(pr$Psi0 * .im_inv(G)) -
    (pr$a_sigma + 1) * log(s2) - pr$b_sigma / s2
}
