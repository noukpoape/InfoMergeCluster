# =========================================================================== #
#  Sufficient information statistics (Proposition 1)                          #
#  V_i = Z_i G Z_i' + s2 I; batched Cholesky, W = L^{-1} [B_i | Z_i | y_i]    #
#  A_i = B'V^-1 B, b_i = B'V^-1 y, and for the M-step of G:                   #
#  ZVB_i = Z'V^-1 B, ZVZ_i = Z'V^-1 Z, ZVy_i = Z'V^-1 y                       #
# =========================================================================== #

.im_subject_stats <- function(S, G, s2) {
  N <- S$N
  p <- S$p
  q <- S$q
  ip <- seq_len(p)
  iq <- p + seq_len(q)
  A <- matrix(0, N, p * p)
  b <- matrix(0, N, p)
  ZVB <- matrix(0, N, q * p)
  ZVZ <- matrix(0, N, q * q)
  ZVy <- matrix(0, N, q)
  c0 <- numeric(N)
  for (g in S$groups) {
    n <- g$n
    m <- length(g$I)
    ZG <- lapply(g$Z, function(Zj) Zj %*% G)
    V <- array(0, c(m, n, n))
    for (j in seq_len(n)) for (l in seq_len(j)) {
      v <- rowSums(ZG[[j]] * g$Z[[l]])
      if (j == l) v <- v + s2
      V[, j, l] <- v
      V[, l, j] <- v
    }
    L <- .im_bchol(V)
    W <- .im_bforward(L, g$X)
    Ag <- matrix(0, m, p * p)
    bg <- matrix(0, m, p)
    zvb <- matrix(0, m, q * p)
    zvz <- matrix(0, m, q * q)
    zvy <- matrix(0, m, q)
    yVy <- numeric(m)
    ld <- numeric(m)
    for (j in seq_len(n)) {
      wb <- W[[j]][, ip, drop = FALSE]
      wz <- W[[j]][, iq, drop = FALSE]
      wy <- W[[j]][, p + q + 1]
      Ag <- Ag + .im_rowouter(wb, wb)
      zvb <- zvb + .im_rowouter(wz, wb)
      zvz <- zvz + .im_rowouter(wz, wz)
      bg <- bg + wb * wy
      zvy <- zvy + wz * wy
      yVy <- yVy + wy^2
      ld <- ld + 2 * log(L[, j, j])
    }
    A[g$I, ] <- Ag
    b[g$I, ] <- bg
    ZVB[g$I, ] <- zvb
    ZVZ[g$I, ] <- zvz
    ZVy[g$I, ] <- zvy
    c0[g$I] <- -0.5 * (n * log(2 * pi) + ld + yVy)
  }
  list(A = A, b = b, c0 = c0, ZVB = ZVB, ZVZ = ZVZ, ZVy = ZVy)
}

# =========================================================================== #
#  EM on the sufficient statistics (Proposition 2)                            #
# =========================================================================== #

.im_estep <- function(st, beta, pik) {
  K <- ncol(beta)
  N <- nrow(st$b)
  L <- matrix(0, N, K)
  for (k in seq_len(K)) {
    bk <- beta[, k]
    lp <- if (is.matrix(pik)) log(pmax(pik[, k], 1e-300)) else log(max(pik[k], 1e-300))
    L[, k] <- lp + st$c0 + as.vector(st$b %*% bk) -
      0.5 * as.vector(st$A %*% as.vector(tcrossprod(bk)))
  }
  mx <- L[cbind(seq_len(N), max.col(L, ties.method = "first"))]
  E <- exp(L - mx)
  s <- rowSums(E)
  list(tau = E / s, ll = sum(mx + log(s)))
}

.im_mstep_beta <- function(st, tau, lam = 0, P = NULL) {
  p <- ncol(st$b)
  K <- ncol(tau)
  AT <- crossprod(st$A, tau)                    # p^2 x K
  bT <- crossprod(st$b, tau)                    # p   x K
  ridge <- 1e-8 * max(mean(diag(matrix(colMeans(st$A), p, p))), 1e-12)
  beta <- matrix(0, p, K)
  H <- vector("list", K)
  for (k in seq_len(K)) {
    H[[k]] <- matrix(AT[, k], p, p) + diag(ridge, p)
    if (!is.null(P) && lam > 0) H[[k]] <- H[[k]] + lam * P
    beta[, k] <- solve(H[[k]], bT[, k])
  }
  list(beta = beta, pi = colMeans(tau), H = H, AT = AT)
}

# M-step of the concomitant weights pi_k(x_i) = softmax(W_i alpha)_k:
#   maximizes sum_ik tau_ik log pi_ik - ridge/2 |slopes|^2 (cluster 1 =
#   reference), by Newton steps with step halving; W: N x (1 + r).
.im_mstep_pi <- function(tau, W, ridge = 1e-4, iters = 25, tol = 1e-9) {
  N <- nrow(tau)
  K <- ncol(tau)
  dw <- ncol(W)
  A <- matrix(0, dw, K)
  if (K == 1) return(list(alpha = A, P = matrix(1, N, 1)))
  nk <- colSums(tau)
  A[1, ] <- log(pmax(nk, 1e-10) / max(nk[1], 1e-10))
  A[1, 1] <- 0
  sm <- function(A) {
    E <- W %*% A
    E <- E - apply(E, 1, max)
    Pm <- exp(E)
    Pm / rowSums(Pm)
  }
  obj <- function(A) {
    E <- W %*% A
    mx <- apply(E, 1, max)
    sum(tau * E) - sum(rowSums(tau) * (mx + log(rowSums(exp(E - mx))))) -
      0.5 * ridge * sum(A[-1, -1]^2)
  }
  Km <- K - 1
  pen <- rep(c(0, rep(ridge, dw - 1)), Km)
  f_old <- obj(A)
  for (it in seq_len(iters)) {
    Pm <- sm(A)
    R <- tau - Pm
    gA <- crossprod(W, R[, -1, drop = FALSE])
    gA[-1, ] <- gA[-1, ] - ridge * A[-1, -1, drop = FALSE]
    H <- matrix(0, dw * Km, dw * Km)
    for (j in seq_len(Km)) for (l in j:Km) {
      wgt <- Pm[, j + 1] * ((j == l) - Pm[, l + 1])
      blk <- crossprod(W * wgt, W)
      ij <- (j - 1) * dw + seq_len(dw)
      il <- (l - 1) * dw + seq_len(dw)
      H[ij, il] <- blk
      H[il, ij] <- t(blk)
    }
    H <- H + diag(pen + 1e-8, dw * Km)
    step <- matrix(solve(H, as.vector(gA)), dw, Km)
    stp <- 1
    repeat {
      An <- A
      An[, -1] <- A[, -1] + stp * step
      f_new <- obj(An)
      if (is.finite(f_new) && f_new >= f_old - 1e-12) break
      stp <- stp / 2
      if (stp < 1e-4) { An <- A; f_new <- f_old; break }
    }
    A <- An
    done <- abs(f_new - f_old) < tol * (1 + abs(f_old))
    f_old <- f_new
    if (done) break
  }
  list(alpha = A, P = sm(A))
}

# weights used by the E-step: vector pi_k, or N x K matrix pi_k(x_i)
.im_weights <- function(tau, ms, ctrl) {
  if (is.null(ctrl$Wc)) ms$pi else .im_mstep_pi(tau, ctrl$Wc)$P
}

# update of the smoothing parameter (Fellner-Schall):
#   lambda <- sum_k [rank(P) - lambda tr(H_k^{-1} P)] / sum_k beta_k' P beta_k
.im_lam_update <- function(beta, H, lam, P, rk) {
  num <- 0
  den <- 0
  for (k in seq_len(ncol(beta))) {
    num <- num + rk - lam * sum(diag(solve(H[[k]], P)))
    den <- den + sum(beta[, k] * (P %*% beta[, k]))
  }
  min(max(max(num, 1e-8) / max(den, 1e-300), 1e-8), 1e10)
}

# sum_{j,l} M[, j + (l-1)p] W[, j] W[, l]
.im_rowquad <- function(M, W) {
  p <- ncol(W)
  out <- numeric(nrow(W))
  for (l in seq_len(p)) {
    cols <- (l - 1L) * p + seq_len(p)
    out <- out + W[, l] * rowSums(M[, cols, drop = FALSE] * W)
  }
  out
}

# M-step for (G, s2), vectorized
#   E[u_i | y_i, z_i = k] = G (Z'V^-1 y - Z'V^-1 B beta_k),  Var = G - G Z'V^-1 Z G
#   posterior mode:  G <- (N Sigma_G + Psi0) / (N + nu0 + q + 1),
#                    s2 <- (sse + 2 b) / (sum(n_i) + 2 a + 2)  (> 0 by construction)
#   then numerical safeguard: eigenvalues of G >= floors$G.
.im_mstep_cov <- function(S, st, tau, beta, G, floors, prior) {
  q <- S$q
  N <- S$N
  Iq <- diag(q)
  Mg <- matrix(0, q, q)
  sse <- 0
  for (k in seq_len(ncol(beta))) {
    sel <- which(tau[, k] > 1e-12)
    if (length(sel) == 0) next
    tk <- tau[sel, k]
    bk <- beta[, k]
    g <- st$ZVy[sel, , drop = FALSE] - st$ZVB[sel, , drop = FALSE] %*% kronecker(bk, Iq)
    Mg <- Mg + crossprod(g * tk, g)
    m <- g %*% G
    res <- S$yy[sel] - 2 * as.vector(S$cc[sel, , drop = FALSE] %*% bk) -
      2 * rowSums(S$Zy[sel, , drop = FALSE] * m) +
      as.vector(S$Q[sel, , drop = FALSE] %*% as.vector(tcrossprod(bk))) +
      2 * rowSums((S$BZ[sel, , drop = FALSE] %*% kronecker(Iq, bk)) * m) +
      .im_rowquad(S$ZZ[sel, , drop = FALSE], m)
    sse <- sse + sum(tk * res)
  }
  Euu <- G %*% Mg %*% G + N * G - G %*% matrix(colSums(st$ZVZ), q, q) %*% G
  trGZZ <- as.vector(S$ZZ %*% as.vector(G))
  trGZVZGZZ <- rowSums((st$ZVZ %*% kronecker(G, G)) * S$ZZ)
  sse <- sse + sum(trGZZ - trGZVZGZZ)
  G_new <- (Euu + prior$Psi0) / (N + prior$nu0 + q + 1)
  s2_new <- (sse + 2 * prior$b_sigma) / (sum(S$n) + 2 * prior$a_sigma + 2)
  list(G = .im_psd_floor(G_new, floors$G), s2 = max(s2_new, floors$s2))
}

# parameterization of (G, s2, lambda) for SQUAREM extrapolation
.im_pack <- function(G, s2, lam) {
  L <- t(chol(G))
  c(L[lower.tri(L, diag = TRUE)], log(s2), log(max(lam, 1e-300)))
}

.im_unpack <- function(th, q, floors, penalty) {
  nt <- length(th)
  L <- matrix(0, q, q)
  L[lower.tri(L, diag = TRUE)] <- th[seq_len(nt - 2)]
  lam <- if (penalty) min(max(exp(th[nt]), 1e-8), 1e10) else 0
  list(G = tcrossprod(L) + diag(floors$G, q), s2 = max(exp(th[nt - 1]), floors$s2),
       lam = lam)
}

# one application of the EM map: (G, s2, lambda) -> (G', s2', lambda')
.im_Fmap <- function(S, th, tau, ctrl) {
  cs <- .im_unpack(th, S$q, ctrl$floors, ctrl$penalty)
  st <- .im_subject_stats(S, cs$G, cs$s2)
  ms <- .im_mstep_beta(st, tau, cs$lam, ctrl$P)
  pw <- .im_weights(tau, ms, ctrl)
  for (it in seq_len(ctrl$inner_cap)) {
    es <- .im_estep(st, ms$beta, pw)
    tau <- es$tau
    if (it < ctrl$inner_cap) {
      ms <- .im_mstep_beta(st, tau, cs$lam, ctrl$P)
      pw <- .im_weights(tau, ms, ctrl)
    }
  }
  pen <- 0
  lam_new <- 0
  if (ctrl$penalty) {
    pen <- 0.5 * cs$lam * sum(ms$beta * (ctrl$P %*% ms$beta))
    lam_new <- .im_lam_update(ms$beta, ms$H, cs$lam, ctrl$P, ctrl$rk)
  }
  cv <- .im_mstep_cov(S, st, tau, ms$beta, cs$G, ctrl$floors, ctrl$prior)
  # penalized objective (eq:pen): l_N - roughness penalty + log p0(G, s2)
  obj <- es$ll - pen + .im_log_prior(cs$G, cs$s2, ctrl$prior)
  list(th = .im_pack(cv$G, cv$s2, lam_new), ll = obj, tau = tau)
}

# EM accelerated by SQUAREM (Varadhan & Roland, 2008) with safeguard
.im_em <- function(S, K, tau, G, s2, lam, max_cycles, ctrl, max_iter = 500) {
  th0 <- .im_pack(G, s2, lam)
  ll_old <- -Inf
  n_eval <- 0L
  for (cyc in seq_len(max_cycles)) {
    F1 <- .im_Fmap(S, th0, tau, ctrl)
    F2 <- .im_Fmap(S, F1$th, F1$tau, ctrl)
    r <- F1$th - th0
    v <- F2$th - F1$th - r
    nv <- sqrt(sum(v^2))
    alpha <- if (nv > 0) -sqrt(sum(r^2)) / nv else -1
    alpha <- max(min(alpha, -1), -50)
    F3 <- tryCatch(.im_Fmap(S, th0 - 2 * alpha * r + alpha^2 * v, F2$tau, ctrl),
                   error = function(e) NULL)
    n_eval <- n_eval + 3L
    if (!is.null(F3) && is.finite(F3$ll) && F3$ll >= F2$ll) {
      th0 <- F3$th; tau <- F3$tau; ll <- F3$ll
    } else {
      th0 <- F2$th; tau <- F2$tau; ll <- F2$ll
    }
    if (abs(ll - ll_old) < ctrl$tol * abs(ll)) break
    ll_old <- ll
  }
  # polishing: EM on (beta, pi) with (G, s2, lambda) fixed
  cs <- .im_unpack(th0, S$q, ctrl$floors, ctrl$penalty)
  st <- .im_subject_stats(S, cs$G, cs$s2)
  ms <- .im_mstep_beta(st, tau, cs$lam, ctrl$P)
  pw <- .im_weights(tau, ms, ctrl)
  ll <- -Inf
  for (it in seq_len(max_iter)) {
    es <- .im_estep(st, ms$beta, pw)
    tau <- es$tau
    conv <- abs(es$ll - ll) < 1e-9 * abs(es$ll)
    ll <- es$ll
    if (conv) break
    ms <- .im_mstep_beta(st, tau, cs$lam, ctrl$P)
    pw <- .im_weights(tau, ms, ctrl)
  }
  alpha <- if (is.null(ctrl$Wc)) NULL else .im_mstep_pi(tau, ctrl$Wc)$alpha
  p <- S$p
  edf <- sum(vapply(seq_len(K), function(k)
    sum(diag(solve(ms$H[[k]], matrix(ms$AT[, k], p, p)))), numeric(1)))
  # posterior covariance of the coefficients (penalized B-splines, fixed assignments)
  Vbeta <- lapply(ms$H, function(H) chol2inv(chol((H + t(H)) / 2)))
  # ll: log-likelihood l_N (BIC/ICL); log_prior: log p0 at the posterior mode
  list(K = K, beta = ms$beta, pi = colMeans(tau), G = cs$G, s2 = cs$s2, lam = cs$lam,
       edf = edf, tau = tau, ll = ll, log_prior = .im_log_prior(cs$G, cs$s2, ctrl$prior),
       alpha = alpha,
       n_alpha = if (is.null(ctrl$Wc)) 0 else (K - 1) * (ncol(ctrl$Wc) - 1),
       n_eval = n_eval, Vbeta = Vbeta)
}
