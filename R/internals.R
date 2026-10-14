# Part of the InfoMergeCluster package.
# Authors: Komlan Midodzi Noukpoape <noukpoape@gmail.com>, Cristian Preda

# Internal functions (not exported): priors, bases, batched linear algebra,
# sufficient statistics, EM, merge tree and criteria.

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

# =========================================================================== #
#  Bases and design statistics                                                #
# =========================================================================== #

.im_make_basis <- function(t, df, degree) {
  ord <- degree + 1L
  n_int <- df - ord
  if (n_int < 0) stop("'df' must be >= degree + 1.")
  lo <- min(t)
  hi <- max(t)
  interior <- numeric(0)
  if (n_int > 0) {
    pr <- seq(0, 1, length.out = n_int + 2)[-c(1, n_int + 2)]
    interior <- unique(stats::quantile(t, probs = pr, names = FALSE))
    interior <- interior[interior > lo & interior < hi]
  }
  knots <- c(rep(lo, ord), interior, rep(hi, ord))
  list(knots = knots, ord = ord, lo = lo, hi = hi, p = length(knots) - ord)
}

.im_basis_eval <- function(basis, t) {
  t <- pmin(pmax(t, basis$lo), basis$hi)
  splines::splineDesign(basis$knots, t, ord = basis$ord, outer.ok = TRUE)
}

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch algorithm)
.im_gauss_legendre <- function(n) {
  if (n == 1L) return(list(x = 0, w = 2))
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

# penalty matrix of the penalized B-splines (O'Sullivan):
#   P_jl = int_0^1 B_j^(m)(u) B_l^(m)(u) du,   u = (t - lo) / (hi - lo)
# On each interval between knots, B_j^(m) B_l^(m) is a polynomial of degree
# 2 (degree - m): a Gauss-Legendre quadrature with degree - m + 1 points
# is therefore exact. The change of variable t -> u multiplies the integral
# in t by (hi - lo)^(2m - 1).
.im_penalty_matrix <- function(basis, m) {
  m <- as.integer(m)
  if (m < 1L || m > basis$ord - 1L)
    stop("'pen_order' must satisfy 1 <= pen_order <= degree.")
  br <- unique(basis$knots)
  a <- br[-length(br)]
  h <- diff(br)
  gl <- .im_gauss_legendre(basis$ord - m)
  ng <- length(gl$x)
  x <- rep(a + h / 2, each = ng) + rep(h / 2, each = ng) * gl$x
  w <- rep(h / 2, each = ng) * gl$w
  Dm <- splines::splineDesign(basis$knots, x, ord = basis$ord,
                              derivs = rep(m, length(x)), outer.ok = TRUE)
  P <- crossprod(Dm * sqrt(w)) * (basis$hi - basis$lo)^(2 * m - 1)
  (P + t(P)) / 2
}

# random-effect design
.im_Z <- function(basis, re, t, B = NULL) {
  if (re == "spline") return(if (is.null(B)) .im_basis_eval(basis, t) else B)
  if (re == "intercept") return(matrix(1, length(t), 1))
  tt <- (pmin(pmax(t, basis$lo), basis$hi) - basis$lo) / (basis$hi - basis$lo) - 0.5
  cbind(1, tt)
}

# row-wise products: vec(U_r W_r') in column order (index j + (l-1) ncol(U))
.im_rowouter <- function(U, W) {
  a <- ncol(U)
  b <- ncol(W)
  U[, rep(seq_len(a), times = b), drop = FALSE] * W[, rep(seq_len(b), each = a), drop = FALSE]
}

# statistics that depend only on the design
.im_design_stats <- function(B, Zm, y, idx, N) {
  agg <- function(M) unname(rowsum(M, idx, reorder = TRUE))
  list(Q  = agg(.im_rowouter(B, B)),
       BZ = agg(.im_rowouter(B, Zm)),
       ZZ = agg(.im_rowouter(Zm, Zm)),
       cc = agg(B * y),
       Zy = agg(Zm * y),
       yy = as.vector(agg(matrix(y^2))),
       n  = tabulate(idx, N))
}

# groups subjects having the same number of measurements (vectorized computation)
.im_groups <- function(B, Zm, y, idx, N) {
  n <- tabulate(idx, N)
  start <- c(1L, cumsum(n)[-N] + 1L)
  lapply(sort(unique(n)), function(nv) {
    I <- which(n == nv)
    rows <- matrix(outer(start[I], 0:(nv - 1L), "+"), nrow = length(I))
    list(n = nv, I = I,
         X = lapply(seq_len(nv), function(j)
           cbind(B[rows[, j], , drop = FALSE], Zm[rows[, j], , drop = FALSE], y[rows[, j]])),
         Z = lapply(seq_len(nv), function(j) Zm[rows[, j], , drop = FALSE]))
  })
}

# =========================================================================== #
#  "Batched" linear algebra (vectorized over subjects)                        #
#  V: m x q x q array of symmetric positive definite matrices                 #
# =========================================================================== #

.im_bchol <- function(V) {
  q <- dim(V)[2]
  L <- array(0, dim(V))
  for (j in seq_len(q)) {
    s <- V[, j, j]
    if (j > 1) for (k in 1:(j - 1)) s <- s - L[, j, k]^2
    if (any(s <= 0)) stop("matrix not positive definite in .im_bchol")
    L[, j, j] <- sqrt(s)
    if (j < q) for (i in (j + 1):q) {
      s <- V[, i, j]
      if (j > 1) for (k in 1:(j - 1)) s <- s - L[, i, k] * L[, j, k]
      L[, i, j] <- s / L[, j, j]
    }
  }
  L
}

# solves L Z = X; X: list of q matrices m x r (row j of each system)
.im_bforward <- function(L, X) {
  q <- length(X)
  Z <- vector("list", q)
  for (j in seq_len(q)) {
    s <- X[[j]]
    if (j > 1) for (k in 1:(j - 1)) s <- s - L[, j, k] * Z[[k]]
    Z[[j]] <- s / L[, j, j]
  }
  Z
}

# solves L' X = Z
.im_bbackward <- function(L, Z) {
  q <- length(Z)
  X <- vector("list", q)
  for (j in q:1) {
    s <- Z[[j]]
    if (j < q) for (k in (j + 1):q) s <- s - L[, k, j] * X[[k]]
    X[[j]] <- s / L[, j, j]
  }
  X
}

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

# =========================================================================== #
#  InfoMergeCluster: merge tree based on the exact ICL                        #
# =========================================================================== #

.im_tree <- function(S, G, s2, lam, K0, K_max, a_dir, nstart, floors, ctrl, G_iter = 5) {
  N <- S$N
  p <- S$p
  st <- .im_subject_stats(S, G, s2)
  
  # (a) shrunken individual estimates (light ridge + smoothing penalty)
  Pl <- lam * ctrl$P
  Atot <- matrix(colSums(st$A), p, p)
  m0 <- solve(Atot + Pl + diag(1e-10 * mean(diag(Atot)), p), colSums(st$b))
  R0 <- diag(0.1 * mean(diag(Atot)) / N, p) + Pl
  r0 <- as.vector(R0 %*% m0)
  M <- array(sweep(st$A, 2, as.vector(R0), "+"), c(N, p, p))
  Lm <- .im_bchol(M)
  rhs <- lapply(seq_len(p), function(j) matrix(st$b[, j] + r0[j], ncol = 1))
  theta <- do.call(cbind, .im_bbackward(Lm, .im_bforward(Lm, rhs)))
  
  # (b) over-segmentation
  n_distinct <- nrow(unique(round(theta, 10)))
  K0 <- max(2L, min(K0, n_distinct - 1L))
  km <- suppressWarnings(stats::kmeans(theta, centers = K0, nstart = nstart,
                                       iter.max = 100))
  cl <- as.integer(factor(km$cluster))
  K0 <- max(cl)
  
  # (c) within micro-cluster covariance
  tau_h <- .im_onehot(cl, K0, 0)
  Gw <- G
  s2w <- s2
  fl_w <- list(G = 1e-10 * mean(diag(G)), s2 = floors$s2)
  for (it in seq_len(G_iter)) {
    stw <- .im_subject_stats(S, Gw, s2w)
    bw <- .im_mstep_beta(stw, tau_h, lam, ctrl$P)$beta
    cv <- .im_mstep_cov(S, stw, tau_h, bw, Gw, fl_w, ctrl$prior)
    Gw <- cv$G
    s2w <- cv$s2
  }
  st <- .im_subject_stats(S, Gw, s2w)
  
  # (d) micro-cluster statistics and empirical Bayes prior
  AC <- rowsum(st$A, cl, reorder = TRUE)
  bC <- rowsum(st$b, cl, reorder = TRUE)
  nC <- tabulate(cl, K0)
  Atot <- matrix(colSums(st$A), p, p)
  m0 <- solve(Atot + Pl + diag(1e-10 * mean(diag(Atot)), p), colSums(st$b))
  R0 <- diag(0.1 * mean(diag(Atot)) / N, p) + Pl
  r0 <- as.vector(R0 %*% m0)
  mu <- matrix(0, K0, p)
  for (k in seq_len(K0)) {
    mu[k, ] <- solve(R0 + matrix(AC[k, ], p, p), bC[k, ] + r0)
  }
  wts <- nC / sum(nC)
  m0 <- colSums(mu * wts)
  Sm <- crossprod(sweep(mu, 2, m0) * sqrt(wts))
  T0 <- Sm + diag(0.05 * mean(diag(Sm)) + 1e-12, p)
  T0inv <- .im_inv(T0)
  T0m <- as.vector(T0inv %*% m0)
  ldT0 <- .im_logdet(T0)
  mTm <- sum(m0 * T0m)
  
  Lfun <- function(Avec, bv) {
    R <- chol(T0inv + matrix(Avec, p, p))
    z <- backsolve(R, bv + T0m, transpose = TRUE)
    -0.5 * (ldT0 + 2 * sum(log(diag(R)))) + 0.5 * sum(z^2) - 0.5 * mTm
  }
  icl_part <- function(K, nv) {
    lgamma(K * a_dir) - K * lgamma(a_dir) + sum(lgamma(nv + a_dir)) - lgamma(N + K * a_dir)
  }
  Lc <- vapply(seq_len(K0), function(k) Lfun(AC[k, ], bC[k, ]), numeric(1))
  dfun <- function(k, l) {
    Lfun(AC[k, ] + AC[l, ], bC[k, ] + bC[l, ]) - Lc[k] - Lc[l] +
      lgamma(nC[k] + nC[l] + a_dir) - lgamma(nC[k] + a_dir) - lgamma(nC[l] + a_dir)
  }
  
  # (e) greedy merges
  D <- matrix(-Inf, K0, K0)
  for (k in 1:(K0 - 1)) for (l in (k + 1):K0) D[k, l] <- dfun(k, l)
  act <- rep(TRUE, K0)
  grp <- seq_len(K0)
  icl_path <- rep(NA_real_, K0)
  parts <- vector("list", K_max)
  K <- K0
  icl_path[K] <- sum(Lc) + icl_part(K, nC)
  if (K <= K_max) parts[[K]] <- grp
  while (K > 1) {
    w <- arrayInd(which.max(D), dim(D))
    k <- w[1, 1]
    l <- w[1, 2]
    AC[k, ] <- AC[k, ] + AC[l, ]
    bC[k, ] <- bC[k, ] + bC[l, ]
    nC[k] <- nC[k] + nC[l]
    nC[l] <- 0
    act[l] <- FALSE
    Lc[k] <- Lfun(AC[k, ], bC[k, ])
    grp[grp == l] <- k
    D[l, ] <- -Inf
    D[, l] <- -Inf
    for (j in which(act)) {
      if (j == k) next
      if (j < k) D[j, k] <- dfun(j, k) else D[k, j] <- dfun(k, j)
    }
    K <- K - 1L
    icl_path[K] <- sum(Lc[act]) + icl_part(K, nC[act])
    if (K <= K_max) parts[[K]] <- grp
  }
  parts <- lapply(parts, function(g) if (is.null(g)) NULL else as.integer(factor(g))[cl])
  
  list(parts = parts, icl_path = icl_path, K_tree = which.max(icl_path), K0 = K0,
       micro = cl, G_within = Gw, s2_within = s2w, T0 = T0, m0 = m0)
}

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
