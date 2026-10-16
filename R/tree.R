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
