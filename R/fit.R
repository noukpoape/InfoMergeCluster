# full fit for a given random-effect structure
.im_fit_one <- function(d, re, K_max, K_fixed = NULL, df, degree, penalty, pen_order,
                        cov_names = character(0), cov_center = numeric(0),
                        weights = "constant", min_size, K0,
                        criterion, a_dir,
                        n_outer, quick_cycles, full_cycles, inner_cap, tol,
                        nstart, verbose) {
  ids <- unique(d$id)
  idx <- match(d$id, ids)
  N <- length(ids)
  
  basis <- .im_make_basis(d$t, df, degree)
  B <- .im_basis_eval(basis, d$t)
  Zm <- .im_Z(basis, re, d$t, B)
  p0 <- ncol(B)                      # spline part
  r <- length(cov_names)             # covariates
  if (r > 0) B <- cbind(B, as.matrix(d[cov_names]))   # augmented design [B | X]
  p <- ncol(B)
  q <- ncol(Zm)
  # subject-level design of the concomitant weights: (1, mean covariates)
  Wc <- NULL
  if (weights == "covariates") {
    Xs <- rowsum(as.matrix(d[cov_names]), idx, reorder = TRUE) / tabulate(idx, N)
    Wc <- cbind(1, unname(Xs))
  }
  
  S <- .im_design_stats(B, Zm, d$y, idx, N)
  S$groups <- .im_groups(B, Zm, d$y, idx, N)
  S$p <- p
  S$q <- q
  S$N <- N
  
  # imposed K: only the model with K_fixed clusters is fitted
  if (!is.null(K_fixed)) K_max <- as.integer(K_fixed)
  if (is.null(K0)) K0 <- min(100L, max(3L * K_max, ceiling(2 * sqrt(N))))
  K0 <- as.integer(min(K0, floor(N / 3)))
  if (!is.null(K_fixed)) K0 <- max(K0, K_max)
  K_max <- as.integer(min(K_max, K0))
  
  vy <- stats::var(d$y)
  P <- matrix(0, p, p)                # penalty diag(P_spline, 0): gamma not smoothed
  if (penalty) {
    P[seq_len(p0), seq_len(p0)] <- .im_penalty_matrix(basis, pen_order)
    rk <- p0 - pen_order
  } else {
    rk <- 0
  }
  # residuals of the one-class model and empirical variogram (nugget, sill)
  Qtot <- matrix(colSums(S$Q), p, p)
  beta0 <- solve(Qtot + diag(1e-8 * mean(diag(Qtot)), p), colSums(S$cc))
  r0 <- d$y - as.vector(B %*% beta0)
  vario <- .im_variogram(idx, d$t, r0)
  pr <- .im_prior(q, vy, vario, Zm)
  # no floor on s2: the priors guarantee s2 > 0
  floors <- list(G = 1e-8 * vy, s2 = 0)
  if (verbose)
    message(sprintf("  variogram: nugget = %.4g, sill = %.4g (%d pairs%s)",
                    pr$nugget, vario$sill, vario$n_pairs,
                    if (pr$nugget_source == "variogram") "" else
                      paste0("; nugget: ", pr$nugget_source)))
  ctrl <- list(inner_cap = inner_cap, tol = tol, floors = floors,
               P = P, rk = rk, penalty = penalty, prior = pr, Wc = Wc)
  score <- function(f) {
    tb <- .im_criteria_table(list(f), q, N, penalty)
    if (criterion == "BIC") tb$BIC else tb$ICL
  }
  
  if (verbose) message(sprintf("InfoMergeCluster: N = %d subjects, %d measurements, p = %d, q = %d, K0 = %d",
                               N, nrow(d), p, q, K0))
  
  # one-class model (beta0 and r0 computed above)
  v_res <- max(mean(r0^2), 1e-8 * vy)  # strictly positive starting value
  fits <- vector("list", K_max)
  lam0 <- if (penalty) 1e-3 * mean(diag(Qtot)[seq_len(p0)]) /
    (v_res * mean(diag(P)[seq_len(p0)])) else 0
  fits[[1]] <- .im_em(S, 1L, matrix(1, N, 1), diag(v_res / 2, q), v_res / 2, lam0,
                      full_cycles, ctrl)
  G <- fits[[1]]$G
  s2 <- fits[[1]]$s2
  lam <- fits[[1]]$lam
  
  pick <- function(tab) {
    crit <- if (criterion == "BIC") tab$BIC else tab$ICL
    tab$K[which.min(crit)]
  }
  
  # fit with minimum cluster size: if a cluster has a posterior size
  # < min_size, it is removed and the EM is rerun with K - 1
  valid <- function(f) all(f$pi * N >= min_size)
  try_fit <- function(K, tau0, G0, s20, lam0_, cycles) {
    f <- .im_em(S, K, tau0, G0, s20, lam0_, cycles, ctrl)
    if (valid(f)) return(f)
    keep <- f$pi * N >= min_size
    if (sum(keep) < 1) return(NULL)
    tau <- f$tau[, keep, drop = FALSE]
    tau <- tau / rowSums(tau)
    f2 <- .im_em(S, sum(keep), tau, f$G, f$s2, f$lam, cycles, ctrl)
    if (valid(f2)) f2 else NULL
  }
  store <- function(fits, f) {
    if (is.null(f)) return(fits)
    k <- f$K
    if (is.null(fits[[k]]) || score(f) < score(fits[[k]])) fits[[k]] <- f
    fits
  }
  
  K_prev <- NA_integer_
  tree <- NULL
  outer_log <- list()
  
  if (!is.null(K_fixed)) {
    # ------------------------------------------------------------------------ #
    #  K imposed by the user: direct clustering with K clusters.               #
    #  No outer loop, no selection of K, no pruning:                           #
    #  the InfoMergeCluster tree gives the initial partition, then one EM.     #
    # ------------------------------------------------------------------------ #
    K_hat <- K_max
    if (K_hat >= 2) {
      tree <- .im_tree(S, G, s2, lam, K0, K_hat, a_dir, nstart, floors, ctrl)
      lab <- tree$parts[[K_hat]]
      if (is.null(lab)) lab <- sample(rep_len(seq_len(K_hat), N))   # fallback
      fits[[K_hat]] <- .im_em(S, K_hat, .im_onehot(lab, K_hat, 0.01),
                              tree$G_within, tree$s2_within, lam, full_cycles, ctrl)
      small <- fits[[K_hat]]$pi * N < min_size
      if (any(small))
        warning(sprintf("Imposed K = %d: %d cluster(s) with posterior size < min_size = %g.",
                        K_hat, sum(small), min_size), call. = FALSE)
    }
    if (verbose) message(sprintf("  clustering with K = %d (imposed)", K_hat))
    tab <- .im_criteria_table(fits, q, N, penalty)
  } else for (o in seq_len(n_outer)) {
    tree <- .im_tree(S, G, s2, lam, K0, K_max, a_dir, nstart, floors, ctrl)
    G_start <- if (o == 1) tree$G_within else G
    s2_start <- if (o == 1) tree$s2_within else s2
    
    if (K_max >= 2) {
      for (K in 2:K_max) {
        lab <- tree$parts[[K]]
        if (is.null(lab)) next
        fits <- store(fits, try_fit(K, .im_onehot(lab, K, 0.01), G_start, s2_start,
                                    lam, quick_cycles))
      }
    }
    K_hat <- pick(.im_criteria_table(fits, q, N, penalty))
    
    for (K in unique(pmin(pmax(K_hat + (-1:1), 1L), K_max))) {
      if (is.null(fits[[K]])) next
      fits <- store(fits, try_fit(K, fits[[K]]$tau, fits[[K]]$G, fits[[K]]$s2,
                                  fits[[K]]$lam, full_cycles))
    }
    tab <- .im_criteria_table(fits, q, N, penalty)
    K_hat <- pick(tab)
    
    outer_log[[o]] <- data.frame(iteration = o, K_hat = K_hat, K_tree = tree$K_tree)
    if (verbose) message(sprintf("  outer iteration %d: K (%s) = %d   [tree K = %d]",
                                 o, criterion, K_hat, tree$K_tree))
    if (!is.na(K_prev) && K_hat == K_prev) break
    K_prev <- K_hat
    G <- fits[[K_hat]]$G
    s2 <- fits[[K_hat]]$s2
    lam <- fits[[K_hat]]$lam
  }
  
  best <- fits[[K_hat]]
  ord <- order(best$pi, decreasing = TRUE)
  beta <- best$beta[seq_len(p0), ord, drop = FALSE]
  gamma <- NULL
  if (r > 0) {
    gamma <- best$beta[p0 + seq_len(r), ord, drop = FALSE]
    dimnames(gamma) <- list(sub("^cov_", "", cov_names), paste0("C", seq_len(K_hat)))
  }
  alpha <- NULL
  if (!is.null(best$alpha)) {
    alpha <- best$alpha[, ord, drop = FALSE]
    alpha <- alpha - alpha[, 1]                 # reference: cluster C1
    dimnames(alpha) <- list(c("(Intercept)", sub("^cov_", "", cov_names)),
                            paste0("C", seq_len(K_hat)))
  }
  tau <- best$tau[, ord, drop = FALSE]
  labs <- paste0("C", seq_len(K_hat))
  dimnames(tau) <- list(as.character(ids), labs)
  colnames(beta) <- labs
  
  out <- list(
    K         = K_hat,
    K_fixed   = !is.null(K_fixed),
    re        = re,
    criterion = criterion,
    cluster   = stats::setNames(max.col(tau, ties.method = "first"), as.character(ids)),
    posterior = tau,
    pi        = stats::setNames(best$pi[ord], labs),
    beta      = beta,
    gamma     = gamma,
    alpha     = alpha,
    weights   = weights,
    covariates = if (r > 0) list(names = cov_names, center = cov_center) else NULL,
    Vbeta     = stats::setNames(best$Vbeta[ord], labs),
    G         = best$G,
    sigma2    = best$s2,
    lambda    = best$lam,
    edf       = best$edf,
    loglik    = best$ll,
    log_prior = best$log_prior,
    prior     = pr,
    variogram = vario,
    table     = tab,
    K_tree    = tree$K_tree,
    icl_path  = tree$icl_path,
    outer     = do.call(rbind, outer_log),
    basis     = basis,
    fits      = fits,
    N         = N,
    n_obs     = nrow(d),
    var_y     = vy
  )
  class(out) <- "infomergecluster"
  out
}
