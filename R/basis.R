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
