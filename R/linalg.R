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
