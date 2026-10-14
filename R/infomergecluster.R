# Part of the InfoMergeCluster package.
# Authors: Komlan Midodzi Noukpoape <noukpoape@gmail.com>, Cristian Preda

#' Cluster sparse longitudinal data
#'
#' Fits a mixture of linear mixed models with penalized B-spline mean
#' curves to sparse longitudinal data, and either selects the number of
#' clusters automatically (BIC or ICL) or uses the number of clusters
#' supplied in \code{K}.
#'
#' Conditionally on its cluster \eqn{z_i = k}, the measurements of subject
#' \eqn{i} follow
#' \deqn{y_i = B_i \beta_k + X_i \gamma_k + Z_i u_i + \varepsilon_i, \quad
#'   u_i \sim N_q(0, G), \quad \varepsilon_i \sim N(0, \sigma^2 I),}
#' with \eqn{P(z_i = k) = \pi_k}. \eqn{B_i} is a B-spline basis evaluated at
#' the times of subject \eqn{i}, \eqn{\beta_k} the coefficients of the mean
#' curve of cluster \eqn{k}, \eqn{X_i} the (optional, centered) covariates
#' with cluster-specific effects \eqn{\gamma_k}, and \eqn{Z_i} the
#' random-effect design (intercept, linear or spline). The mean curves are
#' penalized B-splines (O'Sullivan penalty on the \code{pen_order}-th
#' derivative), with a smoothing parameter estimated within the EM algorithm
#' (Fellner-Schall update); the covariate effects are not penalized.
#'
#' Conjugate priors are always placed on the variance structure,
#' \eqn{G \sim IW_q(q, \Psi_0)} and \eqn{\sigma^2 \sim IG(1, b_\sigma)}.
#' Their hyperparameters are calibrated on the empirical variogram of the
#' residuals of the one-class model: the prior mode of \eqn{\sigma^2} is the
#' nugget, and the prior mode of \eqn{G} attributes half of (sill - nugget)
#' to the random effect. The estimator is the posterior mode; the residual
#' variance cannot collapse to zero, even for data without measurement
#' error. The reported log-likelihood, BIC and ICL do not include the prior.
#'
#' With \code{weights = "covariates"}, the cluster proportions depend on the
#' subject's covariates (concomitant variables),
#' \eqn{\pi_k(x_i) = \exp(w_i'\alpha_k) / \sum_l \exp(w_i'\alpha_l)}, with
#' \eqn{w_i = (1, \bar{x}_i)} (subject means of the centered covariates).
#'
#' When \code{K} is \code{NULL}, an InfoMergeCluster tree (greedy merges of
#' micro-clusters based on the exact integrated completed likelihood)
#' initializes the EM for every \eqn{K \le} \code{K_max}, and an outer loop
#' refines the variance parameters until the selected \eqn{K} is stable.
#' Clusters smaller than \code{min_size} are pruned. When \code{K} is given,
#' clustering is done directly with \code{K} clusters (one tree for the
#' initialization, then a single EM), with no outer loop and no pruning.
#'
#' @param data data, in one of the following forms (columns are identified
#'   by position). Without covariates (\code{n_cov = 0}): a long data.frame
#'   (1st column = subject identifier, second-to-last = time, last =
#'   response; intermediate columns are ignored), or a list of data.frames
#'   or matrices, one element per subject (second-to-last column = time,
#'   last = response; the identifier is the element name, otherwise its
#'   rank). With covariates, see \code{n_cov}.
#' @param K_max maximum number of clusters examined (ignored if \code{K} is given).
#' @param K imposed number of clusters (integer \eqn{\ge 1}). \code{NULL}
#'   (default): \code{K} is chosen automatically between 1 and \code{K_max}.
#'   If \code{K} is given, it is not estimated (\code{n_outer} and
#'   \code{quick_cycles} are then ignored, and no cluster is pruned).
#' @param df number of B-spline functions (p).
#' @param degree spline degree (3 = cubic).
#' @param penalty logical; penalize the integral of the squared derivative of
#'   order \code{pen_order}, with an estimated smoothing parameter
#'   (\code{TRUE} recommended).
#' @param pen_order order m of the penalized derivative
#'   (\eqn{1 \le m \le} \code{degree}; 2 = curvature).
#' @param re random-effect structure: \code{"linear"} (random intercept and
#'   slope, default), \code{"intercept"} (random intercept) or
#'   \code{"spline"} (random effect on the B-spline basis).
#' @param n_cov number of covariates (0 by default) acting on the mean
#'   trajectories with cluster-specific effects \eqn{\gamma_k}. Columns are
#'   identified by position: for a long data.frame, 1st column = identifier,
#'   2nd = time, 3rd = response and columns 4 to \code{3 + n_cov} =
#'   covariates; for a list (one element per subject), 1st column = time,
#'   2nd = response and columns 3 to \code{2 + n_cov} = covariates. One value
#'   per measurement (constant within a subject or time-varying). Factors are
#'   coded by dummy variables and all covariates are centered, so that the
#'   mean curves are those of a subject with mean covariate values.
#' @param weights \code{"constant"} (cluster proportions \eqn{\pi_k}, default)
#'   or \code{"covariates"} (proportions depending on the subject's covariates
#'   through a multinomial logistic model; requires \code{n_cov > 0}).
#' @param min_size minimum cluster size (posterior size \eqn{N \pi_k}).
#' @param K0 number of micro-clusters of the tree (\code{NULL} = automatic).
#' @param criterion \code{"BIC"} (target: mixture components) or
#'   \code{"ICL"} (well-separated classes).
#' @param a_dir Dirichlet prior of the exact ICL; only affects the tree used
#'   to initialize the EM.
#' @param n_outer maximum number of iterations of the outer loop.
#' @param quick_cycles SQUAREM cycles for the quick fit of each K.
#' @param full_cycles SQUAREM cycles for the full fits.
#' @param inner_cap (E, M-beta) steps per evaluation of the EM map.
#' @param tol relative tolerance on the log-likelihood.
#' @param nstart number of starts of the over-segmentation k-means.
#' @param seed optional random seed.
#' @param keep_data keep the data in the result (needed by \code{plot} and
#'   \code{\link{fitted_curves}}).
#' @param verbose print progress messages.
#'
#' @return An object of class \code{"infomergecluster"}, a list containing
#'   in particular:
#'   \item{K}{number of clusters.}
#'   \item{K_fixed}{\code{TRUE} if \code{K} was imposed by the user.}
#'   \item{cluster}{MAP cluster of each subject (named by subject identifier).}
#'   \item{posterior}{matrix of posterior membership probabilities.}
#'   \item{pi}{cluster proportions.}
#'   \item{beta}{B-spline coefficients of the mean curves (one column per cluster).}
#'   \item{gamma}{cluster-specific covariate effects (\code{NULL} without covariates).}
#'   \item{alpha}{coefficients of the multinomial logistic model of the
#'     weights, cluster C1 as reference (\code{NULL} unless
#'     \code{weights = "covariates"}).}
#'   \item{G, sigma2, lambda}{random-effect covariance, residual variance and
#'     smoothing parameter.}
#'   \item{loglik, edf}{log-likelihood and effective degrees of freedom.}
#'   \item{table}{criteria (log-likelihood, BIC, ICL, entropy) of every fitted K.}
#'   \item{re}{random-effect structure used.}
#'   \item{prior, variogram}{hyperparameters of the priors and empirical variogram.}
#'   \item{basis, data}{B-spline basis and data in long format.}
#'
#' @references
#' Biernacki, C., Celeux, G. and Govaert, G. (2000). Assessing a mixture model
#' for clustering with the integrated completed likelihood. IEEE Transactions
#' on Pattern Analysis and Machine Intelligence, 22(7), 719--725.
#'
#' Dayton, C. M. and Macready, G. B. (1988). Concomitant-variable latent-class
#' models. Journal of the American Statistical Association, 83(401), 173--178.
#'
#' Diggle, P. J. (1988). An approach to the analysis of repeated measurements.
#' Biometrics, 44(4), 959--971.
#'
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with B-splines
#' and penalties. Statistical Science, 11(2), 89--121.
#'
#' James, G. M. and Sugar, C. A. (2003). Clustering for sparsely sampled
#' functional data. Journal of the American Statistical Association, 98(462),
#' 397--408.
#'
#' O'Sullivan, F. (1986). A statistical perspective on ill-posed inverse
#' problems. Statistical Science, 1(4), 502--518.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. The Annals of
#' Statistics, 6(2), 461--464.
#'
#' Varadhan, R. and Roland, C. (2008). Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm. Scandinavian Journal
#' of Statistics, 35(2), 335--353.
#'
#' Wood, S. N. and Fasiolo, M. (2017). A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. Biometrics, 73(4), 1071--1081.
#'
#' @seealso \code{\link{plot.infomergecluster}},
#'   \code{\link{predict.infomergecluster}}, \code{\link{fitted_curves}},
#'   \code{\link{infomergecluster_stability}}, \code{\link{f.data.generation}}
#'
#' @examples
#' # scenario 1 with its reference settings: 2 groups of 43 and 57 curves
#' set.seed(1)
#' gen <- f.data.generation(senario = 1, n.k = c(43, 57), n.cts = c(1.8, 2.4), n.min = 2)
#' plot.data(gen$liste, legend = TRUE)
#'
#' # number of clusters imposed
#' fit2 <- infomergecluster(gen$df, K = 2, re = "intercept", verbose = FALSE)
#' fit2
#' table(fit2$cluster, attr(gen$liste, "groupe"))
#'
#' \donttest{
#' # number of clusters selected by BIC
#' fit <- infomergecluster(gen$df, K_max = 5, re = "intercept", verbose = FALSE)
#' plot(fit, envelope = TRUE)
#'
#' # with a covariate acting as a level shift (requires the fda package):
#' # 3 clusters, y = f_k(t) + gamma_k z + error, gamma = (1, 0, 1), f_3 = f_1
#' if (requireNamespace("fda", quietly = TRUE)) {
#'   set.seed(1234)
#'   m <- 5
#'   b <- fda::create.bspline.basis(c(0, 1), nbasis = m)
#'   b1 <- sample(0:10, m, replace = TRUE)
#'   b2 <- sample(0:10, m, replace = TRUE)
#'   f_list <- list(fda::fd(b1, b), fda::fd(b2, b), fda::fd(b1, b))
#'   n1 <- 30; n2 <- 40; n3 <- 30
#'   grp <- rep(1:3, c(n1, n2, n3))
#'   gam <- c(1, 0, 1)[grp]
#'   z_all <- c(rnorm(n1, 1, 1), rnorm(n2, 0, 1), rnorm(n3, -1, 1))
#'   # one matrix per subject: time, response, covariate
#'   fdata.cov <- lapply(seq_along(grp), function(i) {
#'     x <- sort(runif(3 + rbinom(1, 8, 0.5)))
#'     y <- fda::eval.fd(f_list[[grp[i]]], x) + gam[i] * z_all[i] +
#'       rnorm(length(x), 0, 0.5)
#'     cbind(x = x, y = y, z = rep(z_all[i], length(x)))
#'   })
#'   plot.data(fdata.cov, n.k = c(n1, n2, n3))
#'   fit_2 <- infomergecluster(data = fdata.cov, n_cov = 1, re = "intercept",
#'                             df = 20, K_max = 10)
#'   fit_2$gamma
#'   plot(fit_2, envelope = TRUE, level = 0.95, legend = TRUE)
#' }
#' }
#' @export
infomergecluster <- function(data, K_max = 8, K = NULL, df = 15, degree = 3,
                             penalty = TRUE, pen_order = 2,
                             re = c("linear", "intercept", "spline"),
                             n_cov = 0, weights = c("constant", "covariates"),
                             min_size = 2, K0 = NULL,
                             criterion = c("BIC", "ICL"), a_dir = 0.5,
                             n_outer = 3, quick_cycles = 4, full_cycles = 15,
                             inner_cap = 3, tol = 1e-5, nstart = 5, seed = NULL,
                             keep_data = TRUE, verbose = TRUE) {
  
  re <- match.arg(re)
  weights <- match.arg(weights)
  criterion <- match.arg(criterion)
  t_start <- proc.time()[["elapsed"]]
  if (!is.numeric(n_cov) || length(n_cov) != 1 || !is.finite(n_cov) || n_cov < 0 ||
      n_cov != round(n_cov))
    stop("'n_cov' must be an integer >= 0.")
  n_cov <- as.integer(n_cov)
  
  d <- .im_as_long(data, n_cov)
  cov_raw <- if (n_cov > 0) d[paste0("z", seq_len(n_cov))] else NULL
  d <- d[, c("id", "t", "y"), drop = FALSE]
  n_raw <- nrow(d)
  ok <- stats::complete.cases(d)
  if (!is.null(cov_raw)) ok <- ok & stats::complete.cases(cov_raw)
  d <- d[ok, , drop = FALSE]
  if (verbose && nrow(d) < n_raw)
    message(sprintf("%d incomplete or non-numeric measurement(s) ignored.", n_raw - nrow(d)))
  cov <- .im_covariate_design(cov_raw, ok)
  if (weights == "covariates" && length(cov$names) == 0)
    stop("weights = \"covariates\" requires n_cov > 0.")
  if (length(cov$names)) d <- cbind(d, cov$X)
  d <- d[order(d$id, d$t), , drop = FALSE]
  rownames(d) <- NULL
  N_subj <- length(unique(d$id))
  if (N_subj < 10) stop("Too few subjects (N < 10).")
  if (!is.null(K)) {
    if (!is.numeric(K) || length(K) != 1 || !is.finite(K) || K < 1 || K != round(K))
      stop("'K' must be NULL or an integer >= 1.")
    K <- as.integer(K)
    if (K > floor(N_subj / 2))
      stop(sprintf("'K' = %d is too large for N = %d subjects (K <= %d).",
                   K, N_subj, floor(N_subj / 2)))
  }
  
  args <- list(K_max = K_max, K_fixed = K, df = df, degree = degree, penalty = penalty,
               pen_order = pen_order, cov_names = cov$names, cov_center = cov$center,
               weights = weights, min_size = min_size, K0 = K0,
               criterion = criterion, a_dir = a_dir, n_outer = n_outer,
               quick_cycles = quick_cycles, full_cycles = full_cycles,
               inner_cap = inner_cap, tol = tol, nstart = nstart, verbose = verbose)
  
  if (!is.null(seed)) set.seed(seed)
  if (verbose) message(sprintf("--- random effect: %s", re))
  out <- do.call(.im_fit_one, c(list(d = d, re = re), args))
  out$time_sec <- proc.time()[["elapsed"]] - t_start
  out$data <- if (keep_data) d else NULL
  out
}

# conversion of the data to the internal long format (id, t, y), by column
# position:
#   - long data.frame (or matrix): 1st column = subject identifier,
#     second-to-last = time, last = response;
#   - list of data.frames or matrices (one element per subject):
#     second-to-last column = time, last = response; the identifier is the
#     element name if all names are non-empty and distinct, otherwise its rank.
.im_num <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  suppressWarnings(as.numeric(v))
}

#   n_cov > 0 (columns identified by position):
#     - long data.frame: 1st = identifier, 2nd = time, 3rd = response,
#       columns 4..(3 + n_cov) = covariates;
#     - list: 1st column = time, 2nd = response, columns 3..(2 + n_cov) =
#       covariates.
#   The covariates are returned as columns z1, ..., z<n_cov>, with their type.
.im_as_long <- function(data, n_cov = 0L) {
  zn <- if (n_cov > 0) paste0("z", seq_len(n_cov)) else character(0)
  grab <- function(e, cols) {
    cv <- e[cols]
    names(cv) <- zn
    cv
  }
  if (is.data.frame(data) || is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    nc <- ncol(data)
    if (nc < 3)
      stop("A long data.frame must have at least 3 columns: identifier (1st), ",
           "time (second-to-last) and response (last).")
    id <- data[[1]]
    if (is.factor(id)) id <- as.character(id)
    if (n_cov > 0) {
      if (nc < 3L + n_cov)
        stop(sprintf("With n_cov = %d, a long data.frame must have at least %d columns: identifier, time, response and the covariates.", n_cov, 3L + n_cov))
      out <- data.frame(id = id, t = .im_num(data[[2]]), y = .im_num(data[[3]]),
                        stringsAsFactors = FALSE)
      return(cbind(out, grab(data, 3L + seq_len(n_cov))))
    }
    return(data.frame(id = id, t = .im_num(data[[nc - 1]]), y = .im_num(data[[nc]]),
                      stringsAsFactors = FALSE))
  }
  if (is.list(data)) {
    if (length(data) == 0) stop("The list 'data' is empty.")
    ok <- vapply(data, function(e) is.data.frame(e) || is.matrix(e), logical(1))
    if (!all(ok))
      stop("Each element of the list must be a data.frame or a matrix ",
           "(offending elements: ", paste(which(!ok), collapse = ", "), ").")
    nm <- names(data)
    ids <- if (!is.null(nm) && all(nzchar(nm)) && !anyDuplicated(nm)) nm else seq_along(data)
    parts <- lapply(seq_along(data), function(i) {
      e <- as.data.frame(data[[i]], stringsAsFactors = FALSE)
      nc <- ncol(e)
      if (nc < 2)
        stop(sprintf("Element %d must have at least 2 columns: time (second-to-last) and response (last).", i))
      if (nrow(e) == 0) return(NULL)
      if (n_cov > 0) {
        if (nc < 2L + n_cov)
          stop(sprintf("With n_cov = %d, element %d must have at least %d columns: time, response and the covariates.", n_cov, i, 2L + n_cov))
        out_i <- data.frame(id = ids[i], t = .im_num(e[[1]]), y = .im_num(e[[2]]),
                            stringsAsFactors = FALSE)
        return(cbind(out_i, grab(e, 2L + seq_len(n_cov))))
      }
      data.frame(id = ids[i], t = .im_num(e[[nc - 1]]), y = .im_num(e[[nc]]),
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, parts)
    if (is.null(out)) stop("No measurements in 'data'.")
    return(out)
  }
  stop("'data' must be a data.frame, a list of data.frames or a list of matrices.")
}

# numeric design of the covariates (dummy coding, centering) on the kept rows
.im_covariate_design <- function(cov_raw, ok) {
  if (is.null(cov_raw)) return(list(X = NULL, names = character(0), center = numeric(0)))
  cv <- cov_raw[ok, , drop = FALSE]
  cv[] <- lapply(cv, function(v) if (is.character(v)) factor(v) else v)
  cv <- droplevels(cv)
  X <- stats::model.matrix(~ ., data = cv)[, -1, drop = FALSE]
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  keep <- apply(X, 2, function(v) stats::var(v) > 0)
  if (any(!keep))
    warning("Covariate column(s) without variation removed: ",
            paste(colnames(X)[!keep], collapse = ", "), call. = FALSE)
  X <- X[, keep, drop = FALSE]
  if (ncol(X) == 0) return(list(X = NULL, names = character(0), center = numeric(0)))
  center <- colMeans(X)
  X <- sweep(X, 2, center)
  nm <- paste0("cov_", colnames(X))
  colnames(X) <- nm
  list(X = as.data.frame(X), names = nm, center = stats::setNames(center, nm))
}

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
