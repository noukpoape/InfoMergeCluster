# =========================================================================== #
#  Main function                                                             #
# =========================================================================== #

#' Clustering of sparse longitudinal data (InfoMergeCluster)
#'
#' Fits a mixture of linear mixed models whose cluster mean curves are
#' penalised B-splines (O'Sullivan), optionally with covariates, and selects
#' the number of clusters automatically, unless it is fixed by the user
#' through \code{K}. A greedy merge tree based on the exact integrated
#' completed likelihood (ICL), built on an over-segmentation into
#' micro-clusters, initialises an EM algorithm (accelerated by SQUAREM) for
#' each number of clusters; the retained number minimises BIC (or ICL).
#'
#' \strong{Model.} \eqn{y_i | z_i = k \sim N(B_i \beta_k + X_i \gamma_k, V_i)},
#' with \eqn{V_i = Z_i G Z_i' + \sigma^2 I}, where \eqn{B_i} is the B-spline
#' basis, \eqn{X_i} the (centred) covariates, \eqn{\gamma_k} their
#' cluster-specific effects, \eqn{Z_i} the random-effect design and \eqn{G}
#' its covariance (shared by all clusters).
#'
#' \strong{Data layout.} Columns are identified by position, without names.
#' Without covariates (\code{n_cov = 0}):
#' \itemize{
#'   \item long data frame: first column = subject identifier, second to
#'     last = time, last = response (columns in between, e.g. a group label,
#'     are ignored);
#'   \item list of data frames or matrices (one element per subject): second
#'     to last column = time, last = response; the identifier is the element
#'     name when all names are non-empty and distinct, otherwise its index.
#' }
#' With covariates (\code{n_cov > 0}):
#' \itemize{
#'   \item long data frame: column 1 = identifier, 2 = time, 3 = response,
#'     columns 4 to \code{3 + n_cov} = covariates;
#'   \item list: column 1 = time, 2 = response, columns 3 to
#'     \code{2 + n_cov} = covariates.
#' }
#' Covariates take one value per measurement (constant within a subject or
#' time-varying). They are named \code{z1}, ..., factors and character
#' columns are coded by dummy variables, and all covariates are centred, so
#' that \eqn{\beta_k} describes the mean trajectory at the mean covariate
#' values. With \code{weights = "covariates"}, the cluster proportions depend
#' on the subject's mean covariates through a multinomial logistic model
#' (concomitant variables).
#'
#' \strong{Priors.} The variance structure always has conjugate priors,
#' \eqn{G \sim IW_q(q, \Psi_0)} and \eqn{\sigma^2 \sim IG(1, b_\sigma)}, so
#' that the estimator is the posterior mode and no lower bound on
#' \eqn{\sigma^2} is needed. The hyperparameters are not user-settable: they
#' are calibrated on the empirical variogram (Diggle, 1988) of the residuals
#' of the one-class model. The prior mode of \eqn{\sigma^2} is the nugget,
#' and the prior mode of \eqn{G} is \eqn{g I_q}, where \eqn{g} assigns half
#' of (sill - nugget) to the random effect. The reported log-likelihood
#' (\code{loglik}, BIC, ICL) does not include the prior.
#'
#' @param data data: long data frame, list of data frames or list of
#'   matrices (see Details).
#' @param K_max maximum number of clusters considered when the number of
#'   clusters is selected automatically (ignored when \code{K} is given).
#' @param K number of clusters. \code{NULL} (default): selected
#'   automatically by BIC (or ICL) among 1, ..., \code{K_max}. A positive
#'   integer: clustering is done directly with \code{K} clusters and K is
#'   not estimated (no outer loop, \code{n_outer} and \code{quick_cycles}
#'   are ignored, no pruning; a warning is issued if a cluster is smaller
#'   than \code{min_size}). Must not exceed \eqn{N / 2}.
#' @param df number of B-spline basis functions (\eqn{p}).
#' @param degree spline degree (3 = cubic).
#' @param penalty if \code{TRUE} (recommended), integrated squared
#'   derivative penalty of order \code{pen_order}, with the smoothing
#'   parameter estimated by Fellner-Schall updates.
#' @param pen_order order \eqn{m} of the penalised derivative
#'   (\code{1 <= pen_order <= degree}; 2 = curvature).
#' @param re random-effect structure: \code{"linear"} (1, t),
#'   \code{"intercept"} (1) or \code{"spline"} (\eqn{Z_i = B_i}).
#' @param n_cov number of covariates (0 by default) acting on the mean
#'   trajectories with cluster-specific effects; see Details for the column
#'   layout.
#' @param weights \code{"constant"} (cluster proportions \eqn{\pi_k},
#'   default) or \code{"covariates"} (proportions depending on the subject's
#'   covariates through a multinomial logistic model; requires
#'   \code{n_cov > 0}).
#' @param min_size minimum cluster size (posterior count \eqn{N \pi_k}); in
#'   automatic selection, a smaller cluster is pruned and the model refitted
#'   with \eqn{K - 1} clusters.
#' @param K0 number of micro-clusters in the over-segmentation (\code{NULL}
#'   = automatic).
#' @param criterion \code{"BIC"} (targets mixture components) or
#'   \code{"ICL"} (targets well-separated classes).
#' @param a_dir Dirichlet parameter of the exact ICL; only affects the tree
#'   used to initialise EM.
#' @param n_outer maximum number of outer iterations.
#' @param quick_cycles SQUAREM cycles of the quick fit of each K.
#' @param full_cycles SQUAREM cycles of the full fit of K-1, K and K+1.
#' @param inner_cap (E, M-beta) steps per evaluation of the EM map.
#' @param tol relative tolerance on the log-likelihood.
#' @param nstart number of k-means starts for the over-segmentation.
#' @param seed optional random seed.
#' @param keep_data keep the data in the returned object (needed by
#'   \code{plot(..., show_data = TRUE)} and \code{\link{fitted_curves}}).
#' @param verbose print progress messages.
#'
#' @return An object of class \code{"infomergecluster"}, a list containing in
#'   particular:
#' \item{K}{number of clusters (selected, or fixed by the user).}
#' \item{K_fixed}{\code{TRUE} when \code{K} was given by the user.}
#' \item{cluster}{MAP assignment of each subject (vector named by
#'   identifier).}
#' \item{posterior}{posterior probabilities (subjects x clusters).}
#' \item{pi}{mean cluster proportions (decreasing order).}
#' \item{beta}{B-spline coefficients of the mean curves (one column per
#'   cluster).}
#' \item{gamma}{cluster-specific covariate effects (covariates x clusters),
#'   or \code{NULL}.}
#' \item{alpha}{coefficients of the multinomial logistic model of the
#'   proportions (reference cluster C1) when
#'   \code{weights = "covariates"}, otherwise \code{NULL}.}
#' \item{covariates}{names and centring values of the covariates, or
#'   \code{NULL}.}
#' \item{Vbeta}{posterior covariances of the coefficients.}
#' \item{G, sigma2, lambda}{random-effect covariance, residual variance,
#'   smoothing parameter.}
#' \item{edf}{total effective degrees of freedom.}
#' \item{loglik, log_prior}{log-likelihood and log prior density at the
#'   mode.}
#' \item{table}{criteria (loglik, BIC, ICL, entropy, ...) for each fitted
#'   K.}
#' \item{prior, variogram}{hyperparameters used and empirical variogram.}
#' \item{basis, fits, data}{B-spline basis, fits for each K, data in long
#'   format.}
#'
#' @references
#' Diggle, P. J. (1988). An approach to the analysis of repeated
#' measurements. \emph{Biometrics}, 44, 959--971.
#'
#' Varadhan, R. and Roland, C. (2008). Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm.
#' \emph{Scandinavian Journal of Statistics}, 35, 335--353.
#'
#' @seealso \code{\link{predict.infomergecluster}},
#'   \code{\link{plot.infomergecluster}}, \code{\link{fitted_curves}},
#'   \code{\link{infomergecluster_stability}}.
#'
#' @examples
#' set.seed(1)
#' g <- f.data.generation(senario = 1, n.k = c(30, 40), n.cts = c(0, 2))
#'
#' # number of clusters selected automatically (list format)
#' fit <- infomergecluster(g$liste, K_max = 4, verbose = FALSE)
#' fit
#' adjusted_rand(fit$cluster, attr(g$liste, "groupe"))
#' plot(fit, envelope = TRUE)
#'
#' # number of clusters fixed by the user (long data frame)
#' fit2 <- infomergecluster(g$df, K = 2, verbose = FALSE)
#' fit2$K
#'
#' # with a covariate: columns id, time, response, covariate
#' d <- g$df[, c("id", "x", "y")]
#' age <- stats::rnorm(70, 50, 10)
#' d$age <- age[d$id]
#' fit3 <- infomergecluster(d, K = 2, n_cov = 1, verbose = FALSE)
#' fit3$gamma
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
