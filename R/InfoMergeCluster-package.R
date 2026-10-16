###############################################################################
#  InfoMergeCluster (v1): clustering of sparse longitudinal data
#                   with automatic selection of the number of clusters
#
#  Model    y_i | z_i = k  ~  N( B_i beta_k , V_i ),   V_i = Z_i G Z_i' + s2 I
#           B_i : B-spline basis (mean trajectory of the cluster)
#           Z_i : random-effect design, one of:
#                 "intercept" (1), "linear" (1, t), "spline" (Z_i = B_i)
#           G   : random-effect covariance (shared by all clusters)
#
#  Version 1
#    Conjugate priors on the variance structure,
#        G ~ IW_q(nu0, Psi0)   and   s2 ~ IG(a_sigma, b_sigma),
#    (section "Priors and penalized likelihood" of the document).
#    The estimator is the posterior mode, with z and u integrated out, in
#    the (G, s2) parameterization. The M-step for (G, s2) is
#        G  <- (N Sigma_G + Psi0) / (N + nu0 + q + 1),
#        s2 <- (sum(n_i) s2_EM + 2 b_sigma) / (sum(n_i) + 2 a_sigma + 2),
#    with no floor on s2: the M-step gives
#        s2 >= 2 b_sigma / (sum(n_i) + 2 a_sigma + 2) > 0  and  G >= Psi0 / (N + nu0 + q + 1),
#    and the posterior mode exists without floors (proposition "existence of
#    the posterior mode without floors" of the document).
#    The objective tracked by SQUAREM includes log p0(G, s2); the reported
#    log-likelihood (loglik, BIC, ICL) remains l_N, without the prior,
#    which does not enter the number of parameters.
#    The priors are always active and their hyperparameters cannot be set
#    by the user (nu0 = q, a_sigma = 1, g_share = 0.5; nugget and sill
#    estimated from the variogram).
#    Hyperparameters computed from the empirical variogram (Diggle,
#    1988) of the residuals r of the one-class model: for each pair of
#    measurements of the same subject, v = (r_ij - r_ik)^2 / 2 and
#    u = |t_ij - t_ik|. Under the model, the smooth components (class
#    deviations, random effect) vanish as u -> 0, so that the nugget
#    gamma(0+) = s2. It is estimated by the intercept of a regression of v
#    on u restricted to the smallest lags. The sill is var(r). Then
#        prior mode of s2 = nugget              (b_sigma = (a_sigma + 1) nugget)
#        prior mode of G  = g * I_q, with g such that the mean variance
#            of Z(t)'u at the observed times equals g_share (sill - nugget)
#            (Psi0 = (nu0 + q + 1) g I_q).
#    Sill minus nugget also contains the heterogeneity between classes:
#    g_share = 0.5 sets the share attributed a priori to the random effect.
#    The random-effect structure (argument re) is chosen by the user; a
#    single fit is performed.
#    Covariates (argument n_cov, columns identified by position): model
#        y_i | z_i = k = B_i beta_k + X_i gamma_k + Z_i u_i + e_i,
#    with cluster-specific effects gamma_k. The basis is augmented,
#    [B_i | X_i], the penalty is diag(P, 0) (gamma_k is not smoothed), and
#    the information statistics, the EM and the exact merge cost are
#    unchanged in dimension p + r. Covariates are centered (on all
#    measurements), so that beta_k describes the mean trajectory at the
#    mean covariate values. Covariates are given in the data, one value per
#    measurement, and may be constant within a subject or time-varying.
#    Argument weights: "constant" (pi_k) or "covariates" (concomitant
#    variables: pi_k(x_i) = exp(w_i' alpha_k) / sum_l exp(w_i' alpha_l),
#    w_i = (1, subject mean of the centered covariates)); the M-step of
#    alpha is a weighted multinomial logistic regression (Newton steps,
#    small ridge on the slopes).
#
#  Dependencies: base packages (stats, splines, graphics, grDevices).
###############################################################################

#' InfoMergeCluster: clustering of sparse longitudinal data
#'
#' Clustering of sparse longitudinal trajectories with a mixture of mixed
#' models whose mean curves are penalised B-splines, optionally with
#' covariates, with automatic choice of the number of clusters (exact-ICL
#' merge tree, SQUAREM-accelerated EM, BIC/ICL) or a number of clusters
#' fixed by the user.
#'
#' Main functions: \code{\link{infomergecluster}},
#' \code{\link{infomergecluster_stability}},
#' \code{\link{predict.infomergecluster}},
#' \code{\link{plot.infomergecluster}}, \code{\link{fitted_curves}}.
#' Simulation: \code{\link{f.data.generation}}.
#'
#' @keywords internal
"_PACKAGE"
