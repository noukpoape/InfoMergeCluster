# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi NOUKPOAPE <noukpoape@gmail.com>

#' InfoMergeCluster: clustering of sparse longitudinal data
#'
#' Model-based clustering of sparse and irregular longitudinal data, with
#' penalized B-spline mean trajectories, a random effect, conjugate priors on
#' the variance structure calibrated on the empirical variogram, and
#' automatic (BIC/ICL) or user-imposed choice of the number of clusters.
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{infomergecluster}}: fit the model;
#'   \item \code{\link{plot.infomergecluster}}, \code{\link{predict.infomergecluster}}:
#'     mean curves and envelopes;
#'   \item \code{\link{fitted_curves}}: predicted individual trajectories;
#'   \item \code{\link{infomergecluster_stability}}: stability across seeds;
#'   \item \code{\link{f.data.generation}}, \code{\link{plot.data}}:
#'     simulated data (scenarios 1 to 4) and their plot.
#' }
#'
#' @keywords internal
#' @importFrom stats predict rnorm runif rbinom
#' @importFrom graphics lines legend
#' @importFrom grDevices hcl.colors
"_PACKAGE"
