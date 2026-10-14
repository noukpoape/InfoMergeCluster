# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi Noukpoape <noukpoape@gmail.com>

#' InfoMergeCluster: clustering of sparse longitudinal data
#'
#' Model-based clustering of sparse and irregular longitudinal data with a
#' mixture of linear mixed models: penalized B-spline mean trajectories, a
#' random effect, conjugate priors on the variance structure calibrated on
#' the empirical variogram, optional covariates (on the trajectories and on
#' the cluster proportions), and automatic (BIC/ICL) or user-imposed choice
#' of the number of clusters.
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{infomergecluster}}: fit the model;
#'   \item \code{\link{plot.infomergecluster}}, \code{\link{predict.infomergecluster}}:
#'     mean curves and envelopes;
#'   \item \code{\link{fitted_curves}}: predicted individual trajectories;
#'   \item \code{\link{infomergecluster_stability}}: stability across seeds;
#'   \item \code{\link{f.data.generation}}, \code{\link{plot.data}}:
#'     simulated data (scenarios 1 to 4) and their plot;
#'   \item \code{\link{infomergecluster_citation}}: how to cite the package.
#' }
#'
#' @section Citation:
#' To cite InfoMergeCluster in publications, use
#' \code{infomergecluster_citation()}, which prints the BibTeX entry:
#' \preformatted{
#' @@manual{noukpoape2026infomergecluster,
#'   title={InfoMergeCluster: Clustering of Sparse Longitudinal Data with
#'          Automatic Selection of the Number of Clusters},
#'   author={Noukpoape, Komlan Midodzi},
#'   year={2026},
#'   note={R package version 1.0.0},
#'   url={https://github.com/noukpoape/InfoMergeCluster}
#' }
#' }
#'
#' @keywords internal
#' @importFrom stats predict rnorm runif rbinom
#' @importFrom graphics lines legend
#' @importFrom grDevices hcl.colors
"_PACKAGE"
