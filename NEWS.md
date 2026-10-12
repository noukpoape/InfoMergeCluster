# InfoMergeCluster 1.0.0

* First release.
* `infomergecluster()`: clustering of sparse longitudinal data with penalized
  B-spline mean curves, random effect (intercept, linear, spline or automatic
  choice), conjugate priors on the variance structure calibrated on the
  empirical variogram, and automatic (BIC/ICL) or imposed (`K`) number of
  clusters.
* Methods `print()`, `plot()` and `predict()` (confidence and prediction
  envelopes), `fitted_curves()` and `infomergecluster_stability()`.
* Data generators `simulate_sparse_longitudinal()` and `f.data.generation()`.
