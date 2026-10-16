# InfoMergeCluster 1.0.0

* First release.
* `infomergecluster()`: clustering of sparse longitudinal data with a mixture
  of linear mixed models (penalized B-spline mean curves, random effect
  `"linear"`, `"intercept"` or `"spline"`), automatic (BIC/ICL) or imposed
  (`K`) number of clusters.
* Conjugate priors on the variance structure, always active and calibrated on
  the empirical variogram of the residuals (no floor on the residual variance).
* Covariates (`n_cov`) with cluster-specific effects on the trajectories, and
  cluster proportions depending on the covariates (`weights = "covariates"`);
  example with a covariate acting as a level shift (README and
  `?infomergecluster`, uses the suggested package `fda`).
* Methods `print()`, `plot()` and `predict()` (confidence and prediction
  envelopes, curves for given covariate values), `fitted_curves()` and
  `infomergecluster_stability()`.
* Data generator `f.data.generation()` (scenarios 1 to 4, reference settings
  documented in its help page) and `plot.data()`.
* `infomergecluster_citation()`: BibTeX entry (Google Scholar format) or
  plain-text reference of the package; the entry is also installed as
  `CITATION.bib`.
