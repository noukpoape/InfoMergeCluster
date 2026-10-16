# InfoMergeCluster 1.0.0

Package rebuilt from version 1 of the InfoMergeCluster code.

## Breaking changes

* The main function is renamed `infomergecluster()` (formerly `infomerge()`);
  fitted objects have class `"infomergecluster"`, and `infomerge_stability()`
  becomes `infomergecluster_stability()`.
* The priors on the variance structure are always active and calibrated
  automatically on the variogram: the arguments `prior` and `s2_min` are
  removed (no lower bound on `sigma2` is needed any more).
* `re = "auto"` is removed: the random-effect structure (`"linear"`,
  `"intercept"` or `"spline"`) is chosen by the user and a single fit is
  performed.
* `simulate_sparse_longitudinal()` and `plot_data_cluster()` are removed;
  use `f.data.generation()` to simulate test data.
* An imposed `K` may now go up to N / 2 (formerly N / 3).

## New features

* Covariates (`n_cov`): covariates with cluster-specific effects on the mean
  trajectories (`gamma`), read by column position, centred, with dummy
  coding of factors. `predict()` gains an argument `covariates` to evaluate
  the mean curves at given covariate values.
* Covariate-dependent cluster proportions (`weights = "covariates"`):
  multinomial logistic model on the subject's mean covariates (`alpha`).
* `f.data.generation()` / `fsenario()`: scenario parameters reorganised
  (same random draws for a given seed); `plot_data(legend = TRUE)` now draws
  a legend of the groups.
* Messages and printed output are in English.

# InfoMergeCluster 0.9.1

* `infomerge()` gains an argument `K`. When `K` is given, the model is fitted
  with exactly `K` clusters and the number of clusters is no longer selected
  (`K_max` is then ignored, clusters are not pruned by `min_size`, and a
  warning is issued if a cluster is smaller than `min_size`). `K = NULL`
  (default) keeps the automatic selection.
* `print()` and `plot()` indicate when K was fixed by the user.
* Citation information: `citation("InfoMergeCluster")` (file `inst/CITATION`),
  `CITATION.cff` for GitHub's "Cite this repository" button, and a
  Citation section in the README.

# InfoMergeCluster 0.9.0

First release as a package (code version v9).

* Code split into thematic files in `R/`, roxygen documentation and `man/`
  files, `NAMESPACE` with registered S3 methods (`print`, `predict`, `plot`
  for class `infomerge`).
* `plot.data()` and `plot.data.cluster()` are renamed `plot_data()` and
  `plot_data_cluster()`: with a dot in their names, R treated them as
  `plot()` methods for classes "data" and "data.cluster".
* `simulate_sparse_longitudinal()` now returns its columns in the order
  `id, true_class, time, y`. Previously `true_class` came last, and
  `infomerge()`, which reads the response from the last column, would have
  clustered subjects on the true class instead of `y`.
* Accented characters in R string literals are written as `\uXXXX` escapes
  (a portability requirement of `R CMD check`); printed messages are
  unchanged.
* `testthat` tests added.

# Code history (before the package)

* v9: conjugate priors `G ~ IW(nu0, Psi0)`, `s2 ~ IG(a, b)`, with
  hyperparameters calibrated on the empirical variogram of the residuals.
* v8: penalised B-splines (O'Sullivan) with an exact integrated penalty.
* v7: `data` accepts a long data frame or a list of data frames/matrices.
* v6: minimum cluster size (`min_size`).
* v5: confidence and prediction envelopes.
* v4: `infomerge_stability()`, warning when the residual variance is close
  to zero.
* v3: P-splines with estimated smoothing parameter (Fellner-Schall).
* v2: adjustable random-effect structure (`re`), lower bound `s2_min`.
