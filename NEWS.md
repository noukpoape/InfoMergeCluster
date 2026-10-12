# InfoMergeCluster 1.0.0

* First release.
* `infomergecluster()`: clustering of sparse longitudinal data with penalized
  B-spline mean curves, random effect (intercept, linear, spline or automatic
  choice), conjugate priors on the variance structure calibrated on the
  empirical variogram, and automatic (BIC/ICL) or imposed (`K`) number of
  clusters.
* Methods `print()`, `plot()` and `predict()` (confidence and prediction
  envelopes), `fitted_curves()` and `infomergecluster_stability()`.
* Data generator `f.data.generation()` (scenarios 1 to 4, each with reference
  group sizes and shifts used by default) and `plot.data()`.
* `infomergecluster_citation()`: BibTeX entry (Google Scholar format) or plain-text
  reference of the package; the entry is also installed as `CITATION.bib`.
