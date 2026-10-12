# InfoMergeCluster

[![R-CMD-check](https://github.com/noukpoape/InfoMergeCluster/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noukpoape/InfoMergeCluster/actions/workflows/R-CMD-check.yaml)

**Clustering of sparse longitudinal data with automatic selection of the number of clusters.**

`InfoMergeCluster` groups subjects followed over time — with few, irregular
measurements per subject — into clusters of similar trajectories.

- Each cluster has a smooth mean trajectory (penalized B-splines, smoothing
  parameter estimated automatically).
- Individual deviations are captured by a random effect: intercept, linear or
  spline (`re = "auto"` picks the best).
- Conjugate priors on the variance structure, calibrated on the empirical
  variogram of the residuals, stabilize the fit.
- The number of clusters is selected by BIC or ICL, or imposed with `K`.
- Confidence and prediction envelopes, predicted individual trajectories and a
  stability diagnostic across random seeds are included.

The package only depends on base R (`stats`, `splines`, `graphics`,
`grDevices`, `utils`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("noukpoape/InfoMergeCluster")
```

## Usage

The example below uses scenario 2 of the built-in data generator:
50 curves in 2 groups (20 and 30 curves), each observed at a few random times.

```r
library(InfoMergeCluster)

# scenario 2: 2 groups of 20 and 30 curves
set.seed(2)
gen <- f.data.generation(senario = 2, n.k = c(20, 30))
true_groups <- attr(gen$liste, "groupe")

# trajectories with their true groups
plot.data(gen$liste, legend = TRUE)

# automatic selection of the number of clusters (BIC)
# gen$df is a long data.frame: 1st column = id, second-to-last = time, last = response
fit <- infomergecluster(gen$df, K_max = 5)
fit                                   # selected K, proportions, criteria
plot(fit, envelope = TRUE)            # trajectories, mean curves, 95% band

# number of clusters imposed by the user (no selection)
fit2 <- infomergecluster(gen$df, K = 2)

# the list format (one matrix per curve) gives the same kind of fit
fit2_list <- infomergecluster(gen$liste, K = 2)

# agreement with the true groups
table(fit2$cluster, true_groups)
adjusted_rand(fit2$cluster, true_groups)

# mean curves with prediction envelope, predicted individual curves
head(predict(fit2, interval = "prediction"))
curves <- fitted_curves(fit2)

# stability across seeds
st <- infomergecluster_stability(gen$df, K_max = 5, seeds = 1:5)
st$stability$mean_ari
```

Larger simulated data sets (several hundred subjects, up to 6 clusters) can be
generated with `simulate_sparse_longitudinal()`.

## Main arguments of `infomergecluster()`

| Argument    | Meaning                                                              |
|-------------|----------------------------------------------------------------------|
| `K`         | imposed number of clusters (`NULL` = automatic selection)            |
| `K_max`     | largest number of clusters examined when `K = NULL`                  |
| `criterion` | `"BIC"` or `"ICL"`                                                   |
| `re`        | random effect: `"linear"`, `"intercept"`, `"spline"` or `"auto"`     |
| `df`        | number of B-spline functions                                         |
| `prior`     | priors on the variance structure (`NULL` = maximum likelihood)       |
| `min_size`  | minimum cluster size                                                 |

See `?infomergecluster` for the full list.

## Author

Komlan Midodzi NOUKPOAPE — <noukpoape@gmail.com>

## License

MIT © Komlan Midodzi NOUKPOAPE
