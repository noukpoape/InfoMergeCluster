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

```r
library(InfoMergeCluster)

# simulated data: 300 subjects, 3 true clusters, 2 to 5 measurements each
dat <- simulate_sparse_longitudinal(N = 300, K = 3, seed = 1)

# long data.frame: 1st column = subject id, second-to-last = time, last = response
fit <- infomergecluster(dat[, c("id", "time", "y")], K_max = 6)
fit                                   # selected K, proportions, criteria
plot(fit, envelope = TRUE)            # trajectories, mean curves, 95% band

# number of clusters imposed by the user (no selection)
fit3 <- infomergecluster(dat[, c("id", "time", "y")], K = 3)

# mean curves with prediction envelope, predicted individual curves
head(predict(fit3, interval = "prediction"))
curves <- fitted_curves(fit3)

# agreement with the true classes
adjusted_rand(fit3$cluster, dat$true_class[!duplicated(dat$id)])

# stability across seeds
st <- infomergecluster_stability(dat[, c("id", "time", "y")], K_max = 6, seeds = 1:5)
st$stability$mean_ari
```

Data can also be supplied as a list with one data.frame or matrix per subject
(second-to-last column = time, last column = response):

```r
set.seed(1)
gen <- f.data.generation(senario = 1, n.k = c(20, 30))
plot.data(gen$liste)
fit <- infomergecluster(gen$liste, K = 2)
```

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
