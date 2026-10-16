# InfoMergeCluster

[![R-CMD-check](https://github.com/noukpoape/InfoMergeCluster/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noukpoape/InfoMergeCluster/actions/workflows/R-CMD-check.yaml)

Clustering of sparse longitudinal data, with automatic choice of the number
of clusters.

Model: `y_i | z_i = k ~ N(B_i beta_k + X_i gamma_k, Z_i G Z_i' + s2 I)`.
Cluster mean curves are penalised B-splines (O'Sullivan); optional
covariates have cluster-specific effects `gamma_k`, and the cluster
proportions may depend on the covariates; the random effect is a random
intercept, a random line or a spline. Conjugate priors on `(G, s2)` are
calibrated automatically on the empirical variogram of the residuals. A
greedy merge tree based on the exact ICL initialises a SQUAREM-accelerated
EM algorithm; K is selected by BIC (or ICL), or fixed by the user.

Author: Komlan Midodzi Noukpoape (<noukpoape@gmail.com>)

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("noukpoape/InfoMergeCluster")
```

The package only depends on base R packages (stats, splines, graphics,
grDevices, utils); no compilation is needed.

## Example

```r
library(InfoMergeCluster)

# simulated curves in two groups (list and data-frame formats)
set.seed(1)
g <- f.data.generation(senario = 1, n.k = c(30, 40), n.cts = c(0, 2))
plot_data(g$liste, legend = TRUE)

# number of clusters selected automatically
fit <- infomergecluster(g$liste, K_max = 6)
fit
plot(fit, envelope = TRUE, envelope_type = "prediction")
head(predict(fit, interval = "confidence"))
adjusted_rand(fit$cluster, attr(g$liste, "groupe"))

# number of clusters fixed by the user
fit2 <- infomergecluster(g$df, K = 2)

# with a covariate (columns: id, time, response, covariate)
d <- g$df[, c("id", "x", "y")]
age <- rnorm(70, 50, 10)
d$age <- age[d$id]
fit3 <- infomergecluster(d, K = 2, n_cov = 1)
fit3$gamma                                   # cluster-specific effects
predict(fit3, times = 0.5, covariates = 60)  # mean curves at age 60

# cluster proportions depending on the covariate
fit4 <- infomergecluster(d, K = 2, n_cov = 1, weights = "covariates")
fit4$alpha

# stability across seeds
st <- infomergecluster_stability(g$liste, K_max = 5, seeds = 1:5)
st$stability$K_freq
```

## Data format

Columns are identified by position, without names.

Without covariates (`n_cov = 0`, default):

- **long data frame**: first column = subject identifier, second to last =
  time, last = response (columns in between are ignored);
- **list** of data frames or matrices, one element per subject: second to
  last column = time, last = response.

With `n_cov` covariates:

- **long data frame**: column 1 = identifier, 2 = time, 3 = response,
  columns 4 to `3 + n_cov` = covariates;
- **list**: column 1 = time, 2 = response, columns 3 to `2 + n_cov` =
  covariates.

Covariates take one value per measurement (constant within a subject or
time-varying); factors and character columns are coded by dummy variables,
and all covariates are centred.

## Exported functions

| Function | Purpose |
|---|---|
| `infomergecluster()` | fit the model; K selected automatically or fixed with `K` |
| `print()`, `plot()`, `predict()` | methods for `infomergecluster` objects |
| `fitted_curves()` | predicted individual trajectories (BLUP) |
| `infomergecluster_stability()` | rerun over several seeds, agreement between partitions |
| `f.data.generation()`, `fsenario()` | simulated test data (scenarios 1 to 4) |
| `plot_data()` | plot of simulated curves with their true groups |
| `adjusted_rand()` | adjusted Rand index |

## Citation

If you use InfoMergeCluster in your work, please cite it. In R:

```r
citation("InfoMergeCluster")
```

which gives:

> Noukpoape, K. M. (2026). *InfoMergeCluster: Clustering of Sparse
> Longitudinal Data with Automatic Choice of the Number of Clusters*.
> R package version 1.0.0. https://github.com/noukpoape/InfoMergeCluster

BibTeX entry:

```bibtex
@Manual{InfoMergeCluster,
  title  = {{InfoMergeCluster}: Clustering of Sparse Longitudinal Data with
            Automatic Choice of the Number of Clusters},
  author = {Komlan Midodzi Noukpoape},
  year   = {2026},
  note   = {R package version 1.0.0},
  url    = {https://github.com/noukpoape/InfoMergeCluster},
}
```

GitHub also offers these references through the **Cite this repository**
button (file `CITATION.cff`).

## Documentation

The documentation is written as roxygen comments in `R/`. After editing it,
regenerate `man/` and `NAMESPACE` with:

```r
roxygen2::roxygenise("InfoMergeCluster")
```
