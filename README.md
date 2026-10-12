# InfoMergeCluster

[![R-CMD-check](https://github.com/noukpoape/InfoMergeCluster/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noukpoape/InfoMergeCluster/actions/workflows/R-CMD-check.yaml)

**Clustering of sparse longitudinal data with automatic selection of the number of clusters.**

`InfoMergeCluster` groups subjects followed over time — with few, irregular
measurements per subject — into clusters of similar trajectories.

- Each cluster has a smooth mean trajectory (penalized B-splines, smoothing
  parameter estimated automatically).
- Individual deviations are captured by a random effect: intercept, linear
  (default) or spline.
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

The example below uses scenario 1 of the built-in data generator with its
reference settings: 50 curves in 2 groups (20 and 30 curves), each observed
at a few random times.

```r
library(InfoMergeCluster)

# scenario 1: 2 groups of 20 and 30 curves
set.seed(1)
gen <- f.data.generation(senario = 1)
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

## Simulated scenarios

`f.data.generation()` provides four scenarios. When `n.k` and `n.cts` are not
supplied, each scenario uses its reference settings:

| Scenario | Groups | Group sizes `n.k`                  | Shifts `n.cts`                  |
|:--------:|:------:|------------------------------------|---------------------------------|
| 1        | 2      | `c(20, 30)` (50 curves)            | `c(1.8, 2.4)`                   |
| 2        | 3      | `c(20, 30, 25)` (75 curves)        | `c(1.8, 2.4, 3.0)`              |
| 3        | 4      | `c(20, 30, 25, 35)` (110 curves)   | `c(1.8, 2.4, 3.0, 3.6)`         |
| 4        | 5      | `c(20, 30, 25, 35, 40)` (150 curves) | `c(1.8, 2.4, 3.0, 3.6, 4.2)`  |

```r
gen3 <- f.data.generation(senario = 3)                  # reference settings
gen3b <- f.data.generation(senario = 3, n.k = c(50, 50, 50, 50))   # custom sizes
```

## Main arguments of `infomergecluster()`

| Argument    | Meaning                                                              |
|-------------|----------------------------------------------------------------------|
| `K`         | imposed number of clusters (`NULL` = automatic selection)            |
| `K_max`     | largest number of clusters examined when `K = NULL`                  |
| `criterion` | `"BIC"` or `"ICL"`                                                   |
| `re`        | random effect: `"linear"` (default), `"intercept"` or `"spline"`     |
| `df`        | number of B-spline functions                                         |
| `prior`     | priors on the variance structure (`NULL` = maximum likelihood)       |
| `min_size`  | minimum cluster size                                                 |

See `?infomergecluster` for the full list.

## Citation

If you use `InfoMergeCluster` in a publication, please cite it as:

> Noukpoape, K. M. (2026). *InfoMergeCluster: Clustering of Sparse
> Longitudinal Data with Automatic Selection of the Number of Clusters.*
> R package version 1.0.0. https://github.com/noukpoape/InfoMergeCluster

BibTeX:

```bibtex
@manual{noukpoape2026infomergecluster,
  title={InfoMergeCluster: Clustering of Sparse Longitudinal Data with Automatic Selection of the Number of Clusters},
  author={Noukpoape, Komlan Midodzi},
  year={2026},
  note={R package version 1.0.0},
  url={https://github.com/noukpoape/InfoMergeCluster}
}
```

The package itself gives this reference:

```r
library(InfoMergeCluster)
infomergecluster_citation()                     # BibTeX entry above
infomergecluster_citation("text")               # plain-text reference
infomergecluster_citation(file = "refs.bib")    # append it to a .bib file
citation("InfoMergeCluster")                    # standard R citation

# the .bib file installed with the package
system.file("CITATION.bib", package = "InfoMergeCluster")
```

It is also available on GitHub through the **Cite this repository** button
(file `CITATION.cff`).

## Author

Komlan Midodzi Noukpoape — <noukpoape@gmail.com>

## License

MIT © Komlan Midodzi Noukpoape
