# RNOmni: Rank Normal Omnibus Association Test

Zachary R. McCaw <br>
Updated: 2026-02-24

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RNOmni)](https://cran.r-project.org/package=RNOmni)
[![](https://cranlogs.r-pkg.org/badges/grand-total/RNOmni)](https://CRAN.R-project.org/package=RNOmni)

Genetic association tests based on the **inverse normal transformation (INT)**. Recommended for continuous traits with non-normally distributed residuals. INT-based tests control type I error when standard linear regression does not (e.g. skewed or kurtotic residuals) and typically **outperform** linear regression in power.

- **D-INT** (direct): the phenotype is rank-normalized, then tested.
- **I-INT** (indirect): phenotypic residuals are rank-normalized, then tested.
- **O-INT** (omnibus): combines D-INT and I-INT via Cauchy combination for a single robust test.

Reference: [McCaw et al., *Biometrics* (2020)](https://doi.org/10.1111/biom.13214).

## Installation


``` r
# CRAN.
install.packages("RNOmni")

# GitHub.
remotes::install_github("zrmacc/RNOmni", build_vignettes = TRUE)
```

## Main functions

| Function | Description |
|----------|-------------|
| `OINT()` | Omnibus INT test |
| `DINT()` | Direct INT test |
| `IINT()` | Indirect INT test |
| `BAT()` | Basic association test (no transformation) |
| `RankNorm()` | Rank-based inverse normal transform |
| `OmniP()` | Cauchy combination of p-values |

## Quick example


``` r
library(RNOmni)
set.seed(100)
n <- 500
X <- cbind(1, rnorm(n))
G <- replicate(100, rbinom(n, 2, 0.25))
storage.mode(G) <- "numeric"
y <- exp(X %*% c(1, 0.5) + rnorm(n))   # skewed phenotype
p <- OINT(y = y, G = G, X = X, simple = TRUE)
```

## Vignette

[Package vignette](https://htmlpreview.github.io/?https://github.com/zrmacc/RNOmni/blob/master/vignettes/RNOmni.html).
