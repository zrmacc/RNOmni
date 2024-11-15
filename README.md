# Rank Normal Omnibus Association Test
Zachary R. McCaw <br>
Updated: 2024-11-09

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RNOmni)](https://cran.r-project.org/package=RNOmni)
[![](https://cranlogs.r-pkg.org/badges/grand-total/RNOmni)](https://CRAN.R-project.org/package=RNOmni)

This package implements genetic association tests based on the inverse normal transform (INT). These tests are recommend for continuous traits with non-normally distributed residuals. INT-based tests robustly control the type I error in settings where standard linear regression does not, as when the residual distribution exhibits excess skew or kurtosis. Moreover, INT-based tests dominate standard linear regression in terms of power. These tests may be classified into two types. In direct INT (D-INT), the phenotype is itself transformed. In indirect INT (I-INT), phenotypic residuals are transformed. The omnibus test (O-INT) adaptively combines D-INT and I-INT into a single robust and statistically powerful approach. See [Operating characteristics of the rank-based inverse normal transformation for quantitative trait analysis in genome-wide association studies](https://doi.org/10.1111/biom.13214).

# Installation

```r
devtools::install_github(repo = 'zrmacc/RNOmni')
```

# Package Vignette

Link to package [vignette](https://htmlpreview.github.io/?https://github.com/zrmacc/RNOmni/blob/master/vignettes/RNOmni.html).
