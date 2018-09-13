---
title: "README"
author: "Zachary McCaw"
date: "2018-09-13"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Motivating Example](#motivating-example)
* [Data](#data)
* [Omnibus Test](#omnibus-test)
* [Basic Association Test](#basic-association-test)
* [Direct Inverse Normal Transformation](#direct-inverse-normal-transformation)
* [Indirect Inverse Normal Transformation](#indirect-inverse-normal-transformation)
* [Notes](#notes)

# Motivating Example

Consider genetic association testing with a continuous trait whose residual distribution is far from normal. Departures from normality may take the form of skew towards a particular direction or an excess of extreme (outlying) residuals. A commonly right-skewed trait is body mass index. Other traits with potentially non-normal residuals include lung function measurements and polysomnography signals. When the departure from normality is severe, direct application of standard association tests, even those that only depend on asymptotic normality, can lead to an excess of false positive associations, and thereby invalidate inference. 

The rank based inverse normal transformation (INT) has been proposed to counteract departures from normality. During INT, the sample measurements are first mapped to the probability scale by replacing the observed values with fractional ranks. Next, these probabilities are then transformed into Z-scores using the probit function. To demonstrate, in the following a sample of size $n=1000$ is drawn from the $\chi_{1}^{2}$ distribution. After transformation, the distribution of the measurements in the sample is indistinguishable from normal. 

INT of the outcome in an association model is motivated by the fact that test statistics from distributions that are nearly normal converge in distribution faster than test statistics from distributions that are far from normal. By imposing normality on the outcome in the association model, the residual distribution, while not necessarily normal, cannot deviate too excessively. 


```r
library(RNOmni);
# Sample from the chi-1 square distribution
y = rchisq(n=1000,df=1);
# Rank-normalize
z = rankNorm(y);
```
<img src="README_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

# Data
## Simulated Data
In the following, data are simulated for $n=10^{3}$ subjects. Genotypes are drawn for $10^{3}$ loci in linkage equilibrium with minor allele frequency $0.25$. The model matrix `X` contains an intercept, four standard normal covariates `Z`, and the first four principal components of the genetic relatedness matrix. The intercept is set to one, and the remaining regression coefficients are simulated as random effects. The proportion of residual variation explained by covariates is 20%, while the proportion of residual variation explained by principle components is 5%. Two phenotypes with additive residuals are simulated. The first `yn` has standard normal residuals, while the second `yt` has $t_{3}$ residuals, scaled to have unit variance. 


```r
set.seed(100);
# Sample size
n = 1e3;
## Simulate genotypes
G = replicate(rbinom(n,size=2,prob=0.25),n=1e3);
storage.mode(G) = "numeric";
# Genetic principal components
S = svd(scale(G))$u[,1:4];
S = scale(S);
# Covariates
Z = scale(matrix(rnorm(n*4),nrow=n));
# Overall design
X = cbind(1,Z,S);
# Coefficient
b = c(1,rnorm(n=4,sd=1/sqrt(15)),rnorm(n=4,sd=1/sqrt(60)));
# Linear predictor
h = as.numeric(X%*%b);
# Normal phenotype
yn = h+rnorm(n);
# T-3 phenotype
yt = h+rt(n,df=3)/sqrt(3);
```
<img src="README_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

## Data Formatting
The outcome `y` is expected as a numeric vecotr. Genotypes `G` are expected as a numeric matrix, with subjects are rows. If adjusting for covariates or population structure, `X` is expected as a numeric matrix, with an intercept included. Factors and interactions are expanded in advance, e.g. using `model.matrix`. Missingness is not expected in either the outcome vector `y` or the model matrix `X`. Observations missing genotype information `G` are excluded from association testing only at those loci where the genotype is missing. 

# Omnibus Test
`RNOmni` performs an adaptive test of association between the loci in $G$ and the phenotype $y$, while adjusting for the model matrix $X$. Internally, `RNOmni` conducts two tests of association. In direct INT [DINT](#direct-inverse-normal-transformation), the INT is applied directly to the phenotype, and association testing is performed using the transformed phenotype. In indirect INT [IINT](#indirect-inverse-normal-transformation), the phenotype is regressed on the model matrix, and association testing is conducted on the INT-transformed phenotypic residuals. The omnibus statistic is the maximum of the DINT and IINT statistics. Hence, the omnibus statistics select which of DINT and IINT provide more evidence against the null. Synthesizing these complementary approaches affords the omnibus test robustness to the distribution of phenotypic residuals. In simulations against non-normal phenotypes, the omnibus test controlled the type I error in the absence of genetic associations, and improved power in the presence of genetic associations. Under the same settings, standard linear regression variously failed to control the type I error in the absence of associations, and was underpowered in the presence of associations. 

Assigning a $p$-value to the omnibus statistic requires an estimate of the correlation $\rho$ between the test statistics provided by `DINT` and `IINT` under the null. When many loci are under consideration, a computationally efficient strategy is to estimate $\rho$ by taking the empirical correlation between the test statistics across loci. Alternatively, the bootstrap is available when there are fewer loci under consideration, or when locus-specific correlation estimates are desired. In simulations, the value of $\rho$ estimated by taking the correlation across loci agreed with the average of the locus-specific estimates obtained using bootstrap. Finally, the user may estimate the correlation externally, then manually specify the value of $\rho$. 

By default, the output of `RNOmni` is a numeric matrix of $p$-values, with rows corresponding to the loci (columns) of $G$. The columns are the $p$-values from the DINT, the IINT, and the omnibus tests, respectively. If `keep.stats=T`, the test statistics are retained. If `keep.rho=T`, the estimated correlation between the DINT and IINT test statistics is retained. 


```r
cat("Omnibus Test, Normal Phenotype, Average Correaltion Method\n");
pn = RNOmni(y=yn,G=G,X=X,method="AvgCorr");
round(head(pn),digits=3);
cat("\n");
cat("Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method\n");
pn = RNOmni(y=yn,G=G[,1:10],X=X,method="Bootstrap",B=100);
round(head(pn),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Average Correaltion Method\n");
pt = RNOmni(y=yt,G=G,X=X,method="AvgCorr");
round(head(pt),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method\n");
pt = RNOmni(y=yt,G=G[,1:10],X=X,method="Bootstrap",keep.rho=T,B=100);
round(head(pt),digits=3);
cat("\n");
cat("Replicate the Omnibus Test on the T3 Phenotype, Manually Specifying Correlation\n");
pt = RNOmni(y=yt,G=G,X=X,method="Manual",set.rho=pt[,"Corr"],keep.rho=T);
round(head(pt),digits=3);
cat("\n");
```

```
## Omnibus Test, Normal Phenotype, Average Correaltion Method
##    DINT  IINT  Omni
## 1 0.471 0.481 0.494
## 2 0.454 0.431 0.453
## 3 0.720 0.718 0.737
## 4 0.342 0.356 0.363
## 5 0.056 0.060 0.063
## 6 0.888 0.844 0.858
## 
## Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method
##    DINT  IINT  Omni
## 1 0.471 0.481 0.485
## 2 0.454 0.431 0.445
## 3 0.720 0.718 0.737
## 4 0.342 0.356 0.354
## 5 0.056 0.060 0.060
## 6 0.888 0.844 0.856
## 
## Omnibus Test, T3 Phenotype, Average Correaltion Method
##    DINT  IINT  Omni
## 1 0.708 0.624 0.673
## 2 0.290 0.341 0.335
## 3 0.400 0.418 0.450
## 4 0.579 0.748 0.629
## 5 0.148 0.122 0.148
## 6 0.398 0.501 0.448
## 
## Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method
##    DINT  IINT  Omni  Corr
## 1 0.708 0.624 0.642 0.993
## 2 0.290 0.341 0.304 0.995
## 3 0.400 0.418 0.416 0.995
## 4 0.579 0.748 0.591 0.997
## 5 0.148 0.122 0.131 0.994
## 6 0.398 0.501 0.426 0.984
## 
## Replicate the Omnibus Test on the T3 Phenotype, Manually Specifying Correlation
##    DINT  IINT  Omni  Corr
## 1 0.708 0.624 0.642 0.993
## 2 0.290 0.341 0.304 0.995
## 3 0.400 0.418 0.416 0.995
## 4 0.579 0.748 0.591 0.997
## 5 0.148 0.122 0.131 0.994
## 6 0.398 0.501 0.426 0.984
```

# Basic Association Test
`BAT` regresses the untransformed phenotype `y` on genotype at each locus in `G`, adjusting for the model matrix `X`. A $p$-value assessing the null hypothesis of no genetic association is calculated using a score test. The output is a numeric matrix, including the score statistic and $p$-value for each locus in `G`.

```r
# Basic Association Test, Normal Phenotype
pn = BAT(y=yn,G=G,X=X);
round(head(pn),digits=3);
# Basic Association Test, T3 Phenotype
pt = BAT(y=yt,G=G,X=X);
round(head(pt),digits=3);
```

```
##   Score     P
## 1 0.514 0.474
## 2 0.496 0.481
## 3 0.143 0.705
## 4 0.840 0.360
## 5 3.646 0.056
## 6 0.045 0.832
##   Score     P
## 1 0.095 0.757
## 2 0.498 0.481
## 3 0.427 0.513
## 4 0.104 0.747
## 5 2.098 0.148
## 6 1.852 0.174
```

# Direct Inverse Normal Transformation
`DINT` regresses the INT-transformed phenotype `y` on genotype at each locus in `G`, adjusting for the model matrix `X`. A $p$-value assessing the null hypothesis of no genetic association is calculated using a score test. The output is a numeric matrix, including the score statistic and $p$-value for each locus in `G`.


```r
# Direct INT Test, Normal Phenotype
pn = DINT(y=yn,G=G,X=X);
round(head(pn),digits=3);
# Direct INT Test, T3 Phenotype
pt = DINT(y=yt,G=G,X=X);
round(head(pt),digits=3);
```

```
##   Score     P
## 1 0.519 0.471
## 2 0.561 0.454
## 3 0.129 0.720
## 4 0.904 0.342
## 5 3.650 0.056
## 6 0.020 0.888
##   Score     P
## 1 0.140 0.708
## 2 1.119 0.290
## 3 0.708 0.400
## 4 0.309 0.579
## 5 2.091 0.148
## 6 0.716 0.398
```

# Indirect Inverse Normal Transformation
`IINT` implements a two-stage association test. In the first stage, the untransformed phenotype `y` is regressed on the model matrix `X` to obtain phenotypic residuals. Likewise, the genotype matrix `G` is regressed on the model matrix `X` to obtain genotypic residuals. In the second stage, the INT-transformed phenotypic residuals are regressed on the genotypic residuals. The output is a numeric matrix, including the Wald statistic and $p$-value for each locus in `G`.


```r
# Indirect INT Test, Normal Phenotype
pn = IINT(y=yn,G=G,X=X);
round(head(pn),digits=3);
# Indirect INT Test, T3 Phenotype
pt = IINT(y=yt,G=G,X=X);
round(head(pt),digits=3);
```

```
##    Wald     P
## 1 0.497 0.481
## 2 0.620 0.431
## 3 0.130 0.718
## 4 0.852 0.356
## 5 3.547 0.060
## 6 0.039 0.844
##    Wald     P
## 1 0.240 0.624
## 2 0.908 0.341
## 3 0.656 0.418
## 4 0.104 0.748
## 5 2.387 0.122
## 6 0.454 0.501
```

# Notes

## Definition of the Rank Based Inverse Normal Transformation
Suppose that a continuous measurement $u_{i}$ is observed for each of $n$ subjects. Let $\text{rank}(u_{i})$ denote the sample rank of $u_{i}$ when the measurements are placed in ascending order. The rank based inverse normal transformation is defined as:

$$
\text{INT}(u_{i}) = \Phi^{-1}\left[\frac{\text{rank}(u_{i})-k}{n-2k+1}\right] 
$$

Here $\Phi^{-1}$ is the probit function, and $k\in(0,1/2)$ is an adjustable offset. By default, the Blom offset of $k=3/8$ is adopted.

## Parallelization
All association tests have the option of being run in parallel. To do so, register a parallel backend, e.g. `doMC::registerDoMC(cores=4)`, then specify the `parallel=T` option.
