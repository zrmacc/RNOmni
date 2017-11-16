# Rank Normal Omnibus Association Test
Zachary McCaw  
Updated: 10/27/18  



## Contents

* [Motivating Example](#motivating-example)
* [Rank Normal Omnibus Test](#omnibus-association-test)
* [Additional Association Tests](#additional-association-tests)
* [Comparison of Association Tests](#comparison-of-association-tests)
* [Additional Details](#additional-details)

## Motivating Example
Motivation for investigating the rank based inverse normal transformation (INT) came from the study of obstructive sleep apnea (OSA). The gold standard measurement for diagnosing OSA is the apnea hypopnea index (AHI), a continuous, positively skewed trait with a distribution that resembles $\chi^{2}$. The residuals obtained by regressing AHI on genotype and covariates are often non-normal, and application of standard association tests leads to an excess of false positive associations, i.e. inflated type I error. To counteract departures from normality, INT transformation of AHI prior to genetic association testing has been proposed. As demonstrated in the following figure, INT of a continuous measurement causes the distribution of that measurement in the sample to appear normal. 


```r
# Chi-1 data
y = rchisq(n=1000,df=1);
# Rank-normalize
z = RNOmni::rankNormal(y);
```
<img src="Figs/A01-1.png" style="display: block; margin: auto;" />

#### Simulated Data
Within `RNOmni`, simulated data is available for $10^{3}$ subjects. Covariates include `Age` and `Sex`. Structure adjustments include the first two principal components, `pc1` and `pc2`, of the centered and scaled subject by locus genotype matrix. Genotypes at $10^{3}$ loci on chromosome one are also included. All loci are common, with sample minor allele frequency falling in the range $[0.230, 0.403]$. Two independent phenotypes were generated under the null hypothesis of no genotypic effect. `YN` has normally distributed residuals, while $`YT3$ has heavy tailed residuals from a $t_{3}$ distribution. The residual distributions were scaled to have unit variance. 


```r
library(RNOmni);
## Example data
X = RNOmni::X;
cat("Covariates\n");
round(head(X),digits=2);
cat("\n");
cat("Structure Adjustments\n");
S = RNOmni::S;
round(head(S),digits=2);
cat("\n");
cat("Genotype Matrix\n");
G = RNOmni::G;
G[1:6,1:6];
cat("\n");
cat("Sample Minor Allele Frequency\n");
summary(apply(G,MARGIN=2,FUN=mean)/2);
cat("\n");
cat("Phenotypes\n");
Y = RNOmni::Y;
round(head(Y),digits=2);
```

```
## Covariates
##          a s
## [1,] 48.33 0
## [2,] 45.02 1
## [3,] 52.74 1
## [4,] 50.27 0
## [5,] 50.91 1
## [6,] 48.08 1
## 
## Structure Adjustments
##        pc1   pc2
## [1,] -1.90 11.54
## [2,]  5.43 -1.75
## [3,] 16.76 -3.23
## [4,] 19.29 13.15
## [5,]  2.92  1.49
## [6,] 10.92 -7.55
## 
## Genotype Matrix
##      s1 s2 s3 s4 s5 s6
## [1,]  0  0  0  1  0  0
## [2,]  1  0  2  0  1  0
## [3,]  0  0  0  1  1  0
## [4,]  0  0  2  2  1  1
## [5,]  2  1  0  2  1  1
## [6,]  0  0  0  0  0  0
## 
## Sample Minor Allele Frequency
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2295  0.3050  0.3225  0.3217  0.3390  0.4030 
## 
## Phenotypes
##        YN  YT3
## [1,] 2.68 2.32
## [2,] 3.42 4.43
## [3,] 3.80 4.95
## [4,] 4.44 3.47
## [5,] 3.79 4.39
## [6,] 3.78 3.54
```
#### Normal QQ Plots for Phenotypic Residuals
The following quantile-quantile (QQ) plots demonstrate the effect of INT on model residuals. In blue are the residuals for standard linear regression of the normal phenotype on covariates $X$ and structure adjustments $S$. In dark green are the residual for standard linear regression of the $t_{3}$ phenotype on $X$ and $S$. The heavier tails of the $t_{3}$ distribution are evidenced by the departures of the observed quantiles from their expectations away from the origin.

In light green are the residuals from [direct INT](#direct-inverse-normal-transformation) (DINT) and [partially indirect INT](#partially-indirect-inverse-normal-transformation) (PIINT) of the $t_{3}$ phenotype. The validity of the INT-based association tests derives from the closer adherence of the observed quantiles to their expectations in these models. The omnibus test synthesizes the association statistics from the DINT and PIINT models.  

<img src="Figs/A03-1.png" style="display: block; margin: auto;" />

## Rank Normal Omnibus Test
`RNOmni` implements an adaptive test of association between the loci in $G$ and the phenotype $y$, while adjusting for covariates $X$ and population structure $S$. Internally, `RNOmni` conducts two association tests, `DINT` and `PIINT`, described below, then calculates an omnibus statistic based on whichever approach provides more evidence against the null hypothesis. Synthesizing two complementary, INT-based approaches, affords the omnibus test robustness to the distribution of phenotypic residuals. In the absence of a genotypic effect, `RNOmni` routinely controls the type I error. In the presence of a genotypic effect, `RNOmni` provides power comparable to the better of `DINT` and `PIINT`. 

Estimation of a $p$-value for the omnibus statistic requires an estimate of the correlation $\rho$ between the test statistics provided by `DINT` and `PIINT`. When the sample size and number of loci are both relatively large, a computationally efficient estimate of $\rho$ is obtained by averaging across loci. If either the sample size or the number of loci is relatively small, bootstrap can provide a locus-specific estimates of $\rho$. To accelerate the bootstrap, registered a parallel packend, e.g. `doMC::registerDoMC(cores=2)`, then pass the `parallel=T` option to `RNOmni`. 

The output of `RNOmni` is a numeric matrix of $p$-values, with rows corresponding to the rows of $G$. The columns are the $p$-values from the `DINT`, the `PIINT`, and the omnibus tests, respectively. Note that, without additional adjustment for multiple testing, taking the minimum $p$-value across each row would not result in a valid test of association. 

```r
cat("Omnibus Test, Normal Phenotype, Average Correaltion Method\n");
p1.omni.avg = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr");
round(head(p1.omni.avg),digits=3);
cat("\n");
cat("Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method\n");
set.seed(100);
p1.omni.boot = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100);
round(head(p1.omni.boot),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Average Correaltion Method\n");
p2.omni.avg = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="AvgCorr");
round(head(p2.omni.avg),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method\n");
p2.omni.boot = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Bootstrap",B=100);
round(head(p2.omni.boot),digits=3);
cat("\n");
```

```
## Omnibus Test, Normal Phenotype, Average Correaltion Method
##       DINT PIINT RNOmni
## [1,] 0.626 0.653  0.698
## [2,] 0.727 0.764  0.789
## [3,] 0.165 0.167  0.212
## [4,] 0.866 0.910  0.907
## [5,] 0.469 0.509  0.544
## [6,] 0.584 0.568  0.642
## 
## Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method
##       DINT PIINT RNOmni
## [1,] 0.626 0.653  0.661
## [2,] 0.727 0.764  0.757
## [3,] 0.165 0.167  0.188
## [4,] 0.866 0.910  0.884
## [5,] 0.469 0.509  0.503
## [6,] 0.584 0.568  0.598
## 
## Omnibus Test, T3 Phenotype, Average Correaltion Method
##       DINT PIINT RNOmni
## [1,] 0.751 0.690  0.736
## [2,] 0.540 0.531  0.583
## [3,] 0.201 0.249  0.238
## [4,] 0.192 0.197  0.228
## [5,] 0.329 0.332  0.376
## [6,] 0.462 0.431  0.482
## 
## Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method
##       DINT PIINT RNOmni
## [1,] 0.751 0.690  0.718
## [2,] 0.540 0.531  0.568
## [3,] 0.201 0.249  0.233
## [4,] 0.192 0.197  0.219
## [5,] 0.329 0.332  0.359
## [6,] 0.462 0.431  0.472
```
Since the phenotype was simulated under the null hypothesis of no genotypic effect, the expected false positive rate at $\alpha$ level $0.05$ is $5\%$. For both the normal and heavy tailed $t_{3}$ phenotypes, the $95\%$ confidence interval for the type I error includes the expected value of $0.05$. As shown in the [comparison of association tests](#comparison-of-association-tests), naively applying the [basic association test](#basic-association-test) leads to an excess of false positive associations in the latter case. 


```
## Type I Error of Rank Normal Omnibus Test:
##   Phenotype    Method  Size     L     U
## 1    Normal   AvgCorr 0.038 0.026 0.050
## 2    Normal Bootstrap 0.052 0.038 0.066
## 3        T3   AvgCorr 0.057 0.042 0.072
## 4        T3 Bootstrap 0.062 0.047 0.077
```

## Additional Association Tests
In addition to the omnibus test, three genetic association tests are implemented as part of `RNOmni`. These are the basic association test `BAT`, the direct INT method `DINT`, and the partially indirect INT method `PIINT`.

#### Basic Association Test
`BAT` regresses the untransformed phenotype $y$ on genotype at each locus in $G$, adjusting for covariates $X$ and population structure $S$. A $p$-value assessing the null hypothesis of no genotypic effect is estimated using the Wald statistic. The output is a numeric vector, with one $p$-value per row of $G$.

```r
# Basic Association Test
p1.bat = RNOmni::BAT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.bat),digits=3);
```

```
## [1] 0.679 0.737 0.165 0.901 0.477 0.568
```

#### Direct Inverse Normal Transformation
`DINT` regresses the transformed phenotype $\text{INT}(y)$ on genotype at each locus in $G$, adjusting for covariates $X$ and population structure $S$. A $p$-value assessing the null hypothesis of no genotypic effect is estimated via the Wald statistic. The output is a numeric vector, with one $p$-value per row of $G$.


```r
# Direct INT Test
p1.dint = RNOmni::DINT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.dint),digits=3);
```

```
## [1] 0.626 0.727 0.165 0.866 0.469 0.584
```

#### Partially Indirect Inverse Normal Transformation
`PIINT` implements a two-stage association test. In the first stage, the untransformed phenotype $y$ is regressed on covariates $X$ to obtain residuals $e$. In the second stage, the transformed residuals $\text{INT}(e)$ are regressed on genotype at each locus in $G$, adjusting for population structure $S$. A $p$-value assessing the null hypothesis of no genotypic effect is estimated via the Wald statistic. The output is a numeric vector, with one $p$-value per row of $G$.


```r
# Partially Indirect INT Test
p1.piint = RNOmni::PIINT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.piint),digits=3);
```

```
## [1] 0.653 0.764 0.167 0.910 0.509 0.568
```

In naming `PIINT`, "indirect" refers to the fact that residuals are formed prior to INT, rather than directly transforming the phenotype. "Partially" refers to the fact that residuals are formed w.r.t. covariates $X$, but not structure adjustments $S$. Performance of a fully indirect association test, in which $y$ is regressed on both $X$ and $S$ during residual formation, was investigated. The fully indirect association test did not consistently provide valid inference. 

## Comparison of Association Tests
The following figure depicts the estimated type I error for association tests against the normal and $t_{3}$ phenotypes at $\alpha$ level $0.05$. Point estimates are obtained by averaging an indicator of rejection, and the error bars provide $95\%$ confidence intervals. Although all methods perform comparable against the normal phenotype, the false positive rate is inflated when the basic association approach is applied in the presence of heavy tailed $t_{3}$ residuals. This problem is exacerbated when considering increasingly small $\alpha$ levels.  


<img src="Figs/C05-1.png" style="display: block; margin: auto;" />

## Additional Details

#### Run time
During package development, the `BAT`, `DINT`, and `PIINT` each took a median of $25$ to $30$ ms to perform $10^2$ association tests for $10^3$ subjects. `RNOmni` using average correlation, which internally performs both `DINT` and `PIINT`, required a median of $65$ to $70$ ms. Using bootstrap to calculate position specific correlations increased the run time of `RNOmni` by a factor of $9$ to $10$ while running $12$ cores in parallel. Using `RNOmni` with the `parallel=T` option is advised if using the bootstrap approach. 


```r
# Subset to 100 loci
H = G[1:100,];
# Time performance
library(microbenchmark);
microbenchmark(BAT(y=Y[,1],G=H,X=X,S=S),DINT(y=Y[,1],G=H,X=X,S=S),PIINT(y=Y[,1],G=H,X=X,S=S),
               RNOmni(y=Y[,1],G=H,X=X,S=S,method="AvgCorr"));
microbenchmark(RNOmni(y=Y[,1],G=H,X=X,S=S,method="Bootstrap",B=100),times=10);
```

```
## Unit: milliseconds
##                                                         expr      min
##                         BAT(y = Y[, 1], G = H, X = X, S = S) 21.40669
##                        DINT(y = Y[, 1], G = H, X = X, S = S) 21.89634
##                       PIINT(y = Y[, 1], G = H, X = X, S = S) 18.70127
##  RNOmni(y = Y[, 1], G = H, X = X, S = S, method = "AvgCorr") 57.16463
##        lq     mean   median       uq      max neval
##  23.37066 26.49970 24.89899 25.82675 188.3917   100
##  23.80037 26.50177 24.86587 25.96194 188.9738   100
##  20.36414 23.11901 21.44712 22.80345 183.9723   100
##  61.14118 65.02569 63.18766 66.12928 225.9110   100
## Unit: seconds
##                                                                         expr
##  RNOmni(y = Y[, 1], G = H, X = X, S = S, method = "Bootstrap",      B = 100)
##       min      lq   mean   median       uq      max neval
##  4.114643 4.13913 4.1983 4.167833 4.282667 4.348065    10
```

#### Missingness
Observations missing either the phenotype $y$ or the the structure adjustments $S$ are excluded. 
Missing covariates $X$ are imputed to the median of the observed values. An observation missing genotype information $G$ is excluded from association testing only at those loci were genotype is unobserved. 


```r
# Introduce Missingness
y.m = Y[,1];
y.m[sample(length(y.m),size=10,replace=F)] = NA;
G.m = G;
G.m[sample(length(G.m),size=10000,replace=F)] = NA;
X.m = X;
X.m[sample(length(X.m),size=100,replace=F)] = NA;
S.m = S;
S.m[sample(length(S.m),size=10,replace=F)] = NA;
# Association Testing after Missingness
pm.bat = RNOmni::BAT(y=y.m,G=G.m,X=X.m,S=S.m);
pm.dint = RNOmni::DINT(y=y.m,G=G.m,X=X.m,S=S.m);
pm.piint = RNOmni::PIINT(y=y.m,G=G.m,X=X.m,S=S.m);
pm.omni.avg = RNOmni::RNOmni(y=y.m,G=G.m,X=X.m,S=S.m);
pm.omni.boot = RNOmni::RNOmni(y=y.m,G=G.m,X=X.m,S=S.m,method="Bootstrap",B=100);
```
Below, size estimates for the normal phenotype in the absence and presence of missingness are tabulated:

```
## Normal Phenotype in the Absence of Missingness, Combined Results
##         Method  Size     L     U
## 1          BAT 0.052 0.038 0.066
## 2         DINT 0.061 0.046 0.076
## 3        PIINT 0.056 0.041 0.071
## 4 Omni.AvgCorr 0.038 0.026 0.050
## 5    Omni.Boot 0.052 0.038 0.066
## 
## Normal Phenotype in the Presence of Missingness, Combined Results
##      Method  Size     L     U
## 1       BAT 0.036 0.024 0.048
## 2      DINT 0.045 0.032 0.058
## 3     PIINT 0.039 0.027 0.051
## 4  Omni.Avg 0.023 0.014 0.032
## 5 Omni.Boot 0.034 0.023 0.045
```