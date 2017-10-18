# Rank Normal Omnibus Association Test
Zachary McCaw  
Updated: 10/17/17  

<style>
body {
text-align: justify}
</style>


## Example Data
Simulated data is available for $10^{3}$ subjects. Covariates include `Age` and `Sex`. Structure adjustments include the first two principal components, `pc1` and `pc2`, of the centered and scaled subject by locus genotype matrix. Genotypes at $10^{2}$ loci on chromosome one are included. Minor allele frequency for each locus falls in the range $[0.110, 0.360]$. Two independent phenotypes are available. `YN` has normally distributed residuals, while the log of `YL` has normally distributed residuals. 

```r
library(RNOmni);
## Example data
# Covariates
X = RNOmni::X;
# Structure adjustments
S = RNOmni::S;
# Genotypes
G = RNOmni::Geno;
# Summarize minor allele frequency
cat("Summary of Minor Allele ")
```

```
## Summary of Minor Allele
```

```r
summary(apply(G,MARGIN=2,FUN=mean)/2);
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.110   0.185   0.215   0.219   0.250   0.360
```

```r
# Phenotypes
Y = RNOmni::Y;
```

## Genetic Association Testing

#### Basic Association Test
`BAT` regresses the untransformed genotype $g$, covariates $X$, and adjustments for substructure $S$. The Wald statistic is used to calculate a $p$-value assessing the null hypothesis of no genotypic effect. The output is a numeric vector, with one $p$-value per locus in $G$.

```r
# Basic Association Test
p1.bat = RNOmni::BAT(y=Y[,1],G=G,X=X,S=S);
microbenchmark::microbenchmark(RNOmni::BAT(y=Y[,1],G=G,X=X,S=S));
```

```
## Unit: milliseconds
##                                          expr      min       lq     mean
##  RNOmni::BAT(y = Y[, 1], G = G, X = X, S = S) 12.47437 12.61138 17.90871
##    median       uq      max neval
##  13.05883 13.29421 170.6918   100
```
#### Direct Inverse Normal Transformation
`DINT` directly applies the rank-based inverse normal transformation (INT) to the phenotype $y$. The transformed phenotype $\text{INT}(y)$ is regressed on genotype $g$, covariates $X$, and adjustments for population structure $S$. The Wald statistic is used to calculate a $p$-value assessing the null hypothesis of no genotypic effect. The output is a numeric vector, with one $p$-value per locus in $G$.

```r
# Direct INT Test
p1.dint = RNOmni::DINT(y=Y[,1],G=G,X=X,S=S);
microbenchmark::microbenchmark(RNOmni::DINT(y=Y[,1],G=G,X=X,S=S));
```

```
## Unit: milliseconds
##                                           expr      min       lq     mean
##  RNOmni::DINT(y = Y[, 1], G = G, X = X, S = S) 12.93472 13.12559 17.00106
##    median      uq      max neval
##  13.60456 13.8557 125.2654   100
```

#### Partially Indirect Inverse Normal Transformation
`PIINT` implements a two-stage, INT-based association test. In the first stage, the phenotype $y$ is regressed on covariates $X$ to obtain residuals $e$. In the second stage, the transformed residuals $\text{INT}(e)$ are regressed on genotype $g$ and adjustments for population structure $S$. The Wald statistic is used to calculate a $p$-value assessing the null hypothesis of no genotypic effect. The output is a numeric vector, with one $p$-value per locus in $G$. 

```r
# Partially Indirect INT Test
p1.piint = RNOmni::PIINT(y=Y[,1],G=G,X=X,S=S);
microbenchmark::microbenchmark(RNOmni::PIINT(y=Y[,1],G=G,X=X,S=S));
```

```
## Unit: milliseconds
##                                            expr      min       lq     mean
##  RNOmni::PIINT(y = Y[, 1], G = G, X = X, S = S) 11.14058 11.25516 13.99609
##    median       uq      max neval
##  11.48239 11.89021 123.7649   100
```

#### Omnibus Inverse Normal Transformation Method
`RNOmni` implements an adaptive, INT-based association test. In the omnibus test, `DINT` and `PIINT` are each applied, then an omnibus statistic is calculated based on whichever approach provides more evidence against the null. Estimation of a $p$-value for the omnibus statistic requires an estimate of the correlation $\rho$ between the test statistics provided by `DINT` and `PIINT`. When the sample size and number of loci are both relatively large, an efficient estimate of $\rho$ is obtained by averaging across loci. When the sample size or number of loci are smaller, bootstrap can be used to provide a locus-specific estimate of $\rho$. 

```r
# Omnibus Test using Average Correlation
p1.omni.avg = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr");
microbenchmark::microbenchmark(RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr"));
```

```
## Unit: milliseconds
##                                                                 expr
##  RNOmni::RNOmni(y = Y[, 1], G = G, X = X, S = S, method = "AvgCorr")
##       min      lq     mean   median       uq     max neval
##  40.00689 41.1458 48.88991 41.42179 41.91861 156.401   100
```

```r
# Omnibus Test using Bootstrap Correlation
p1.omni.boot = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100,cores=12);
microbenchmark::microbenchmark(RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100,cores=12));
```

```
## Unit: milliseconds
##                                                                                             expr
##  RNOmni::RNOmni(y = Y[, 1], G = G, X = X, S = S, method = "Bootstrap",      B = 100, cores = 12)
##       min       lq     mean   median       uq     max neval
##  418.2567 435.1381 549.4566 599.1073 608.9738 960.517   100
```

## Comparison of p-values
#### Normal Phenotype
The following table compares the estimated $p$-values for the normal phenotype obtained by each association test. Note that for the normal phenotype, all approaches provided valid inference in simulation.

```r
# Combined Results
P = cbind("BAT"=p1.bat,"DINT"=p1.dint,"PIINT"=p1.piint,
          "Omni.Avg"=p1.omni.avg[,"RNOmni"],"Omni.Boot"=p1.omni.boot[,"RNOmni"]);
show(round(head(P),digits=4));
```

```
##         BAT   DINT  PIINT Omni.Avg Omni.Boot
## [1,] 0.6536 0.6323 0.6745   0.7595    0.6558
## [2,] 0.4070 0.3799 0.3890   0.5084    0.4087
## [3,] 0.7697 0.7400 0.7315   0.8418    0.7552
## [4,] 0.5547 0.5683 0.5427   0.6770    0.5760
## [5,] 0.8766 0.8454 0.8520   0.9234    0.8733
## [6,] 0.7219 0.7391 0.7092   0.8242    0.7379
```

```r
# Empirical size
cat("Empirical size")
```

```
## Empirical size
```

```r
apply((P<=0.05),MARGIN=2,FUN=mean);
```

```
##       BAT      DINT     PIINT  Omni.Avg Omni.Boot 
##      0.02      0.02      0.02      0.02      0.02
```
#### Log-Normal Phenotype
The following table compares the estimated $p$-values for the log-normal phenotype obtained by each association test. Note that the basic association test did not consistently control the type I error in simulation. 

```r
# Estimate p-values
p2.bat = RNOmni::BAT(y=Y[,2],G=G,X=X,S=S);
p2.omni.avg = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="AvgCorr");
p2.omni.boot = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Bootstrap",B=100,cores=12);
# Combined Results
P2 = cbind("BAT"=p2.bat,"DINT"=p2.omni.avg[,"DINT"],"PIINT"=p2.omni.avg[,"PIINT"],
          "Omni.Avg"=p2.omni.avg[,"RNOmni"],"Omni.Boot"=p2.omni.boot[,"RNOmni"]);
show(round(head(P2),digits=4));
```

```
##         BAT   DINT  PIINT Omni.Avg Omni.Boot
## [1,] 0.8343 0.3242 0.8252   0.4798    0.4542
## [2,] 0.3922 0.5092 0.4808   0.6552    0.6051
## [3,] 0.6859 0.8284 0.5796   0.7504    0.7083
## [4,] 0.5644 0.6298 0.6206   0.7863    0.7354
## [5,] 0.8994 0.7919 0.8406   0.9129    0.8757
## [6,] 0.1769 0.2082 0.4896   0.3292    0.2765
```

```r
# Empirical size
cat("Empirical size")
```

```
## Empirical size
```

```r
apply((P2<=0.05),MARGIN=2,FUN=mean);
```

```
##       BAT      DINT     PIINT  Omni.Avg Omni.Boot 
##      0.01      0.03      0.01      0.01      0.01
```
