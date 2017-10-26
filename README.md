README
================
Zachary McCaw
10/18/2017

Purpose
=======

Genetic association tests for continuous phenotypes often assume normally distributed residuals. When the distribution of phenotypic residuals is skewed or heavy tailed, naive application of methods assuming normality can lead to an excess of false positive associations under the null, and loss of power under the alternative. Rank based inverse normal transformation (INT) has been applied in such cases to mitigate departures from normality. `RNOmni` provides an omnibus INT-based association test, synthesizing two approaches found to robustly control the type I error. In simulations, the omnibus test both provided valid inference, and achieved power comparable to the more powerful of the component methods.

Package Vignette
================

Contents
--------

-   [Example Data](#example-data)
-   [Genetic Association Testing](#genetic-association-testing)
-   [Run Time](#run-time)
-   [Missingness](#missingness)

Example Data
------------

Simulated data is available for 10<sup>3</sup> subjects. Covariates include `Age` and `Sex`. Structure adjustments include the first two principal components, `pc1` and `pc2`, of the centered and scaled subject by locus genotype matrix. Genotypes at 10<sup>2</sup> loci on chromosome one are included. Sample minor allele frequency for each locus falls in the range \[0.110, 0.360\]. Two independent phenotypes are available. `YN` has normally distributed residuals, while the log of `YL` has normally distributed residuals.

``` r
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
G = RNOmni::Geno;
G[1:6,1:6];
cat("\n");
cat("Sample Minor Allele Frequency\n");
summary(apply(G,MARGIN=2,FUN=mean)/2);
cat("\n");
cat("Phenotypes\n");
Y = RNOmni::Y;
round(head(Y),digits=2);
```

    ## Covariates
    ##        Age Sex
    ## [1,] 50.75   0
    ## [2,] 46.92   1
    ## [3,] 52.81   0
    ## [4,] 48.42   0
    ## [5,] 48.05   1
    ## [6,] 47.75   1
    ## 
    ## Structure Adjustments
    ##        pc1    pc2
    ## [1,] -1.90 -11.54
    ## [2,]  5.43   1.75
    ## [3,] 16.76   3.23
    ## [4,] 19.29 -13.15
    ## [5,]  2.92  -1.49
    ## [6,] 10.92   7.55
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
    ##   0.110   0.185   0.215   0.219   0.250   0.360 
    ## 
    ## Phenotypes
    ##        YN    YL
    ## [1,] 2.42 27.29
    ## [2,] 4.15 42.35
    ## [3,] 2.65 64.34
    ## [4,] 3.13  7.97
    ## [5,] 5.08 97.49
    ## [6,] 2.14 17.84

Genetic Association Testing
---------------------------

#### Basic Association Test

`BAT` regresses the untransformed phenotype *y* on genotype *g*, covariates *X*, and adjustments for substructure *S*. The Wald statistic is used to calculate a *p*-value assessing the null hypothesis of no genotypic effect. The output is a numeric vector, with one *p*-value per locus in *G*.

``` r
# Basic Association Test
p1.bat = RNOmni::BAT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.bat),digits=3);
```

    ## [1] 0.654 0.407 0.770 0.555 0.877 0.722

#### Direct Inverse Normal Transformation

`DINT` directly applies the rank-based inverse normal transformation (INT) to the phenotype *y*. The transformed phenotype INT(*y*) is regressed on genotype *g*, covariates *X*, and adjustments for population structure *S*. The Wald statistic is used to calculate a *p*-value assessing the null hypothesis of no genotypic effect. The output is a numeric vector, with one *p*-value per locus in *G*.

``` r
# Direct INT Test
p1.dint = RNOmni::DINT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.dint),digits=3);
```

    ## [1] 0.632 0.380 0.740 0.568 0.845 0.739

#### Partially Indirect Inverse Normal Transformation

`PIINT` implements a two-stage, INT-based association test. In the first stage, the phenotype *y* is regressed on covariates *X* to obtain residuals *e*. In the second stage, the transformed residuals INT(*e*) are regressed on genotype *g* and adjustments for population structure *S*. The Wald statistic is used to calculate a *p*-value assessing the null hypothesis of no genotypic effect. The output is a numeric vector, with one *p*-value per locus in *G*.

``` r
# Partially Indirect INT Test
p1.piint = RNOmni::PIINT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.piint),digits=3);
```

    ## [1] 0.674 0.389 0.731 0.543 0.852 0.709

#### Omnibus Inverse Normal Transformation Method

`RNOmni` implements an adaptive, INT-based association test. In the omnibus test, `DINT` and `PIINT` are each applied, then an omnibus statistic is calculated based on whichever approach provides more evidence against the null. Estimation of a *p*-value for the omnibus statistic requires an estimate of the correlation *ρ* between the test statistics provided by `DINT` and `PIINT`. When the sample size and number of loci are both relatively large, an efficient estimate of *ρ* is obtained by averaging across loci. When the sample size or number of loci are smaller, bootstrap can be used to provide a locus-specific estimate of *ρ*.

``` r
cat("Omnibus Test, Average Correaltion Method\n");
p1.omni.avg = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr");
round(head(p1.omni.avg),digits=3);
cat("\n");
cat("Omnibus Test, Bootstrap Correaltion Method\n");
p1.omni.boot = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100,cores=12);
round(head(p1.omni.boot),digits=3);
cat("\n");
```

    ## Omnibus Test, Average Correaltion Method
    ##       DINT PIINT RNOmni
    ## [1,] 0.632 0.674  0.759
    ## [2,] 0.380 0.389  0.508
    ## [3,] 0.740 0.731  0.842
    ## [4,] 0.568 0.543  0.677
    ## [5,] 0.845 0.852  0.923
    ## [6,] 0.739 0.709  0.824
    ## 
    ## Omnibus Test, Bootstrap Correaltion Method
    ##       DINT PIINT RNOmni
    ## [1,] 0.632 0.674  0.662
    ## [2,] 0.380 0.389  0.404
    ## [3,] 0.740 0.731  0.754
    ## [4,] 0.568 0.543  0.569
    ## [5,] 0.845 0.852  0.862
    ## [6,] 0.739 0.709  0.740

#### Comparison of p-values

The following table compares the estimated *p*-values for the normal phenotype obtained by each association test. Note that for the normal phenotype, all approaches provided valid inference in simulation.

``` r
cat("Normal Phenotype, Combined Results\n");
P = cbind("BAT"=p1.bat,"DINT"=p1.dint,"PIINT"=p1.piint,
          "Omni.Avg"=p1.omni.avg[,"RNOmni"],"Omni.Boot"=p1.omni.boot[,"RNOmni"]);
show(round(head(P),digits=4));
cat("\n");
cat("Empirical Size\n");
apply((P<=0.05),MARGIN=2,FUN=mean);
```

    ## Normal Phenotype, Combined Results
    ##         BAT   DINT  PIINT Omni.Avg Omni.Boot
    ## [1,] 0.6536 0.6323 0.6745   0.7595    0.6615
    ## [2,] 0.4070 0.3799 0.3890   0.5084    0.4036
    ## [3,] 0.7697 0.7400 0.7315   0.8418    0.7541
    ## [4,] 0.5547 0.5683 0.5427   0.6770    0.5694
    ## [5,] 0.8766 0.8454 0.8520   0.9234    0.8615
    ## [6,] 0.7219 0.7391 0.7092   0.8242    0.7399
    ## 
    ## Empirical Size
    ##       BAT      DINT     PIINT  Omni.Avg Omni.Boot 
    ##      0.02      0.02      0.02      0.02      0.02

The following table compares the estimated *p*-values for the log-normal phenotype obtained by each association test. Note that the basic association test did not consistently control the type I error in simulation.

``` r
# Estimate p-values
p2.bat = RNOmni::BAT(y=Y[,2],G=G,X=X,S=S);
p2.omni.avg = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="AvgCorr");
p2.omni.boot = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Bootstrap",B=100,cores=12);
cat("Log Normal Phenotype, Combined Results\n");
P2 = cbind("BAT"=p2.bat,"DINT"=p2.omni.avg[,"DINT"],"PIINT"=p2.omni.avg[,"PIINT"],
          "Omni.Avg"=p2.omni.avg[,"RNOmni"],"Omni.Boot"=p2.omni.boot[,"RNOmni"]);
show(round(head(P2),digits=4));
cat("\n");
cat("Empirical Size\n");
apply((P2<=0.05),MARGIN=2,FUN=mean);
```

    ## Log Normal Phenotype, Combined Results
    ##         BAT   DINT  PIINT Omni.Avg Omni.Boot
    ## [1,] 0.8343 0.3242 0.8252   0.4798    0.4451
    ## [2,] 0.3922 0.5092 0.4808   0.6552    0.6035
    ## [3,] 0.6859 0.8284 0.5796   0.7504    0.7135
    ## [4,] 0.5644 0.6298 0.6206   0.7863    0.7114
    ## [5,] 0.8994 0.7919 0.8406   0.9129    0.8751
    ## [6,] 0.1769 0.2082 0.4896   0.3292    0.2866
    ## 
    ## Empirical Size
    ##       BAT      DINT     PIINT  Omni.Avg Omni.Boot 
    ##      0.01      0.03      0.01      0.01      0.01

Run Time
--------

During package development, the `BAT`, `DINT`, and `PIINT` each took a median of ∼25 ms to perform 10<sup>2</sup> association tests for 10<sup>3</sup> subjects. `RNOmni` using average correlation, which interval performs both `DINT` and `PIINT`, required a median of ∼65 ms. Using bootstrap to calculate position specific correlations increased the run time of `RNOmni` by a factor of 10.

``` r
library(RNOmni);
microbenchmark::microbenchmark(BAT(y=Y[,1],G=G,X=X,S=S),DINT(y=Y[,1],G=G,X=X,S=S),PIINT(y=Y[,1],G=G,X=X,S=S),
                               RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr"),
                               RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100,cores=12));
```

    ## Unit: milliseconds
    ##                                                                                     expr
    ##                                                     BAT(y = Y[, 1], G = G, X = X, S = S)
    ##                                                    DINT(y = Y[, 1], G = G, X = X, S = S)
    ##                                                   PIINT(y = Y[, 1], G = G, X = X, S = S)
    ##                              RNOmni(y = Y[, 1], G = G, X = X, S = S, method = "AvgCorr")
    ##  RNOmni(y = Y[, 1], G = G, X = X, S = S, method = "Bootstrap",      B = 100, cores = 12)
    ##        min        lq      mean    median        uq      max neval
    ##   19.95516  21.55670  31.58464  23.28569  29.17856 196.7335   100
    ##   20.29768  21.33566  31.04566  23.03815  28.35499 185.2453   100
    ##   17.06783  18.75149  28.60262  20.11432  22.68837 181.8650   100
    ##   53.71182  55.84626  69.73768  57.92123  63.89067 231.3433   100
    ##  551.28151 567.08678 661.95787 719.77288 740.73287 755.2562   100

Missingness
-----------

Observations missing either the phenotype *y* or the the substructure adjustments *S* are excluded. Missing covariates *X* are imputed to the median of the observed values. An observation missing genotype information *G* is excluded from association testing only at those loci were genotype is unobserved.

``` r
# Introduce Missingness
y.m = Y[,1];
y.m[sample(length(y.m),size=10,replace=F)] = NA;
G.m = G;
G.m[sample(length(G.m),size=1000,replace=F)] = NA;
X.m = X;
X.m[sample(length(X.m),size=100,replace=F)] = NA;
S.m = S;
S.m[sample(length(S.m),size=10,replace=F)] = NA;
# Association Testing after Missingness
p1.bat.m = RNOmni::BAT(y=y.m,G=G.m,X=X.m,S=S.m);
p1.dint.m = RNOmni::DINT(y=y.m,G=G.m,X=X.m,S=S.m);
p1.piint.m = RNOmni::PIINT(y=y.m,G=G.m,X=X.m,S=S.m);
p1.omni.avg.m = RNOmni::RNOmni(y=y.m,G=G.m,X=X.m,S=S.m);
p1.omni.boot.m = RNOmni::RNOmni(y=y.m,G=G.m,X=X.m,S=S.m,method="Bootstrap",B=100);
```

Tabulation of *p*-values in the presence of missingness:
