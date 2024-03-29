% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IINT.R
\name{IINT}
\alias{IINT}
\title{Indirect-INT}
\usage{
IINT(y, G, X = NULL, k = 0.375, ties.method = "average", simple = FALSE)
}
\arguments{
\item{y}{Numeric phenotype vector.}

\item{G}{Genotype matrix with observations as rows, SNPs as columns.}

\item{X}{Model matrix of covariates and structure adjustments. Should include
an intercept. Omit to perform marginal tests of association.}

\item{k}{Offset applied during rank-normalization. See
\code{\link{RankNorm}}.}

\item{ties.method}{Method of breaking ties, passed to \code{base::rank}.}

\item{simple}{Return the p-values only?}
}
\value{
If \code{simple = TRUE}, returns a vector of p-values, one for each column
  of \code{G}. If \code{simple = FALSE}, returns a numeric matrix, including the
  Wald or Score statistic, its standard error, the Z-score, and the p-value.
}
\description{
Two-stage association testing procedure. In the first stage, phenotype 
\code{y} and genotype \code{G} are each regressed on the model matrix
\code{X} to obtain residuals. The phenotypic residuals are transformed
using \code{\link{RankNorm}}. In the next stage, the INT-transformed
residuals are regressed on the genotypic residuals.
}
\examples{
set.seed(100)
# Design matrix
X <- cbind(1, stats::rnorm(1e3))
# Genotypes
G <- replicate(1e3, stats::rbinom(n = 1e3, size = 2, prob = 0.25))
storage.mode(G) <- "numeric"
# Phenotype
y <- exp(as.numeric(X \%*\% c(1,1)) + stats::rnorm(1e3))
# Association test
p <- IINT(y = y, G = G, X = X)
}
\seealso{
\itemize{
  \item Basic association test \code{\link{BAT}}.
  \item Direct INT test \code{\link{DINT}}.
  \item Omnibus INT test \code{\link{OINT}}.
}
}
