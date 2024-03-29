% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DINT.R
\name{DINT}
\alias{DINT}
\title{Direct-INT}
\usage{
DINT(
  y,
  G,
  X = NULL,
  k = 0.375,
  test = "Score",
  ties.method = "average",
  simple = FALSE
)
}
\arguments{
\item{y}{Numeric phenotype vector.}

\item{G}{Genotype matrix with observations as rows, SNPs as columns.}

\item{X}{Model matrix of covariates and structure adjustments. Should include
an intercept. Omit to perform marginal tests of association.}

\item{k}{Offset applied during rank-normalization. See
\code{\link{RankNorm}}.}

\item{test}{Either Score or Wald.}

\item{ties.method}{Method of breaking ties, passed to \code{base::rank}.}

\item{simple}{Return the p-values only?}
}
\value{
If \code{simple = TRUE}, returns a vector of p-values, one for each column
  of \code{G}. If \code{simple = FALSE}, returns a numeric matrix, including the
  Wald or Score statistic, its standard error, the Z-score, and the p-value.
}
\description{
Applies the rank-based inverse normal transformation (\code{\link{RankNorm}})
to the phenotype \code{y}. Conducts tests of association between the loci in
\code{G} and transformed phenotype, adjusting for the model matrix \code{X}.
}
\examples{
set.seed(100)
# Design matrix
X <- cbind(1, stats::rnorm(1e3))
# Genotypes
G <- replicate(1e3, stats::rbinom(n = 1e3, size = 2, prob = 0.25))
storage.mode(G) <- "numeric"
# Phenotype
y <- exp(as.numeric(X \%*\% c(1, 1)) + stats::rnorm(1e3))
# Association test
p <- DINT(y = y, G = G, X = X)
}
\seealso{
\itemize{
  \item Basic association test \code{\link{BAT}}.
  \item Indirect INT test \code{\link{IINT}}.
  \item Omnibus INT test \code{\link{OINT}}.
}
}
