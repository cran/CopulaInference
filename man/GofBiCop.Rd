% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GofBiCop.R
\name{GofBiCop}
\alias{GofBiCop}
\title{Goodness-of-fit for bivariate copula-based models with arbitrary distributions}
\usage{
GofBiCop(
  data = NULL,
  family,
  rotation = 0,
  Fx = NULL,
  Fxm = NULL,
  Fy = NULL,
  Fym = NULL,
  B = 100,
  n_cores = 1
)
}
\arguments{
\item{data}{Matrix or data frame with 2 columns (X,Y). Can be pseudo-observations. If NULL, Fx and Fy must be provided.}

\item{family}{Copula family: "gaussian", "t", "clayton", "frank", "gumbel", "joe", "plackett'', "bb1", "bb6", "bb7","bb8","ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett".}

\item{rotation}{Rotation: 0 (default value), 90, 180, or 270.}

\item{Fx}{marginal cdf function applied to X (default is NULL).}

\item{Fxm}{left limit of marginal cdf  function applied to X default is NULL).}

\item{Fy}{marginal cdf function  applied to  Y (default is NULL).}

\item{Fym}{left limit of marginal cdf function applied to  Y (default is NULL).}

\item{B}{Number of bootstrap samples (default 100)}

\item{n_cores}{Number of cores to be used for parallel computing (default is 1).}
}
\value{
\item{pvalueSn}{Pvalue of Sn in percent}

\item{pvalueTn}{Pvalue of Tn in percent}

\item{pvalueRn}{Pvalue of Rn in percent}

\item{Sn}{Value of Cramer-von Mises statistic Sn}

\item{Tn}{Value of Kendall's statistic Tn}

\item{Rn}{Value of Spearman's statistic Rn}

\item{cpar}{Copula parameters}

\item{family}{Copula family}

\item{rotation}{Rotation value}

\item{tauth}{Kendall's tau (from the multilinear theoretical copula)}

\item{tauemp}{Empirical Kendall's tau (from the multilinear empirical copula)}

\item{rhoth}{Spearman's rho (from the multilinear theoretical copula)}

\item{rhoemp}{Empirical Spearman's rho (from the multilinear empirical copula)}

\item{parB}{Bootstrapped parameters}

\item{loglik}{Log-likelihood}

\item{aic}{AIC value}

\item{bic}{BIC value}
}
\description{
Goodness-of-fit tests for copula-based models for data with arbitrary distributions. The tests statistics are the Cramer-von Mises statistic (Sn), the difference between the empirical Kendall's tau and the theoretical one, and the difference between the empirical Spearman's rho and the theoretical one.
}
\examples{
data = rvinecopulib::rbicop(10,"gumbel",rotation=0,2)
out=GofBiCop(data,family="gumbel",B=10)


}
\references{
Nasri & Remillard (2023). Identifiability and inference for copula-based semiparametric models for random vectors with arbitrary marginal distributions. arXiv 2301.13408.

Nasri & Remillard (2023). Goodness-of-fit and bootstrapping for copula-based random vectors with arbitrary marginal distributions.

Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
}
