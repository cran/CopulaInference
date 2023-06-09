% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pplac.R
\name{pplac}
\alias{pplac}
\title{Cdf for Plackett copula}
\usage{
pplac(data, rotation = 0, par)
}
\arguments{
\item{data}{Matrix (x,y) of size n x 2}

\item{rotation}{Rotation: 0 (default value), 90, 180, or 270.}

\item{par}{Copula parameter >0}
}
\value{
\item{cdf}{Value of cdf}
}
\description{
This function computes the distribution function of the Plackett copula with parameter par>0.
}
\examples{
pplac(c(0.5,0.8),270,3)


}
