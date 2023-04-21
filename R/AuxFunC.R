#'@title Auxiliary functions using C
#'
#'@description This function computes the empirical margins, their left-limits, Kendall's tau and Spearman's rho for arbitrary data
#'@param data  Matrix (x,y) of size n x 2
#'
#'
#'@return \item{tau}{Kendall's tau}
#'@return \item{rho}{Spearman's rho}
#'@return \item{Fx}{Empirical cdf of x}
#'@return \item{Fxm}{Left-limit of the empirical cdf of x}
#'@return \item{Fy}{Empirical cdf of y}
#'@return \item{Fym}{Left-limit of the empirical cdf of y}
#'
#'@references Nasri (2022). Test of serial dependence for arbitrary distributions. JMVA
#'@references Nasri & Remillard (2023). Tests of independence and randomness for arbitrary data using copula-based covariances, arXiv 2301.07267.
#'
#'@examples
#'data(simgumbel)
#'out=AuxFunC(simgumbel)
#'
#'
#'@export



AuxFunC =function(data){
  dim0=dim(data)
  n = dim0[1]
  d = dim0[2]

  x = data[,1]
  y = data[,2]
  out0 = .C("estdep",
            as.double(x),
            as.double(y),
            as.integer(n),
            tau = double(1),
            rho = double(1),
            s2 = double(1),
            Fx = double(n),
            Fxm = double(n),
            Fy = double(n),
            Fym = double(n),
            Ix = integer(n),
            Iy = integer(n),
            PACKAGE = "CopulaInference"

  )
  return(out0[-c(1:3)])
}
