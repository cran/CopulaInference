#'@title Kendall's tau and Spearman's rho
#'
#'@description This function computes Kendall's tau and Spearman's rho for arbitrary data. These are invariant by increasing mappings.
#'@param data  Matrix or data frame with 2 columns (X,Y). Can be pseudo-observations.
#'
#'@return \item{tau}{Kendall's tau}
#'@return \item{rho}{Spearman's rho}
#'
#'@references Nasri (2022). Test of serial dependence for arbitrary distributions. JMVA
#'@references Nasri & Remillard (2023). Tests of independence and randomness for arbitrary data using copula-based covariances, arXiv 2301.07267.
#'
#'@examples
#'data(simgumbel)
#'out=EstDep(simgumbel)
#'
#'
#'@export



EstDep =function(data){
  n = dim(data)[1]

  x = data[,1]
  y = data[,2]

  out0 = .C("est_dep",
            as.double(x),
            as.double(y),
            as.integer(n),
            tau = double(1),
            rho = double(1),
            PACKAGE = "CopulaInference"

  )
  return(out0[-c(1:3)])
}
