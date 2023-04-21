#'@title Quantile function
#'
#'@description This function computes the inverse of the cdf of a finite distribution for a vector of probabilities.
#'@param u   Vector of probabilities
#'@param y   Ordered values
#'@param Fn  Cdf
#'
#'@return \item{x}{Vector of quantiles}
#'
#'@examples
#'y=c(0,1,2)
#'Fn = c(0.5,0.85,1)
#'out=CdfInv(c(1:9)/10,y,Fn)
#'
#'
#'@export
#'
#'

CdfInv = function(u,y,Fn){
  n = length(u)
  x = rep(0,n)

  for(i in 1:n){
    k = min(which(Fn>=u[i]))
    x[i] = y[k];
  }
 return(x)
}

