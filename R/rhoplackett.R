#'@title Spearman's rho for Plackett copula
#'
#'@description Computes the theoretical Spearman's rho for Plackett copula
#'@param cpar  Copula parameter; can be a vector.
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'
#'
#'@return \item{rho}{Spearman's rho}
#'@references Remillard (2013). Statistical Methods for Financial Engineering. CRC Press
#'
#'@examples
#'rhoplackett(3,rotation=90)
#'
#'
#'@export
#'
#'
rhoplackett = function(cpar,rotation=0){

  d=length(cpar)
  ind = (cpar!=1)
  rho=rep(0,d)
  cpar1 = cpar[ind]
  rho[ind] = (cpar1+1)/(cpar1-1)-2*cpar1*log(cpar1)/(cpar1-1)^2
  if(rotation ==90 || rotation==270){rho=-rho}
  return(rho)
}
