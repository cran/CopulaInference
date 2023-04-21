#'@title Simulation of non-central squared copula
#'
#'@description This function computes generates a bivariate sample from a non-central squared copula (ncs)  associated with a one-parameter copula with parameter cpar, and parameters a1, a2 >0 .
#'@param n  Number of observations
#'@param family Copula family: "ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett''.
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par vector of copula parameter  and non-centrality parameter a1,a2 >0

#'
#'
#'
#'@return \item{U}{Observations}
#'
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#'rncs(100,"ncs-clayton",par=c(2,1,2))
#'
#'
#'@export
#'
#'


rncs=function(n,  family,  rotation=0, par){
  cpar1=c("ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett")
  if(family %in% cpar1)
  {
  U = matrix(0,ncol=2,nrow=n)
  family0=gsub('ncs-', '', family)
  cpar=par[1]
  a1 = par[2]
  a2 = par[3]


  d1 = a1^2;
  d2 = a2^2;
  if(family0=="plackett"){
    V = rplac(n,rotation,cpar)}else{
    V = rvinecopulib::rbicop(n,family0,rotation,cpar)
  }
  X = (qnorm(V)+par[2:3])^2
  U1 = pchisq(X[,1],1,d1)
  U2 = pchisq(X[,2],1,d2)

  U[,1] <- U1
  U[,2] <- U2
  return (U)

  }else{warning("Family not available")
    return(NULL)}
}

