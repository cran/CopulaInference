#'@title Function to perform parametric bootstrap for goodness-of-fit tests
#'
#'@description This function simulates Cramer-von Mises statistic and Kendall's statistic using parametric bootstrap for arbitrary data.
#'@param object object of class `statcvm`.
#
#'
#'
#'@return \item{Sn}{Simulated value of the Cramer-von Mises statistic}
#'@return \item{Tn}{simulated value of the Kendall's statistic}
#'@return \item{Rn}{simulated value of the Spearman's statistic}
#'@return \item{parB}{Estimated parameter}
#'@references Nasri & Remillard (2023). Identifiability and inference for copula-based semiparametric models for random vectors with arbitrary marginal distributions. arXiv 2301.13408.
#'@references Nasri & Remillard (2023). Goodness-of-fit and bootstrapping for copula-based random vectors with arbitrary marginal distributions.
#'
#'
#'@keywords internal
#'
#'@export
#'
#'
bootstrapfun = function(object){
  ncs_fam =c("ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett")

  n        = object$n
  family   = object$family
  rotation = object$rotation
  par      = object$par

  if(family=="plackett")
  {
   U = rplac(n,rotation,par)
  }else if(family %in% ncs_fam){
    U= rncs(n,family,rotation,par)
    }else{
      U= rvinecopulib::rbicop(n,family,rotation,par)
  }

Y1 = object$Y1
Y2 = object$Y2
F1 = object$F1
F2 = object$F2

Bn1   = EmpCdf(U[,1])$Fx
Bn2   = EmpCdf(U[,2])$Fx
X1B   = CdfInv(Bn1,Y1,F1)
X2B   = CdfInv(Bn2,Y2,F2)
dataB = cbind(X1B,X2B)
out0  = EstBiCop(dataB,family,rotation)
outB  = statcvm(out0)


out = list(Sn=outB$Sn,Tn = outB$Tn, Rn=outB$Rn, parB=outB$par)
return(out)
}
