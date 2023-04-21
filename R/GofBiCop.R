#'@title Goodness-of-fit for bivariate copula-based models with arbitrary distributions
#'@description  Goodness-of-fit tests for copula-based models for data with arbitrary distributions. The tests statistics are the Cramer-von Mises statistic (Sn), the difference between the empirical Kendall's tau and the theoretical one, and the difference between the empirical Spearman's rho and the theoretical one.
#'@param data  Matrix or data frame with 2 columns (X,Y). Can be pseudo-observations. If NULL, Fx and Fy must be provided.
#'@param family Copula family: "gaussian", "t", "clayton", "frank", "gumbel", "joe", "plackett'', "bb1", "bb6", "bb7","bb8","ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett".
#'@param rotation Rotation: 0 (default value), 90, 180, or 270.
#'@param Fx      marginal cdf function applied to X (default is NULL).
#'@param Fxm     left limit of marginal cdf  function applied to X default is NULL).
#'@param Fy      marginal cdf function  applied to  Y (default is NULL).
#'@param Fym     left limit of marginal cdf function applied to  Y (default is NULL).
#'@param B Number of bootstrap samples (default 100)
#'@param n_cores Number of cores to be used for parallel computing (default is 1).
#'
#'@return \item{pvalueSn}{Pvalue of Sn in percent}
#'@return \item{pvalueTn}{Pvalue of Tn in percent}
#'@return \item{pvalueRn}{Pvalue of Rn in percent}
#'@return \item{Sn}{Value of Cramer-von Mises statistic Sn}
#'@return \item{Tn}{Value of Kendall's statistic Tn}
#'@return \item{Rn}{Value of Spearman's statistic Rn}
#'@return \item{cpar}{Copula parameters}
#'@return \item{family}{Copula family}
#'@return \item{rotation}{Rotation value}
#'@return \item{tauth}{Kendall's tau (from the multilinear theoretical copula)}
#'@return \item{tauemp}{Empirical Kendall's tau (from the multilinear empirical copula)}
#'@return \item{rhoth}{Spearman's rho (from the multilinear theoretical copula)}
#'@return \item{rhoemp}{Empirical Spearman's rho (from the multilinear empirical copula)}
#'@return \item{parB}{Bootstrapped parameters}
#'@return \item{loglik}{Log-likelihood}
#'@return \item{aic}{AIC value}
#'@return \item{bic}{BIC value}
#'
#'@references Nasri & Remillard (2023). Identifiability and inference for copula-based semiparametric models for random vectors with arbitrary marginal distributions. arXiv 2301.13408.
#'@references Nasri & Remillard (2023). Goodness-of-fit and bootstrapping for copula-based random vectors with arbitrary marginal distributions.
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#'data = rvinecopulib::rbicop(10,"gumbel",rotation=0,2)
#'out=GofBiCop(data,family="gumbel",B=10)
#'
#'
#'@export
#'
#'

GofBiCop = function(data=NULL, family, rotation=0,
                    Fx=NULL, Fxm=NULL, Fy=NULL, Fym=NULL,
                    B=100, n_cores=1)
{

object= EstBiCop(data,family,rotation,Fx,Fxm,Fy,Fym)


out0= statcvm(object)
Sn = out0$Sn
Tn = out0$Tn
Rn = out0$Rn
par = object$par

SnB = rep(0,B)
TnB = SnB
RnB=TnB
parB = matrix(0,nrow=B,ncol=length(par))

cl <- parallel::makePSOCKcluster(n_cores)
doParallel::registerDoParallel(cl)

result <- foreach::foreach(i=1:B, .packages="CopulaInference") %dopar% bootstrapfun(out0)

for (i in 1:B){
  SnB[i] = result[[i]]$Sn
  TnB[i] = result[[i]]$Tn
  RnB[i] = result[[i]]$Rn
  parB[i,] = result[[i]]$par
}

pvalSn = 100*mean(SnB >= Sn)
pvalTn = 100*mean(abs(TnB) >= abs(Tn))
pvalRn = 100*mean(abs(RnB) >= abs(Rn))

statB =list(Sn=B,Tn=TnB,RnB)
out = list(pvalueSn=pvalSn, pvalueTn=pvalTn,  pvalueRn=pvalRn, Sn=Sn, Tn=Tn, Rn=Rn, par=par,
           family=family, rotation=rotation, tauth = out0$tauth,
           tauemp = out0$tauemp, rhoth = out0$rhoth,
           rhoemp = out0$rhoemp, statB=statB, parB=parB, loglik=object$loglik,
           aic=object$aic, bic=object$bic)


parallel::stopCluster(cl)
out
}

