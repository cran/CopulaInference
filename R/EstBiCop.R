#'@title Parameter estimation for bivariate copula-based models with arbitrary distributions
#'
#'@description Computes the estimation of the parameters of a copula-based model with arbitrary distributions, i.e, possibly mixtures of discrete and continuous distributions. Parametric margins are allowed. The estimation is based on a pseudo-likelihood adapted to ties.
#'@param data  Matrix or data frame with 2 columns (X,Y). Can be pseudo-observations. If NULL, Fx and Fy must be provided.
#'@param family Copula family: "gaussian", "t", "clayton", "frank", "gumbel", "joe", "plackett'', "bb1", "bb6", "bb7","bb8","ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett".
#'@param rotation Rotation: 0 (default value), 90, 180, or 270.
#'@param Fx      Marginal cdf function applied to X (default is NULL).
#'@param Fxm     Left-limit of marginal cdf  function applied to X default is NULL).
#'@param Fy      Marginal cdf function  applied to  Y (default is NULL).
#'@param Fym     Left-limit of marginal cdf function applied to  Y (default is NULL).
#'
#'@return \item{par}{Copula parameters}
#'@return \item{family}{Copula family}
#'@return \item{rotation}{Rotation value}
#'@return \item{tauth}{Kendall's tau corresponding to the estimated parameter}
#'@return \item{tauemp}{Empirical Kendall's tau (from the multilinear empirical copula)}
#'@return \item{rhoSth}{Spearman's rho corresponding to the estimated parameter}
#'@return \item{rhoSemp}{Empirical Spearman's tau (from the multilinear empirical copula)}
#'@return \item{loglik}{Log-likelihood}
#'@return \item{aic}{Aic value}
#'@return \item{bic}{Bic value}
#'@return \item{data}{Matrix of values (could be (Fx,Fy))}
#'@return \item{F1}{Cdf of X (Fx if provided, empirical otherwise)}
#'@return \item{F1m}{Left-limit of F1 (Fxm if provided, empirical otherwise)}
#'@return \item{F2}{Cdf of Y (Fy if provided, empirical otherwise)}
#'@return \item{F2m}{Left-limit of F2 (Fym if provided, empirical otherwise)}
#'@return \item{ccdfx}{Conditional cdf of X given Y and it left limit}
#'@return \item{ccdfxm}{Left-limit of ccdfx}
#'@return \item{ccdfy}{Conditional cdf of Y given X and it left limit}
#'@return \item{ccdfym}{Left-limit of ccdfy}
#'
#'@references Nasri & Remillard (2023). Identifiability and inference for copula-based semiparametric models for random vectors with arbitrary marginal distributions. arXiv 2301.13408.
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#'set.seed(2)
#'data = matrix(rpois(20,1),ncol=2)
#'out0=EstBiCop(data,"gumbel")
#'
#'@export
#'
#'

EstBiCop = function(data=NULL, family,  rotation = 0,
                    Fx = NULL, Fxm = NULL, Fy = NULL, Fym = NULL ){

  ncs_fam =c("ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett")
  signe=1

  if(rotation ==90 || rotation ==270){signe = -1}



if(family %in% c("gaussian","t","frank")){rotation=0}

if(is.null(data))
  {
   data=cbind(Fx,Fy)
   F1 = Fx
   if(is.null(Fxm)){
    F1m=F1
    warning("Fx margin assumed to be continuous")
    I1 = rep(0,n)
    }else{
    F1m = Fxm
    I1 = as.numeric( ( (F1-F1m)>0 ) )
  }

  F2 = Fy
  if(is.null(Fym)){
    F2m=F2
    I2 = rep(0,n)
    warning("Fy margin assumed to be continuous")
  }else{
    F2m = Fym
    I2 = as.numeric( ( (F2-F2m)>0 ) )
  }
  out0 =EstDep(data)
}else{

      out0 = AuxFunC(data); #much faster than AuxFun

       # X margins
       if(is.null(Fx)){
        F1  = out0$Fx
        F1m = out0$Fxm
        I1  = out0$Ix
        }else{
          F1 = Fx
          if(is.null(Fxm)){
          F1m=F1
          I1 = rep(0,n)
           }else{
                 F1m = Fxm
                 I1 = as.numeric( ( (F1-F1m)>0 ) )
                }

            }
# Y margins
            if(is.null(Fy)){
              F2  = out0$Fy
              F2m = out0$Fym
              I2  = out0$Iy

            }else{
                  F2 = Fy
                  if(is.null(Fym)){
                     F2m=F2
                     I2 = rep(0,n)
                  }else{
                    F2m = Fym
                    I2 = as.numeric( ( (F2-F2m)>0 ) )
                    }

}


}

n = dim(data)[1]
f1 = I1*(F1-F1m)+1-I1  # always greater than 0!
f2 = I2*(F2-F2m)+1-I2  # always greater than 0!
F11  = cbind(F1, F2)
#F22  = cbind(F1m, F2m)+1E-40;
F22  = cbind(F1m, F2m)
tau  = out0$tau
rhoS = out0$rho

rhoSth=NULL


F12 = cbind(F11[,1],F22[,2]) # x1 is continuous, not x2
F21 = cbind(F22[,1],F11[,2]) # x2 is continuous, not x1

out00 = est_options(family,tau)

LB = out00$LB
UB = out00$UB
start=out00$start


likf = function(par)
{


  if(min(par - LB) < 0 || max(par - UB) > 0) return(1e5)

  if(family=="plackett")
    {
     Cop = function(u) {pplac(u,rotation, par);}
    }else if(family %in% ncs_fam){
      Cop = function(u) {pncs(u,family,rotation, par)}}else{
      Cop = function(u) {rvinecopulib::pbicop(u,family,rotation,par);}
          }

  z = Cop(F11)- Cop(F12)- Cop(F21) + Cop(F22);

  if(family=="plackett")
    {
      C1  =  function(u) {hplac(u,1,rotation, par);}  # x1 is continuous, not x2
      C2  =  function(u) {hplac(u,2,rotation, par);}  # x2 is continuous, not x1
      C12 =  function(u) {dplac(u,rotation, par);}
    }else if(family %in% ncs_fam){
      C1  =  function(u) {hncs(u,1,family,rotation, par);}  # x1 is continuous, not x2
      C2  =  function(u) {hncs(u,2,family,rotation, par);}  # x2 is continuous, not x1
      C12 =  function(u) {dncs(u,family,rotation, par);}
      }else{
         C1  =  function(u) {rvinecopulib::hbicop(u,1,family,rotation,par);}  # x1 is continuous, not x2
         C2  =  function(u) {rvinecopulib::hbicop(u,2,family,rotation,par);}  # x2 is continuous, not x1
         C12 =  function(u) {rvinecopulib::dbicop(u,family,rotation,par);}  # x1 and x2 are both continuous
  }

  p1  = C1(F11) - C1(F12);   # x1 is continuous, not x2
  p2  = C2(F11) - C2(F21);   # x2 is continuous, not x1
  p12 = C12(F11);             # x1 and x2 are both continuous


  z = I1 * I2 * z + I1 * (1-I2) * p2 + (1-I1) * I2 * p1  + (1-I1) * (1-I2) * p12 ;

  LL = log( z + 1E-20);
  return(-sum(LL))
}

mle = nlm(likf,p=start,check.analyticals=F,print.level=2,iterlim=400,hessian=TRUE)

cpar = mle$estimate
InfMat= mle$hessian/n
LL = mle$minimum
#covar = matlib::inv(InfMat)
npar=length(cpar)
aic =  2*npar+2*LL
bic = npar*log(n)+2*LL

if(family=="plackett")
{
  tauth=tauplackett(cpar[1],rotation)
  }else if(family %in% ncs_fam){
  family0=gsub('ncs-', '', family)
  fnum = fnumber(family0)
  if(fnum<7){ tauth = taucop(fnum, cpar[1],rotation);}
  else{tauth = taucop(fnum, rotation, cpar[1])}
  }else
  {
   # fnum = fnumber(family,rotation)
   # tauth = VineCopula::BiCopPar2Tau(fnum, cpar[1],cpar[2])
   fnum = fnumber(family)
   if(fnum<7){ tauth = taucop(fnum, cpar,rotation);}
   else{tauth = taucop(fnum, rotation, cpar);}
  }



if(family=="plackett"){rhoSth = rhoplackett(cpar[1],rotation)}

##### Conditional distributions


if(family=="plackett")
{
  Cop = function(u) {pplac(u,rotation, cpar);}
}else if(family %in% ncs_fam){
  Cop = function(u) {pncs(u,family,rotation, cpar)}}else{
    Cop = function(u) {rvinecopulib::pbicop(u,family,rotation,cpar);}
  }


D1  = (Cop(F11) - Cop(F12)) /f2   #cond with respect to Y
D1m = (Cop(F21) - Cop(F22) )/f2   #cond with respect to Y
D2  = (Cop(F11) - Cop(F21) )/f1;  #cond with respect to X
D2m = (Cop(F12) - Cop(F22) )/f1;  #cond with respect to X

if(family=="plackett")
{
  C1  =  function(u) {hplac(u,1,rotation, cpar);}
  C2  =  function(u) {hplac(u,2,rotation, cpar);}

}else if(family %in% ncs_fam){
  C1  =  function(u) {hncs(u,1,family,rotation, cpar);}
  C2  =  function(u) {hncs(u,2,family,rotation, cpar);}

}else{
  C1  =  function(u) {rvinecopulib::hbicop(u,1,family,rotation,cpar);}
  C2  =  function(u) {rvinecopulib::hbicop(u,2,family,rotation,cpar);}

}


p1  = C2(F11)     # x2 is continuous
p1m = C2(F21)     # x2 is continuous


p2  = C1(F11)     # x1 is continuous
p2m = C1(F12)     # x1 is continuous

ccdfx  = I2* D1 + (1-I2)*p1
ccdfxm = I2* D1m + (1-I2)*p1m

ccdfy  = I1* D2 + (1-I1)*p2
ccdfym = I1* D2m + (1-I1)*p2m

##############################


out = list(par=cpar, family=family, rotation=rotation, tauth=tauth,
           rhoSth=rhoSth, rhoSemp=rhoS, tauemp = tau, loglik=-LL,
           aic=aic,bic=bic, data=data, F1=F1, F1m=F1m, F2=F2,F2m=F2m,
           ccdfx=ccdfx,ccdfxm=ccdfxm,ccdfy=ccdfy, ccdfym=ccdfym)
class(out) <- "EstBiCop"
out
}




