#'@title Goodness-of-fit statistics
#'
#'@description Computation of goodness-of-fit statistics (Cramer-von Mises and the Kendall's tau)
#'@param object Object of class `EstBiCop`.
#'
#'@return \item{Sn}{Cramer-von Mises statistic}
#'@return \item{Tn}{Kendall's statistic}
#'@return \item{Rn}{Spearman's statistic}
#'@return \item{tauemp}{Empirical Kendall's tau}
#'@return \item{tauth}{Kendall's tau of the multilineat theoretical copula}
#'@return \item{rhoemp}{Empirical Spearman's rho}
#'@return \item{rhoth}{Spearman's rho of the multilineat theoretical copula}
#'@return \item{Y1}{Ordered observed values of X1}
#'@return \item{F1}{Empirical cdf of Y1}
#'@return \item{Y2}{Ordered observed values of X2}
#'@return \item{F2}{Empirical cdf of Y2}
#'@return \item{cpar}{Copula parameters}
#'@return \item{family}{Copula family}
#'@return \item{rotation}{Rotation value}
#'@return \item{n}{Sample size}
#'
#'@references Nasri & Remillard (2023). Identifiability and inference for copula-based semiparametric models for random vectors with arbitrary marginal distributions. arXiv 2301.13408.
#'
#'@examples
#' set.seed(2)
#' data = matrix(rpois(20,1),ncol=2)
#' out0 = EstBiCop(data,"gumbel")
#' out = statcvm(out0)
#'
#'
#'@export
#'
#'

statcvm = function(object)
{

ncs_fam =c("ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett")

data = object$data
X1 = data[,1]
X2 = data[,2]
y1 = sort(unique(X1))
y2 = sort(unique(X2))
F1 = sort(unique(object$F1))
F1m = sort(unique(object$F1m))
F2  = sort(unique(object$F2))
F2m = sort(unique(object$F2m))
# fun = function(x){max(0,x)}
# out1 = preparedata(X1)
# out2 = preparedata(X2)


f1 = F1-F1m
f2 = F2-F2m

n1 = length(y1)
if(F1[n1]<1){F1[n1]=1}
n2 = length(y2)
if(F2[n2]<1){F2[n2]=1}

y1m = c(y1[1]-1,y1[1:(n1-1)])
y2m = c(y2[1]-1,y2[1:(n2-1)])

f1rep = rep(f1,n2)
f2rep = rep(f2,each=n1)
f12 = f1rep*f2rep

MempA = BiEmpCdf(data,y1,y2)$cdf
MempB = BiEmpCdf(data,y1,y2m)$cdf
MempC = BiEmpCdf(data,y1m,y2)$cdf
MempD = BiEmpCdf(data,y1m,y2m)$cdf

hemp = MempA-MempB-MempC+MempD

tauemp = -1+sum(hemp*(MempA+MempB+MempC+MempD))
rhoemp = -3+3*sum(f12*(MempA+MempB+MempC+MempD))

a1 = rep(F1,n2)
a2 = rep(F2,each=n1)
aa = cbind(a1,a2)
b1 = rep(F1,n2)
b2 = rep(F2m,each=n1)
bb = cbind(b1,b2)
c1 = rep(F1m,n2)
c2 = rep(F2,each=n1)
cc = cbind(c1,c2)
d1 = rep(F1m,n2)
d2 = rep(F2m,each=n1)
dd = cbind(d1,d2)

family   = object$family
rotation = object$rotation
par      = object$par


if(family=="plackett")
{
  Cop = function(u) {pplac(u,rotation,par);}
  }else if(family %in% ncs_fam){
    Cop = function(u) {pncs(u,family,rotation,par)}
}else{
  Cop = function(u) {rvinecopulib::pbicop(u,family,rotation,par);}
}

#large = 1e10;

MthA = Cop(aa)
MthB = Cop(bb)
MthC = Cop(cc)
MthD = Cop(dd)

hth = MthA-MthB-MthC+MthD

tauth = -1+sum(hth*(MthA+MthB+MthC+MthD))
rhoth = -3+3*sum(f12*(MthA+MthB+MthC+MthD))

MA = MempA-MthA;
MB = MempB-MthB;
MC = MempC-MthC;
MD = MempD-MthD;


V = ( (MA+MB)^2 +(MA+MC)^2 +(MD+MB)^2 +(MD+MC)^2 + MA*MD + MB*MC )/18


n = dim(data)[1]
Sn = n*sum(f12*V) # cvm
Tn = sqrt(n)*(tauemp-tauth) #Kendall
Rn = sqrt(n)*(rhoemp-rhoth)  #Spearman
out0=object
out0$Sn=Sn
out0$Tn=Tn
out0$Rn=Rn
out0$tauemp=tauemp
out0$rhoemp=rhoemp
out0$tauth=tauth
out0$rhoth=rhoth
out0$Y1 = y1
out0$F1 = F1
out0$Y2 = y2
out0$F2 = F2
out0$n=n
class(out0) <- "statcvm"
return(out0)
}
