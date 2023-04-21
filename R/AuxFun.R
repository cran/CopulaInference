#'@title Auxiliary functions
#'
#'@description This function computes the empirical margins, their left-limits, Kendall's tau and Spearman's rho for arbitrary data. Slower than AuxFunC based on C.
#'@param data  Matrix (x,y) of size n x 2
#'
#'
#'@return \item{tau}{Kendall's tau}
#'@return \item{rho}{Spearman's rho}
#'@return \item{Fx}{Empirical cdf of x}
#'@return \item{Fxm}{Left-limit of the empiricial cdf of x}
#'@return \item{Fy}{Empirical cdf of y}
#'@return \item{Fym}{Left-limit of the empiricial cdf of y}
#'
#'@references Nasri (2022). Test of serial dependence for arbitrary distributions. JMVA
#'@references Nasri & Remillard (2023). Tests of independence and randomness for arbitrary data using copula-based covariances, arXiv 2301.07267.
#'@examples
#'data(simgumbel)
#'out=AuxFun(simgumbel)
#'
#'
#'@export
#'
#'




AuxFun = function(data){
      if(is.vector(data)){data=matrix(data,ncol=1)}

          fun1 = function(x,y){ as.numeric(x<=y)}
          fun1m = function(x,y){ as.numeric(x<y)}

       dim0 = dim(data)
       n = dim0[1]
       d = dim0[2]
       Fn  = matrix(0,ncol=d,nrow=n)
       Fnm = matrix(0,ncol=d,nrow=n)
       Vn  = list()
       Vnm = list()
       tau = matrix(0,ncol=d,nrow=d)
       rho = matrix(0,ncol=d,nrow=d)
       for(j in 1:d)
       { x = data[,j]
         Vn[[j]]  = outer(x,x,FUN=fun1)
         Vnm[[j]] = outer(x,x,FUN=fun1m)
         Fn[,j]  = colSums(Vn[[j]])/(n+1)
         Fnm[,j] = colSums(Vnm[[j]])/(n+1)
       }

       if(d==1){tau=NULL; rho = NULL}else{
       for(j in 1:(d-1))
           {
             for(k in ((j+1):d))
               {
               c1 = -1+ colSums( (Vn[[j]] +Vnm[[j]]))/n
               c2 = -1+ colSums( (Vn[[k]] +Vnm[[k]]))/n
               s1 = mean(c1*c1)
               s2 = mean(c2*c2)
               s = sqrt(s1*s2)
               rho[j,k] = mean(c1*c2)/s
               rho[k,j] = rho[j,k]
                tau[j,k] = -1+ mean( mean( (Vn[[j]] +Vnm[[j]]) * (Vn[[k]]+Vnm[[k]]) ) );
                tau[k,j] = tau[j,k]
               }
           tau[j,j]=1
           rho[j,j]=1
          }
}
       out=list(tau=tau,rho=rho, Fn=Fn,Fnm=Fnm,n=n,d=d)
       out

}
