#'@title Cdf for non-central squared copula
#'
#'@description This function computes the distribution function of the  non-central squared copula (ncs)  associated a with one-parameter copula with parameter cpar, and parameters a1, a2 >0 .
#'@param data  Matrix (x,y) of size n x 2
#'@param family Copula family: "ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett''.
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par vector of copula parameter  and non-centrality parameter a1,a2 >0
#'
#'
#'@return \item{cdf}{Value of cdf}
#'
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#'pncs(c(0.5,0.8),"ncs-clayton", par=c(2,1,2),rotation=270)
#'
#'
#'@export
#'
#'


pncs=function(data,family,rotation=0,par){

  cpar1=c("ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett")
  if(family %in% cpar1)
  {
   if(is.vector(data)){data = matrix(data,ncol=2)}
    u = data[,1]
    v = data[,2]

    dim0=dim(data)
    n = dim0[1]


   if(rotation ==0){return(cdf(u,v,family,par))}
    else if(rotation==180){return(cdf(1-u,1-v,family,par))}
    else if(rotation==90){return (v-cdf(1-u,v,family,par))}
    else{ return (u-cdf(u,1-v,family,par))    }
  }else{warning("Family not available")
    return(NULL)}
}

cdf =function(u1,u2,family,par){
  family0=gsub('ncs-', '', family)
  cpar=par[1]
  a1 = par[2]
  a2 = par[3]


  d1 = a1^2;
  d2 = a2^2;


  G1 = sqrt(qchisq(u1,1,d1))  #Ga inverse
  G2 = sqrt(qchisq(u2,1,d2))

  ha = cbind((G1-a1),  (-G1-a1) , (G2-a2) ,(-G2-a2))

  htilde = pnorm(ha)+ 1E-10;

  U11 = cbind(htilde[,1], htilde[,3])
  U12 = cbind(htilde[,1], htilde[,4])
  U21 = cbind(htilde[,2], htilde[,3])
  U22 = cbind(htilde[,2], htilde[,4])


if(family0=="plackett")
{
  Cop = function(u) {pplac(u,rotation=0,cpar);}
}else{
  Cop = function(u) {rvinecopulib::pbicop(u,family0,rotation=0,cpar);}
}


prob = Cop(U11)-Cop(U12)-Cop(U21) + Cop(U22)

return(prob)
}
