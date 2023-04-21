#'@title Density of non-central squared copula
#'
#'@description This function computes the density of the  non-central squared copula (ncs)  associated with a one-parameter copula with parameter cpar, and parameters a1, a2 >0 .
#'@param data  Matrix (x,y) of size n x 2
#'@param family Copula family: "ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett''.
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par vector of copula parameter  and non-centrality parameter a1,a2 >0
#'
#'@return \item{pdf}{Density}
#'
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#'dncs(c(0.5,0.8),"ncs-clayton",par=c(2,1,2))
#'
#'
#'@export
#'
#'


dncs=function(data,  family,  rotation=0, par){
  cpar1=c("ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett")
  if(family %in% cpar1)
  {
  if(is.vector(data)){data = matrix(data,ncol=2)}
  u = data[,1]
  v = data[,2]

  if(rotation ==0){return(pdf(u,v,family,par))}
  else if(rotation==180){return(pdf(1-u,1-v,family,par))}
  else if(rotation==90){return (pdf(1-u,v,family,par))}
  else{ return (pdf(u,1-v,family,par))    }

  }else{warning("Family not available")
    return(NULL)}
}


pdf =function(u1,u2,family,par){

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

  w0  = dnorm(ha);
  U11 = cbind(htilde[,1], htilde[,3])
  U12 = cbind(htilde[,1], htilde[,4])
  U21 = cbind(htilde[,2], htilde[,3])
  U22 = cbind(htilde[,2], htilde[,4])

  w0011  = w0[,1] /(w0[,1]+w0[,2]);
  w0012  = 1-w0011;
  w0021 = w0[,3] /(w0[,3]+w0[,4]);
  w0022  = 1-w0021;


  w11 = w0011 * w0021;
  w12 = w0011 * w0022;
  w21 = w0012 * w0021;
  w22 = w0012 * w0022;

  if(family0=="plackett")
  {
    dCop = function(u) {dplac(u,rotation=0,cpar);}
  }else{
    dCop = function(u) {rvinecopulib::dbicop(u,family0,rotation=0,cpar);}
  }


  f11 = dCop(U11)+ 1E-20
  f12 = dCop(U12)+ 1E-20
  f21 = dCop(U21)+ 1E-20
  f22 = dCop(U22)+ 1E-20

  return(w11 * f11 + w12 * f12 + w21 * f21 + w22 * f22)

}
