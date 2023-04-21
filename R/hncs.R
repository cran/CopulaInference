#'@title Conditional distribution of non-central squared copula
#'
#'@description This function computes the conditional distribution of the  non-central squared copula (ncs)  associated with a one-parameter copula with parameter cpar, and parameters a1, a2 >0 .
#'@param data  Matrix (x,y) of size n x 2
#'@param cond_var Conditioning variable (1 or 2)
#'@param family Copula family: "ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett''.
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par vector of copula parameter  and non-centrality parameter a1,a2 >0
#'
#'
#'
#'@return \item{h}{Conditional  cdf}
#'
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#' hncs(c(0.5,0.8),1,"ncs-clayton",270,c(2,1,2))
#'
#'
#'@export
#'
#'


hncs=function(data, cond_var, family, rotation=0,par){
  cpar1=c("ncs-gaussian", "ncs-clayton", "ncs-frank", "ncs-gumbel", "ncs-joe", "ncs-plackett")
  if(family %in% cpar1)
  {
  if(is.vector(data)){data = matrix(data,ncol=2)}
  u = data[,1]
  v = data[,2]

  if(rotation ==0){return(ccdf(u,v,cond_var,family,par))}
  else if(rotation==90){
    if(cond_var==2){ return(1-ccdf(1-u,v,2,family,par))}else{
      return(ccdf(u,1-v,1,family,par))}}
  else if(rotation==180){
    return(1-ccdf(1-u,1-v,cond_var,family,par))}
  else if(rotation==270){
    if(cond_var==2){
      return(ccdf(u,1-v,2,family,par))
      }else{
      return(1-ccdf(u,1-v,1,family,par))}
  }
  }else{warning("Family not available")
    return(NULL)}
}


ccdf =function(u1,u2,cond_var,family,par){

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
  if(family0=="plackett")
  {
    CCop = function(u) {hplac(u,cond_var,rotation=0,cpar);}
  }else{
    CCop = function(u) {rvinecopulib::hbicop(u,cond_var,family0,rotation=0,cpar);}
  }


  G11 = CCop(U11)
  G12 = CCop(U12)
  G21 = CCop(U21)
  G22 = CCop(U22)

  if(cond_var==1){
    return(w0011 * (G11 - G12) + (1-w0011) * (G21-G22))}else{
    return(w0021 * (G11 - G21) + (1-w0021) * (G12-G22))}

}
