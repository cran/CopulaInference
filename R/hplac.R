#'@title Conditional distribution of Plackett copula
#'
#'@description This function computes the  conditional distribution of the Plackett copula with parameter par>0.
#'@param data  Matrix (x,y) of size n x 2
#'@param cond_var Conditioning variable (1 or 2)
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par Copula parameter >0
#'
#'
#'@return \item{h}{Conditional  cdf}
#'
#'@examples
#' hplac(c(0.5,0.8),1,270,3)
#'
#'
#'@export
#'
#'


hplac=function(data,cond_var,rotation=0,par){

  if(is.vector(data)){data = matrix(data,ncol=2)}
  u = data[,1]
  v = data[,2]

  if(rotation ==0){
     return( hplackett(u,v,cond_var,par))
    }else if(rotation==90){
    if(cond_var==2){
      return(1-hplackett(1-u,v,2,par))
      }else{
        return(hplackett(u,1-v,1,par))}
      }else if(rotation==180){
       return(1-hplackett(1-u,1-v,cond_var,par))
        }else if(rotation==270){
    if(cond_var==2){
      return(hplackett(u,1-v,2,par))}else{
       return(1-hplackett(u,1-v,1,par))}
        }
}


hplackett =function(u,v,cond_var,par){


 n = length(u);


out0 = .C("hpla",
          as.double(u),
          as.double(v),
          as.integer(n),
          as.integer(cond_var),
          as.double(par),
          cpdf= double(n),
          PACKAGE = "CopulaInference"
         )
return(out0$cpdf)
}
