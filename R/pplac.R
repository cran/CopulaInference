#'@title Cdf for Plackett copula
#'
#'@description This function computes the distribution function of the Plackett copula with parameter par>0.
#'@param data  Matrix (x,y) of size n x 2
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par Copula parameter >0
#'
#'
#'@return \item{cdf}{Value of cdf}
#'
#'@examples
#'pplac(c(0.5,0.8),270,3)
#'
#'
#'@export
#'
#'


pplac=function(data,rotation=0, par){


  if(is.vector(data)){data = matrix(data,ncol=2)}
  u = data[,1]
  v = data[,2]

  dim0=dim(data)
  n = dim0[1]


  if(rotation ==0||rotation==180){return(cdfplackett(u,v,par))}
  else if(rotation==90){
    return (v-cdfplackett(1-u,v,par))}
  else{
    return (u-cdfplackett(u,1-v,par))
    }
  }


cdfplackett =function(u,v,par){


 n = length(u);


out0 = .C("ppla",
          as.double(u),
          as.double(v),
          as.integer(n),
          as.double(par),
          cdf= double(n),
          PACKAGE = "CopulaInference"
         )
out0$cdf
}
