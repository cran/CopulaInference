#'@title Density of Plackett copula
#'
#'@description This function computes the density of the Plackett copula with parameter par>0.
#'@param data  Matrix (x,y) of size n x 2
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par Copula parameter >0
#'
#'
#'@return \item{pdf}{Density}
#'
#'@examples
#'dplac(c(0.5,0.8),par=3,rotation=270)
#'
#'
#'@export
#'
#'


dplac=function(data,rotation=0,par){


  if(is.vector(data)){data = matrix(data,ncol=2)}
  u = data[,1]
  v = data[,2]

  dim0=dim(data)
  n = dim0[1]


  if(rotation ==0||rotation==180){return(pdfplackett(u,v,par))}
  else if(rotation==90){
    return (pdfplackett(1-u,v,par))}
  else{
    return (pdfplackett(u,1-v,par))
    }
  }


pdfplackett =function(u,v,par){


 n = length(u);


out0 = .C("dpla",
          as.double(u),
          as.double(v),
          as.integer(n),
          as.double(par),
          pdf= double(n),
          PACKAGE = "CopulaInference"
         )
out0$pdf
}
