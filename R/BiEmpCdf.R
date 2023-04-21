#'@title Empirical bivariate cdf
#'
#'@description This function computes the empirical joint cdf evaluated at all points (y1,y2)
#'@param data  Matrix (x1,x2) of size n x 2
#'@param y1  Vector of size n1
#'@param y2  Vector of size n2
#'@return \item{cdf}{Empirical cdf}
#'
#'@examples
#'data(simgumbel)
#'out=BiEmpCdf(simgumbel,c(0,1),c(-1,0,1))
#'
#'
#'@export
#'
#'


BiEmpCdf =function(data,y1,y2){
  dim0=dim(data)
  n = dim0[1]
  n1 = length(y1)
  n2 = length(y2)
  x1 = data[,1]
  x2 = data[,2]
  out0 = .C("bi_emp_cdf",
            as.double(x1),
            as.double(x2),
            as.integer(n),
            as.integer(n1),
            as.integer(n2),
            as.double(y1),
            as.double(y2),
            cdf = double(n1*n2),
            PACKAGE = "CopulaInference"

  )
  return(out0)
}
