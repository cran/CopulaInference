#'@title Empirical univariate cdf
#'
#'@description This function computes the empirical cdf evaluated at all sample points
#'@param x  Observations
#'@return \item{Fx}{Empirical cdf}
#'@return \item{Fxm}{Left limit of the empirical cdf}
#'@return \item{Ix}{Indicator of atoms}data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
#'
#'@examples
#'data(simgumbel)
#'out=EmpCdf(simgumbel[,1])
#'
#'
#'@export
#'
#'


EmpCdf =function(x){
  n = length(x)

  out0 = .C("emp_cdf",
            as.double(x),
            as.integer(n),
            Fx  = double(n),
            Fxm = double(n),
            Ix  =integer(n),
            PACKAGE = "CopulaInference"

  )
  return(out0)
}
