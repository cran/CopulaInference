#'@title Kendall's tau for Plackettfamily
#'
#'@description This function computes Kendall's tau for Plackett family using numerical integration
#'@param cpar Copula parameter >0
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'
#'
#'@return \item{tau}{Kendall's tau}
#'
#'@examples
#'tauplackett(2,270)
#'
#'
#'@export
#'
#'


tauplackett = function(cpar, rotation=0)
{

out0 = .C("taupla",
          as.double(cpar),
          tau = double(1),
          PACKAGE = "CopulaInference"

)
tau = out0$tau
if(rotation==90||rotation==270){tau=-tau}
tau
}
