#'@title Kendall's tau for a copula family
#'
#'@description This function computes Kendall's tau for a copula family
#'@param family_number  Integer from 1 to 10
#'@param cpar Copula parameters
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'
#'
#'@return \item{tau}{Kendall's tau}
#'
#'@examples
#'taucop(4,2,270) # Gumbel copula
#'
#'
#'@export
#'
#'


taucop = function(family_number, cpar, rotation=0)
{
out0 = .C("taucopula",
          as.integer(family_number),
          as.integer(rotation),
          as.double(cpar),
          tau = double(1),
          PACKAGE = "CopulaInference"

)
out0$tau
}
