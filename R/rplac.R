#'@title Generates observations from the Plackett copula
#'
#'@description This function generates observations from a Plackett copula with parameter par>0.
#'@param n  Number of pairs to be generated
#'@param rotation   Rotation: 0 (default value), 90, 180, or 270.
#'@param par Copula parameter >0
#'
#'
#'@return \item{U}{Matrix of observations}
#'
#'@examples
#'rplac(10,rotation=90,par=2)
#'
#'
#'@export
#'
#'
rplac <- function(n,rotation=0,par){


U1 <- runif(n)
V  <- runif(n)

a <- par
W <- V*(1-V)
A <- a+(a-1)^2*W


C <- W * (1+(a-1)*U1)^2
B <- 2*(a-1)*W * (1-(a+1)*U1)-a

U2 <- 0.5*(-B + sign(V-0.5)*sqrt(B^2-4*A*C))/A

U = matrix(0,n,2)

if(rotation==90){U1=1-U1}
else if(rotation==180){
  U1 = 1-U1
  U2 = 1-U2}
else if(rotation==270){U2=1-U2}

U[,1] <- U1
U[,2] <- U2
return (U)

}
