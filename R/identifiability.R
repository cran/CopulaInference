#'@title Identifiability of  two-parameter copula families
#'
#'@description Determines if a copula family is identifiable with respect to the empirical margins. One-parameter copula families ("gaussian","gumbel","clayton","frank","plackett","joe") are identifiable whatever the margins. The rank of the gradient of the copula on the range of the margins is evaluated at 10000 parameter points within the lower and upper bounds of the copula family.
#'@param data  Matrix or data frame with 2 columns (X,Y). Can be pseudo-observations. If NULL, Fx and Fy must be provided.
#'@param family Copula family: "gaussian", "t", "clayton", "frank", "gumbel", "joe", "plackett'', "bb1", "bb6", "bb7","bb8","ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett".
#'@param rotation Rotation: 0 (default value), 90, 180, or 270.
#'@param Fx      Marginal cdf function applied to X (default is NULL).
#'@param Fy      Marginal cdf function  applied to  Y (default is NULL).
#'
#'@return \item{out}{True or False}
#'
#'@references Nasri & Remillard (2023). Identifiability and inference for copula-based semiparametric models for random vectors with arbitrary marginal distributions. arXiv 2301.13408.
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'
#'@examples
#' set.seed(1)
#' data = matrix(rpois(20,1),ncol=2)
#' out = identifiability(data,"gumbel")
#'
#'
#'@export
#'
#'
identifiability = function(data=NULL, family,rotation=0, Fx=NULL, Fy=NULL)
{
  out = TRUE
  cpar2=c("t", "bb1","bb6","bb7","bb8")
  cpar3 = c("ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett")

  if(is.null(data))
  {
    F1 = sort(unique(Fx))
    F2 = sort(unique(Fy))
    n1 = length(F1)
    n2 = length(F2)
    if(F1[n1]==1){n1=n1-1}
    if(F2[n2]==1){n2=n2-1}

  }else{
    # X margins
    if(is.null(Fx)){
      X1 = data[,1]
      out1 = preparedata(X1)
      F1  = out1$Fn
      n1   = out1$m-1
      }else{F1 = Fx}
    # Y margins
    if(is.null(Fy)){
      X2 = data[,2]
      out2 = preparedata(X2)
      F2  = out2$Fn
      n2   = out2$m-1
      }else{F2 = Fy}
}




if(family %in% cpar2)
  {



    F1 = F1[1:n1]
    F2 = F2[1:n2]
    d = n1*n2

    u1 = rep(F1,n2)
    u2 = rep(F2,each=n1)
    u = cbind(u1,u2)
    Cop = function(par) {rvinecopulib::pbicop(u,family,rotation,par)}
    out0 = est_options(family)
    UB = out0$UB
    LB = out0$LB

    b1 = UB[1]-.001
    b2 = UB[2]-.001
    a1 = LB[1]+.001
    a2 = LB[2]+.001
    h1 = 0.01*(b1-a1)
    h2 = 0.01*(b2-a2)

    for(i in 1:100)
    {
      for(j in 1:100){
       p1 = a1+i*h1
       p = p1-h1
       q1 = a2+j*h2
       q = q1-h2
       c0 = Cop(c(p,q))
       c1 = Cop(c(p1,q))
       c2 = Cop(c(p,q1))
       M = cbind(  c1-c0,c2-c0    )
       R=Matrix::rankMatrix(M)
       if(R[1]<2){break}
      }
    }

    out = (R[1]==2)
  }
  else if(family %in% cpar3)
  {
    X1 = data[,1]
    X2 = data[,2]
    out1 = preparedata(X1)
    out2 = preparedata(X2)
    F1   = out1$Fn
    F2   = out2$Fn
    n1   = out1$m-1
    n2   = out2$m-1
    F1 = F1[1:n1]
    F2 = F2[1:n2]
    d = n1*n2

    u1 = rep(F1,n2)
    u2 = rep(F2,each=n1)
    u = cbind(u1,u2)
    Cop = function(par) {pncs(u,family,rotation,par)}
    out0 = est_options(family)
    UB = out0$UB
    LB = out0$LB

    b1 = UB[1]-.001
    b2 = UB[2]-.005
    b3 = UB[3]-.005
    a1 = LB[1]+.001
    a2 = LB[2]+.005
    a3 = LB[3]+.005
    h1 = 0.01*(b1-a1)
    h2 = 0.1*(b2-a2)
    h3 = 0.1*(b3-a3)

    for(i in 1:100)
    {
      for(j in 1:10){
        for(k in 1:10)
        {
        p1 = a1+i*h1
        p = p1-h1
        q1 = a2+j*h2
        r1 = a3+k*h3
        q = q1-h2
        r = r1-h3
        c0 = Cop(c(p,q,r))
        c1 = Cop(c(p1,q,r))
        c2 = Cop(c(p,q1,r))
        c3 = Cop(c(p,q,r1))
        M = cbind(  c1-c0,c2-c0 ,c3-c0   )
        R=Matrix::rankMatrix(M)
        if(R[1]<3){break}
        }
      }
    }

    out = (R[1]==3)
  }
  out
}
