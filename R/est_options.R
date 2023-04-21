#'@title Options for the estimation of the parameters of bivariate copula-based models
#'
#'@description Sets starting values, upper and lower bounds for the parameters. The bounds are based on those in the rvinecopulib package.
#'@param  family Copula family: "gaussian", "t", "clayton", "frank", "gumbel", "joe", "plackett'', "bb1", "bb6", "bb7","bb8","ncs-gaussian", "ncs-clayton", "ncs-gumbel", "ncs-frank", "ncs-joe","ncs-plackett".
#'@param  tau   Estimated Kendall's tau to compute a starting point (default is 0.5)
#'
#'@return \item{LB}{Lower bound for the parameters}
#'@return \item{UB}{Upper bound for the parameters}
#'@return \item{start}{Starting point for the estimation}
#'
#'@references Nagler & Vatter (2002). rvinecopulib: High Performance Algorithms for Vine Copula Modeling. Version 0.6.2.1.3
#'@references Nasri (2020). On non-central squared copulas. Statistics and Probability Letters.
#'@references Nasri (2022). Test of serial dependence for arbitrary distributions. JMVA.
#'@references Nasri & Remillard (2023). Copula-based dependence measures for arbitrary data, arXiv 2301.07267.
#'
#'@examples
#'out = est_options("bb8")
#'
#'
#'@export
#'
#'

est_options = function(family,tau=0.5)
  {
  switch(family,
      "gaussian"= {LB = -1;
       UB = 1;
       start = sin(pi*tau/2)},

       "t" = { LB = c(-1,2);
       UB = c(1,50);
       start = c(sin(pi*tau/2),15)},

       "clayton" = {ttau = abs(tau);
       LB = 1e-10;
       UB = 28;
       start=2*ttau/(1-ttau);
       },

       "gumbel" = {ttau = abs(tau);
       LB = 1;
       UB = 50;
       start=1/(1-ttau)},

      "frank" = {LB = -35;
      UB = 35;
      start=5},

       "joe" = {LB = 1;
       UB = 30;
       start = 15},

       "bb1" ={LB = c(0,1);
       UB=c(7,7);
       start = c(3,5)},

       "bb6" ={ LB = c(1,1);
       UB = c(6,8);
       start = c(3,5)},

       "bb7"={ LB = c(1,0);
       UB = c(6,25);
       start = c(3,15)},

       "bb8"={ LB = c(1,1e-4);
       UB = c(8,1-1e-4);
       start = c(3,0.5)},

       "plackett" = {LB = 1e-4;
       UB = 1e4;
       start = 10},

      "ncs-gaussian"= {LB = c(-1,0,0);
      UB = c(1,3,3);
      start = c(sin(pi*tau/2),1,2)},


      "ncs-clayton" = {ttau = abs(tau);
      LB = c(1e-10,0,0);
      UB = c(28,3,3);
      start=c(2*ttau/(1-ttau),1,2);
      },

      "ncs-gumbel" = {ttau = abs(tau);
      LB = c(1,0,0);
      UB = c(50,3,3);
      start=c(1/(1-ttau),1,2)},

      "ncs-frank" = {LB = c(-35,0,0);
      UB = c(35,3,3);
      start=c(5,1,2)},

      "ncs-joe" = {LB = c(1,0,0);
      UB = c(30,3,3);
      start = c(15,1,2)},

      "ncs-plackett" = {LB = c(1e-4,0,0);
      UB = c(1e4,3,3);
      start = c(10,1,2)}
)
  out=list(LB=LB, UB=UB, start=start)
}
