#'@title Family number corresponding to VineCopula package
#'
#'@description Computes the  number associated with a copula family  (without rotation)
#'@param family  Copula family: "gaussian", "t", "clayton", "frank", "gumbel", "joe", "plackett'', "bb1", "bb6", "bb7","bb8".
#'
#'
#'@return \item{fnumber}{Number}

#'@references Nagler et al. (2023). VineCopula: Statistical Inference of Vine Copulas, version 2.4.5.
#'
#'@examples
#'fnumber("bb1")
#'
#'@export
#'
#'

fnumber = function(family)
{
switch(family,
       "gaussian"= {fnum=1},

       "t" = { fnum=2},

       "clayton" = { fnum=3},

       "gumbel" = {fnum=4},


       "frank" = {fnum=5},

       "joe" = { fnum=6},

       "bb1" = { fnum=7},

       "bb6" = { fnum=8},

       "bb7" = { fnum=9},

       "bb8" = { fnum=10}


      )

  fnum
}


