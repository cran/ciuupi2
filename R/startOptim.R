startOptim <- function(n.ints, t.alpha){
  # Calculate a vector which provides
  # starting values to optimize the objective
  # function.
  #
  # Inputs:
  # n.ints: number of intervals in (0, d]
  # t.alpha: quantile of the t distribution for 
  #          m and alpha
  #
  # Output: 
  # A vecror of length 2*n.ints - 1.
  #
  # Written by N Ranathunga in September 2020

  out <- c(rep(0, n.ints - 1), rep(t.alpha, n.ints))

}
