compute_sel2_trapez <- function(gam, knots, t.alpha, nodes, weights,
                                s.spl, wvec, psinu.zvec, h2, cons){
  # Compute the value of the scaled expected
  # length 2 for given functions b and s.
  # In other words, this function computes
  #
  # 1 + (1/t_alpha) int_0^infty (ISEL(w,gam)) w fm(w) dw
  #
  # Inputs:
  # gam: parameter
  # knots: location of knots in [0, d]
  # t.alpha: quantile of the t distribution for m and alpha
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # s.spl: s function
  # wvec: g(z)/sqrt(m / (m + 1)) where g(z)=exp(z/2 - exp(-z))
  #       evaluated at z=zvec
  # psinu.zvec: f_m+1(g(z))*d(g(z))/dz evaluated at z=zvec
  # h2: step length
  # cons: sqrt(2/m) * exp(lgamma((m+1)/2) - lgamma(m/2)) where
  #       m is the degrees of freedom
  #
  # Output:
  # The scaled expected length 2 for given functions b and s.
  #
  # Written by N Ranathunga, September 2020


  # Set up a vector to store the results of ISEL2
  ISEL2.zvec <- rep(0, length(wvec))

  for(i in 1:length(wvec)){
    w <- wvec[i]
    ISEL2.zvec[i] <- ISEL_legendre(gam, w, knots, t.alpha,
                                   nodes, weights, s.spl)
  }

  out.int <- h2 * PreciseSums::kahanSum(ISEL2.zvec * psinu.zvec)

  sel2 <- 1 + (cons * out.int) / t.alpha

}
