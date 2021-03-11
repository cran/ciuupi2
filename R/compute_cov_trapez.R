compute_cov_trapez <- function(gam, rho, knots, alpha, t.alpha, 
                               nodes, weights, b.spl, s.spl, wvec, 
                               psinu.zvec, h, cons){
  # Compute the coverage probability of the Kabaila and Giri
  # confidence interval.
  # Apply the transformation (2.6) of Mori (1988), followed  by
  # the trapezoidal rule to the outer integral.
  #
  # Inputs:
  # gam: parameter
  # rho: a known correlation
  # knots: location of knots in [0, d]
  # alpha: nominal coverage is 1 - alpha
  # t.alpha: quantile of the t distribution for m and alpha
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # b.spl: b function
  # s.spl: s function
  # wvec: g(zvec)/sqrt(m / (m + 1)) where g(z)=exp(z/2 - exp(-z))
  # psinu.zvec: f_m+1(g(z))*d(g(z))/dz evaluated at z=zvec
  # h: step length
  # cons: sqrt(2/m) * exp(lgamma((m+1)/2) - lgamma(m/2)) where
  #       m is the degrees of freedom
  #
  # Written by N Ranathunga, September 2020


  # Set up a vector to store the results of ICP
  ICP.zvec <- rep(0, length(wvec))

  for(i in 1:length(wvec)){
  w <- wvec[i]
  ICP.zvec[i] <- ICP_legendre(gam, rho, w, knots, t.alpha,
                              nodes, weights, b.spl, s.spl)
  }

  out.int <- h * PreciseSums::kahanSum(ICP.zvec * psinu.zvec)

  cp <- (1 - alpha) + cons * out.int

}
