constraints_slsqp_trapez <- function(gams, rho, y, n.ints, knots, knots.all, 
                                     alpha, t.alpha, nodes, weights, wvec, 
                                     psinu.zvec, h, cons, natural){
  # This function computes (coverage probability) - (1 - alpha)
  # for a vector of gamma values.
  #
  # Inputs:
  # gams: set of gammas at which the coverage is
  #       required to be greater than or equal to 1 - alpha
  # rho: a known correlation
  # y: contains knots values of the b and s functions
  # n.ints: number of intervals in (0, d]
  # knots: location of knots in [0, d]
  # knots.all: location of knots in [-d, d]
  # alpha: nominal coverage is 1 - alpha
  # t.alpha: quantile of the t distribution for m and alpha
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # wvec: g(zvec)/sqrt(m / (m + 1)) where g(z)=exp(z/2 - exp(-z))
  # psinu.zvec: f_m+1(g(z))*d(g(z))/dz evaluated at z=zvec
  # h: step length
  # cons: sqrt(2/m) * exp(lgamma((m+1)/2) - lgamma(m/2)) where
  #       m is the degrees of freedom
  # natural: equals to 1 for natural cubic spline interpolation
  #          or 0 for clamped cubic spline interpolation
  #
  # Output:
  # A vector of values of (coverage probability) - (1 - alpha)
  #
  # Written by N Ranathunga in September 2020

  len.gams <- length(gams)
  covs <- rep(0, len.gams)

  # Find b and s functions at y
  b.spl <- spline_b(y, n.ints, knots.all, t.alpha, natural)
  s.spl <- spline_s(y, n.ints, knots.all, t.alpha, natural)

  for(i in 1:len.gams){
    covs[i] <- compute_cov_trapez(gams[i], rho, knots, alpha, t.alpha,
                                  nodes, weights, b.spl, s.spl, wvec,
                                  psinu.zvec, h, cons)
  }

  out <- covs - (1 - alpha)

}
