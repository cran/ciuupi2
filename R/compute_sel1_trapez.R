compute_sel1_trapez <- function(gam, knots, t.alpha, nodes, 
                                weights, s.spl, wvec, psinu.zvec, 
                                h1, cons, exp.w){
  # Compute the value of the scaled expected length 1 for 
  # given functions b and s.
  # In other words, this function computes
  #
  # 1 + (1/(t_alpha E(W))) int_0^infty (ISEL1(w,gam)) w^2 fW(w) dw
  #
  # Inputs:
  # gam: parameter
  # knots: location of knots in [0, d]
  # t.alpha: quantile of the t distribution for m and alpha
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # s.spl: s function
  # wvec: g(z)/sqrt(m / (m + 2)) where g(z)=exp(z/2 - exp(-z))
  #       evaluated at z=zvec
  # psinu.zvec: f_m+2(g(z))*d(g(z))/dz evaluated at z=zvec
  # h1: step length
  # cons: (2/m) * exp(lgamma((m+2)/2) - lgamma(m/2)) where
  #       m is the degrees of freedom
  # exp.w: sqrt(2/m) * exp(lgamma((m + 1)/2) - lgamma(m/2))
  #
  # Output:
  # The scaled expected length 1 for given functions b and s.
  #
  # Written by N Ranathunga, September 2020


  # Set up a vector to store the results of ISEL1
  ISEL1.zvec <- rep(0, length(wvec))

  for(i in 1:length(wvec)){
    w <- wvec[i]
    ISEL1.zvec[i] <- ISEL_legendre(gam, w, knots, t.alpha,
                                   nodes, weights, s.spl)
  }

  out.int <- h1 * PreciseSums::kahanSum(ISEL1.zvec * psinu.zvec)

  sel1 <- 1 + (cons * out.int) / (t.alpha * exp.w)

}
