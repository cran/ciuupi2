IISEL <- function(x, gam, w, t.alpha, s.spl){
  # Evaluates the function
  # (s(x) - t_alpha) * (phi(wx-gamma) + phi(wx+gamma))
  # for a vector x.
  #
  # Inputs:
  # x: vector of nodes of the Gauss Legendre quadrature
  # gam: parameter
  # m: degrees of freedom n - p
  # w: a value of the variable of integration in the 
  #    outer integral
  # t.alpha: quantile of the t distribution for m and alpha
  # s.spl: s function
  #
  # Output:
  # A vector with the same dimension as x.
  #
  # Written by N.Ranathunga in September 2020.

  tmp1 <- s.spl(x) - t.alpha
  tmp2 <- stats::dnorm(w*x - gam, 0, 1) +
    stats::dnorm(w*x + gam, 0, 1)
  res <- tmp1 * tmp2

}
