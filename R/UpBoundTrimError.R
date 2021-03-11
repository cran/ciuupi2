UpBoundTrimError <- function(z, d2, nu){
  # Computes the upper bound of the
  # trimming error when approximate
  # the outer integral by the trapezoidal rule.
  #
  # Input:
  # z: a given value
  # d2: (number of evaluations of the outer
  #     integrand) * (step length)
  # nu: degrees of freedom + positive integer
  #
  # Written by N Ranathunga, September 2020

  x1 <- transf(z)
  x2 <- transf(z + d2)
  term1 <- nu * x1^2
  term2 <- nu * x2^2

  out <- stats::pchisq(term1, df=nu) + 1 -
         stats::pchisq(term2, df=nu)

}
