MinUpBndTrErrMinEpsFin <- function(d2, nu, eps){
  # This function minimizes
  # (upper bound of the trimming error) - (10^(-3) * eps)
  # when approximated the outer integral
  # by the trapezoidal rule.
  #
  # Input:
  # d2: (number of evaluations of the outer
  #     integrand) * (step length)
  # nu: degrees of freedom + positive integer
  # eps: a given value for the upper bound of the
  #      approximation error
  #
  # Written by N Ranathunga, September 2020

  temp <- stats::optimize(UpBoundTrimError, d2, nu,
                   interval=c(-4, 4))

  out <- temp$objective - (10^(-3) * eps)

}
