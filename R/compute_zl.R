compute_zl <- function(nu, eps){
  # Compute the length of the interval to apply
  # the trapezoidal rule to the outer integral and
  # compute the location of the first evaluation
  # of the outer integrand.
  #
  # Input:
  # nu: degrees of freedom + positive integer
  # eps: upper bound of the error when approximate the
  #      outer integral by the trapezoidal rule
  #
  # Output:
  # A list of two values
  #
  # Written by N Ranathunga, September 2020

  # Find the length of the interval
  d2 <- stats::uniroot(MinUpBndTrErrMinEpsFin, nu, eps,
                interval = c(0, 10), extendInt="yes")$root

  # Find the location of the first evaluation
  zl <- stats::optimize(UpBoundTrimError, d2, nu, interval=c(-5, 1))$minimum

  out <- list(d2, zl)

}

