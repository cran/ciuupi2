Func_psi_nu <- function(zvec, nu){
  # This function computes a vector of values
  # of f_nu(g(z)) * dg(z)/dz where g(z) = exp(z/2 - exp(-z))
  # and f_nu is the pdf of a random variable with the
  # same distribution as sqrt(Q/nu) where Q ~ chisq(nu)
  # and nu = m + positive integer
  #
  # Input:
  # zvec: a vector where function evaluations are at
  #       when applying the trapezoidal rule
  # nu: degrees of freedom + positive integer
  #
  # Output:
  # A vector of values of f_nu(g(z)) * dg(z)/dz 
  # with the same dimension as zvec.
  #
  # Written by N. Ranathunga in September 2020

  const <- exp( (nu/2) * log(nu) - lgamma(nu/2) -
                  ((nu/2) - 1) * log(2) )
  tmp1 <- exp(-zvec)
  term1 <- exp(nu*zvec/2 - nu*tmp1 - (nu/2)*exp(zvec - 2*tmp1))
  term2 <- 1/2 + tmp1

  psinu.val <- const * term1 * term2

}

