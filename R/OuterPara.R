OuterPara <- function(m, nu, N, eps){
  # Compute a list of values and vectors which will 
  # be used as inputs to compute the outer integrals 
  # of the coverage probability and the two squared scaled 
  # expected lengths.
  #
  # Inputs:
  # m: degrees of freedom n - p
  # nu: m + positive integer
  # N: number of evaluations in the outer integrand
  # eps: upper bound of the approximation error
  #
  # Output:
  # A list of values and vectors.
  #
  # Written by N Ranathunga in September 2020
  
  # Step length to compute the outer integral 
  cm <- sqrt(m / nu)
  d2 <- as.numeric(compute_zl(nu, eps)[1])
  zl <- as.numeric(compute_zl(nu, eps)[2])
  zu <- zl + d2
  h <- d2/(N - 1)
  
  # vectors which will be used to compute
  # the outer integral
  zvec <- seq(zl, zu, by = h)
  wvec <- transf(zvec) / cm
  psinu.zvec <- Func_psi_nu(zvec, nu)
  
  # The constant term in CP for nu=m+1 or
  # the constant term in SEL2 for nu=m+1
  cons1 <- sqrt(2/m) * exp(lgamma(nu/2) - lgamma(m/2))
  
  # The constant term in SEL1 for nu=m+2
  cons2 <- (2/m) * exp(lgamma(nu/2) - lgamma(m/2))
  
  # E(W) term in SEL1
  exp.w <- sqrt(2/m) * exp(lgamma((m + 1)/2) - lgamma(m/2))
  
  out <- list(h=h, cons1=cons1, cons2=cons2, exp.w=exp.w,
              wvec=wvec, psinu.zvec=psinu.zvec)

}


