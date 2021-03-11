integrand_obj1 <- function(x, lambda, m, t.alpha, s.spl){
  # Evaluate the inner integrand of the objective function
  # of the scaled expected length 1 for a vector x.
  # In other words, this function computes
  #
  # (s(x) - t_alpha) (lambda +
  # [((m/(x^2 + m))^((m/2) + 1)  / sqrt(2 pi))] )
  #
  # Inputs:
  # x: vector at which the integrand is to be evaluated
  # lambda: a positive tuning parameter
  # m: degrees of freedom n - p
  # t.alpha: quantile of the t distribution for m and alpha
  # s.spl: s function
  #
  # Output:
  # A vector of values of the inner integrand with the
  # same dimension as x.
  #
  # Written by N Ranathunga in September 2020

  tmp1 <- s.spl(x) - t.alpha
  term1 <- (1/sqrt(2 * pi)) * (m / (x^2 + m)) ^ (m/2 + 1)
  tmp2 <- lambda + term1

  res <- tmp1 * tmp2

}
