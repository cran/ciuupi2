IICP <- function(x, gam, rho, w, t.alpha, b.spl, s.spl){
  # Compute the function
  # (k(x, w, gam, rho) - k_dag(x, w, gam, rho)) *
  # phi(wx - gam) +
  # (k(-x, w, gam, rho) - k_dag(-x, w, gam, rho)) *
  # phi(wx + gam)
  # for a vector x.
  #
  # Inputs:
  # x: vector of nodes of the Gauss Legendre quadrature
  # gam: parameter
  # rho: a known correlation
  # w: a value of the variable of integration in the 
  #    outer integral
  # t.alpha: quantile of the t distribution for m and alpha
  # b.spl: b function
  # s.spl: s function
  #
  # Output:
  # A vector with the same dimension as x.
  #
  # Written by N. Ranathunga in September 2020


  # Finding k_dag(x, w, gam, rho))
  mu1 <- rho * (w*x - gam)
  var <- 1 - rho^2
  k.dag1 <- Psi(-t.alpha * w, t.alpha * w, mu1, var)

  # Finding k(x, w, gam, rho)
  term.a1 <- b.spl(x)
  term.b1 <- s.spl(x)
  lh <- w * (term.a1 - term.b1)
  uh <- w * (term.a1 + term.b1)
  k1 <- Psi(lh, uh, mu1, var)

  # Finding phi(wx - gam)
  term1 <- stats::dnorm(w*x - gam, 0, 1)

  # Finding k_dag(-x, w, gam, rho))
  mu2 <- rho * (-w*x - gam)
  k.dag2 <- Psi(-t.alpha * w, t.alpha * w, mu2, var)

  # Finding k(-x, w, gam, rho)
  term.a2 <- b.spl(-x)
  term.b2 <- s.spl(-x)
  lh2 <- w * (term.a2 - term.b2)
  uh2 <- w * (term.a2 + term.b2)
  k2 <- Psi(lh2, uh2, mu2, var)

  # Finding phi(wx + gam)
  term2 <- stats::dnorm(w*x + gam, 0, 1)

  res <- (k1 - k.dag1) * term1 + (k2 - k.dag2) * term2

}
