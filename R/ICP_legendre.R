ICP_legendre <- function(gam, rho, w, knots, t.alpha,
                         nodes, weights, b.spl, s.spl){
  # Compute the inner integral of the coverage probability
  # of Kabaila and Giri confidence interval.
  # The integral from (0, d) is broken down to integrals 
  # over knots. Each integral is computed using gauss 
  # legendre quadrature. The number of nodes and weights for 
  # the approximation of each integral can be changed.
  #
  # Input:
  # gam: parameter
  # rho: a known correlation
  # w: a value of the variable of integration in the 
  #    outer integral
  # knots: location of knots in [0, d]
  # t.alpha: quantile of the t distribution for m and alpha
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # b.spl: b function
  # s.spl: s function
  #
  # Output:
  # A value for the inner integral of Kabaila and Giri
  # confidence interval.
  #
  # Written by N. Ranathunga in September 2020


  # Set up a vector to store the results
  int <- rep(0, length(knots))

  for(i in 1:(length(knots) - 1)){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    q <- IICP(adj.nodes, gam, rho, w, t.alpha, b.spl, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  ICP <- sum(int)

}
