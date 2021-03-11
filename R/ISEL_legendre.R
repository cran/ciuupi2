ISEL_legendre <- function(gam, w, knots, t.alpha,
                          nodes, weights, s.spl){
  # Computes the inner integral of the scaled expected
  # length 1 and 2 of Kabaila and Giri confidence interval.
  # The integral from (0, d) is broken down to integrals
  # over knots. Each integral is computed using gauss
  # legendre quadrature. The number of nodes and weights
  # for the approximation of each integral can be changed.
  #
  # Inputs:
  # gam: parameter
  # w: a value of the variable of integration in the 
  #    outer integral
  # knots: location of knots in [0, d]
  # t.alpha: quantile of the t distribution for m and alpha
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # s.spl: s function
  #
  # Output:
  # A value for the inner integral of the scaled expected
  # length 1 and 2
  #
  # Written by N.Ranathunga in September 2020.

  # Set up a vector to store the results
  int <- rep(0, length(knots))

  for(i in 1:(length(knots) - 1)){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    q <- IISEL(adj.nodes, gam, w, t.alpha, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  ISEL <- sum(int)

}
