objective2 <- function(y, lambda, m, n.ints, knots, knots.all,
                       t.alpha, nodes, weights, natural){
  # Computes the objective function of the scaled expected
  # length 2 of Kabaila and Giri confidence interval.
  # The integral from (0, d) is broken down to integrals
  # over knots. Each integral is computed using gauss
  # legendre quadrature. 
  #
  # Inputs:
  # y: contains knots values of the b and s functions
  # lambda: a positive tuning parameter
  # m: degrees of freedom n - p
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # t.alpha: quantile of the t distribution for m and alpha
  # knots: location of knots in [0, d]
  # knots.all: location of knots in [-d, d]
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # natural: equals to 1 for natural cubic spline interpolation
  #          or 0 for clamped cubic spline interpolation
  #
  # Output:
  # The value of the objective function.
  #
  # Written by N Ranathunga in September 2020

  s.spl <- spline_s(y, n.ints, knots.all, t.alpha, natural)

  # Set up a vector to store the results
  int <- rep(0, length(knots))

  for(i in 1:(length(knots) - 1)){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    q <- integrand_obj2(adj.nodes, lambda, m, t.alpha, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  out <- sum(int)

}
