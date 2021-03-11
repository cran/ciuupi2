optimize_knots <- function(lambda, rho, alpha, t.alpha, gams, d, m,
                           n.ints, knots, knots.all, nodes, weights,
                           wvec, psinu.zvec, h, cons, natural, obj,
                           start.vec){
  # Find the value of the b s vector
  # that specifies the Kabaila and Giri confidence interval,
  # for a given value of lambda.
  # This vector is found by numerical constrained
  # optimization.
  #
  # Inputs:
  # lambda: a positive tuning parameter
  # rho: a known correlation
  # t.alpha: quantile of the t distribution for m and alpha
  # gams: set of gammas at which the coverage is
  #       required to be greater than or equal to 1 - alpha
  # d: the b and s functions are optimized in the interval (0, d]
  # m: degrees of freedom n - p
  # n.ints: number of intervals in [0,d]
  # knots: location of knots in [0, d]
  # knots.all: location of knots in [-d, d]
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # zvec: a vector of length N where outer integrand is evaluated at
  # wvec: g(zvec)/sqrt(m / (m + 1)) where g(z)=exp(z/2 - exp(-z))
  # psinu.zvec: f_m+1(g(z))*d(g(z))/dz evaluated at z=zvec
  # h: step length
  # cons: sqrt(2/m) * exp(lgamma((m+1)/2) - lgamma(m/2)) where
  #       m is the degrees of freedom
  # natural: equals to 1 for natural cubic spline interpolation
  #          or 0 for clamped cubic spline interpolation
  # obj: 1 for definition 1 of SEL or 2 for  definition 2 of SEL
  # start.vec: a starting vector to the optimization problem
  #
  # Output:
  # The b s vector.
  #
  # Written by N Ranathunga in September 2020


  # Specify lower and upper bounds on the vector of values
  # of the b and s functions evaluated at the knots
  low <- c(rep(-100, n.ints - 1), rep(0.5, n.ints))
  up <- c(rep(100, n.ints - 1), rep(200, n.ints))


  # Make the objective function a function of one argument, y
  if (obj == 1) {
  obj_fun <- functional::Curry(objective1, lambda = lambda, m = m, n.ints = n.ints,
                                 knots = knots, knots.all = knots.all, t.alpha = t.alpha,
                                 nodes = nodes, weights = weights, natural = natural)
  } else {
  obj_fun <- functional::Curry(objective2, lambda = lambda, m = m, n.ints = n.ints,
                                 knots = knots, knots.all = knots.all, t.alpha = t.alpha,
                                 nodes = nodes, weights = weights, natural = natural)
  }


  # Make the constraint function a function of one argument, y
  cons_fun <- functional::Curry(constraints_slsqp_trapez, gams = gams, rho = rho,
                                n.ints = n.ints, knots = knots, knots.all = knots.all,
                                alpha = alpha, t.alpha = t.alpha, nodes = nodes,
                                weights = weights, wvec = wvec, psinu.zvec = psinu.zvec,
                                h = h, cons = cons, natural = natural)

  # Find the values of the knots using the optimization function
  res <- nloptr::slsqp(start.vec, obj_fun, hin = cons_fun, lower = low,
               upper = up, nl.info = FALSE)
  new.par <- res$par

  # Output the vector with knot values which specifies the new
  # confidence interval
  out <- new.par

}
