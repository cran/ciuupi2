compute_lambda <- function(n.iter, rho, alpha, t.alpha, gams, d, m,
                           n.ints, N, eps, knots, knots.all, nodes,
                           weights, natural, obj, start.vec){
  # This program will find the optimized value of lambda for
  # given rho and alpha.
  #
  # Details for finding the optimized value of lambda:
  # For 5 iterations, the bisection method is used to find the
  # value of lambda such that the SEL "loss" is equal to the SEL
  # "gain".  Once these iterations are performed, a final answer
  # for lambda is found by fitting a straight line to the two
  # last iteration values, and finding the x-axis intercept of
  # this straight line.
  #
  # Inputs:
  # n.iter: number of iterations
  # rho: known correlation
  # alpha: 1 - alpha is the minimum coverage probability of the
  #        confidence interval
  # t.alpha: quantile of the t distribution for m and alpha
  # gams: constrain coverage probability at these values
  # d: the b and s functions are optimized in the interval (0, d]
  # m: degrees of freedom n - p
  # n.ints: number of intervals in [0,d]
  # N: number of evaluations in the outer integrand
  # eps: upper bound of the approximation error
  # knots: location of knots in [0, d]
  # knots.all: location of knots in [-d, d]
  # nodes: vector of Gauss Legendre quadrature nodes
  # weights: vector of Gauss Legendre quadrature weights
  # natural: 1 (default) for natural cubic spline interpolation
  #           or 0 for clamped cubic spline interpolation
  # obj: 1 for definition 1 of SEL or 2 for  definition 2 of SEL
  # start.vec: a starting vector to the optimization problem
  #
  # Output:
  # The optimized value of lambda
  #
  # Written by R Mainzer, August 2017
  # Modified by N Ranathunga in September 2020


  # The lower and upper bounds of the initial search interval
  # for the bisection root finding method
  lower <- 0
  upper <- 0.3

  # Set up vectors to store results
  res.lambda <- rep(0, n.iter)
  res.fun <- rep(0, n.iter)

  # Implement the bisection method
  for(i in 1:n.iter){

    tmp.lambda <- (upper + lower)/2
    tmp.fun <- compute_ratio_minus(tmp.lambda, rho, alpha, t.alpha,
                                    gams, d, m, n.ints, N, eps, knots,
                                    knots.all, nodes, weights, natural,
                                    obj, start.vec)

      if(tmp.fun < 0){
      lower <- tmp.lambda
    } else{
      upper <- tmp.lambda
    }

    res.lambda[i] <- tmp.lambda
    res.fun[i] <- tmp.fun

  }

  # Find lambda by a linear interpolation of the last
  # upper and lower values
  x1 <- lower
  if(tmp.lambda == x1){
    y1 <- tmp.fun
  } else {
    y1 <- compute_ratio_minus(x1, rho, alpha, t.alpha,
                               gams, d, m, n.ints, N, eps, knots,
                               knots.all, nodes, weights, natural,
                               obj, start.vec)
  }

  x2 <- upper
  if(tmp.lambda == x2){
    y2 <- tmp.fun
  } else {
    y2 <- compute_ratio_minus(x2, rho, alpha, t.alpha,
                               gams, d, m, n.ints, N, eps, knots,
                               knots.all, nodes, weights, natural,
                               obj, start.vec)
  }

  slope <- (y2 - y1) / (x2 - x1)
  x3 <- x1 - (y1 / slope)

  # Do the linear interpolation one more time
  y3 <- compute_ratio_minus(x3, rho, alpha, t.alpha,
                             gams, d, m, n.ints, N, eps, knots,
                             knots.all, nodes, weights, natural,
                             obj, start.vec)

  if (sign(y1) != sign(y3)){
    slope <- (y3 - y1) / (x3 - x1)
    lambda <- x1 - (y1 / slope)
  } else {
    slope <- (y3 - y2) / (x3 - x2)
    lambda <- x2 - (y2 / slope)
  }

  out <- lambda

}
