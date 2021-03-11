compute_ratio_minus <- function(lambda, rho, alpha, t.alpha, gams,
                                d, m, n.ints, N, eps, knots, knots.all,
                                nodes, weights, natural, obj, start.vec){
  # Compute the ratio ( gain / maximum possible loss) - 1.
  # Another program can then use this program to find the
  # value of lambda which makes (this ratio - 1) = 0.
  #
  # Inputs:
  # lambda: used in specifying the objective function
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
  # The ratio (gain / maximum possible loss) - 1.
  #
  # Written by P.Kabaila in June 2008
  # Rewritten in R by R Mainzer, March 2017
  # Modified by N Ranathunga in September 2020


  # Specify the values of the inputs to compute the
  # outer integral in coverage probability and SEL2
  h <- OuterPara(m, nu=m+1, N, eps)$h
  wvec <- OuterPara(m, nu=m+1, N, eps)$wvec
  psinu.zvec <- OuterPara(m, nu=m+1, N, eps)$psinu.zvec
  cons <- OuterPara(m, nu=m+1, N, eps)$cons1

  new.par <- optimize_knots(lambda, rho, alpha, t.alpha, gams, d, m,
               n.ints, knots, knots.all, nodes, weights,  wvec,
               psinu.zvec, h, cons, natural, obj, start.vec)

  s.spl <- spline_s(new.par, n.ints, knots.all, t.alpha, natural)

  # Compute the required ratio
  if (obj == 1){

    h1 <- OuterPara(m, nu=m+2, N, eps)$h
    wvec1 <- OuterPara(m, nu=m+2, N, eps)$wvec
    psinu.zvec1 <- OuterPara(m, nu=m+2, N, eps)$psinu.zvec
    cons1 <- OuterPara(m, nu=m+2, N, eps)$cons2
    exp.w1 <- OuterPara(m, nu=m+2, N, eps)$exp.w

    sel.max <- stats::optimize(compute_sel1_trapez, c(0, d), maximum = TRUE,
                  knots=knots, t.alpha=t.alpha, nodes=nodes, weights=weights,
                  s.spl=s.spl, wvec=wvec1, psinu.zvec=psinu.zvec1,
                  h1=h1, cons=cons1, exp.w=exp.w1)$objective

    sel.min <- compute_sel1_trapez(gam = 0, knots, t.alpha, nodes, weights,
                   s.spl, wvec=wvec1, psinu.zvec=psinu.zvec1,
                   h1=h1, cons=cons1, exp.w=exp.w1)

  } else if (obj == 2) {

    sel.max <- stats::optimize(compute_sel2_trapez, c(0, d), maximum = TRUE,
                   knots=knots, t.alpha=t.alpha, nodes=nodes, weights=weights,
                   s.spl=s.spl, wvec=wvec, psinu.zvec=psinu.zvec,
                   h2=h, cons=cons)$objective

    sel.min <- compute_sel2_trapez(gam = 0, knots, t.alpha, nodes, weights,
                   s.spl=s.spl, wvec, psinu.zvec,
                   h2=h, cons)
  }

  expected.gain <- 1 - sel.min^2
  max.potential.loss <- sel.max^2 - 1

  # Output the required ratio minus 1
  out <- expected.gain / max.potential.loss - 1

}
