spline_s <- function(y, n.ints, knots.all, t.alpha, natural){
  # Return the value of the s function at a given point x.
  #
  # Inputs:
  # y: contains knot values of the b and s functions
  # n.ints: number of intervals in (0, d]
  # knots.all: location of knots in [-d, d]
  # t.alpha: quantile of the t distribution for m and alpha
  # natural: equals to 1 for natural cubic spline interpolation
  #          or 0 for clamped cubic spline interpolation
  #
  # Written by R Mainzer, March 2017
  # Modified by N Ranathunga in September 2020

  s.vals <- c(y[n.ints:(2 * n.ints - 1)], t.alpha)
  s.vals.all <- c(rev(s.vals), s.vals[2:(n.ints+1)])

  if(natural == 1){
    s.spl <- stats::splinefun(knots.all, s.vals.all, method = "natural")
  } else {
    s.spl.pp <- pracma::cubicspline(knots.all, s.vals.all, endp2nd = TRUE)
    s.spl <- function(x) pracma::ppval(s.spl.pp, x)
  }

  out <- s.spl

}

