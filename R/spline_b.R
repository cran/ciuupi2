spline_b <- function(y, n.ints, knots.all, t.alpha, natural){
  # Return the value of the b function at a given point x.
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

  y.rev <- rev(y[1:(n.ints - 1)])
  b.vals <- c(0, y[1:(n.ints - 1)], 0)
  b.vals.all <- c(0, -y.rev, b.vals)

  # If natural = 1 use natural cubic spline, otherwise use clamped cubic
  # spline
  if(natural == 1){
    b.spl <- stats::splinefun(knots.all, b.vals.all, method = "natural")
  } else {
    b.spl.pp <- pracma::cubicspline(knots.all, b.vals.all, endp2nd = TRUE)
    b.spl <- function(x) pracma::ppval(b.spl.pp, x)
  }

  out <- b.spl

}

