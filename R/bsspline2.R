#'Evaluate the functions b and s at x
#'
#'Evaluate the functions b and s, as specified by the vector
#'(b(d/6),b(2d/6),...,b(5d/6),s(0),s(d/6),...,s(5d/6)) computed using
#'\code{bsciuupi2}, \code{alpha}, \code{m}
#'and \code{natural} at \code{x}.
#'
#'@param x A value or vector of values at which the functions b and s are to be
#'  evaluated
#'@param bsvec The vector \eqn{(b(d/6), b(2d/6), \dots, b(5d/6), s(0), s(d/6),
#'  \dots, s(5d/6))} computed using \code{bsciuupi2}
#'@param alpha The minimum coverage probability is 1 - \code{alpha}
#'@param m Degrees of freedom \code{n - p}
#'@param natural Equal to 1 (default) if the b and s functions are evaluated by
#'  natural cubic spline interpolation or 0 if evaluated by clamped cubic spline
#'  interpolation. This parameter must take the same value as that used in
#'  \code{bsciuupi2}
#'
#'@return A data frame containing \code{x} and the corresponding values of the
#'  functions b and s.
#'
#'
#'@details The function b is an odd continuous function and the function s is an
#'  even continuous function. In addition, b(x)=0 and s(x) is equal to the
#'  \eqn{1 - \alpha/2} quantile of the \eqn{t} distribution with \code{m}
#'  degrees of freedom for all |x| greater than or equal to d, where d is a
#'  sufficiently large positive number (chosen by the function
#'  \code{bsciuupi2}). The values of these functions in the interval
#'  \eqn{[-d,d]} are specified by the vector \eqn{(b(d/6), b(2d/6), \dots,
#'  b(5d/6), s(0), s(d/6), \dots, s(5d/6))} as follows. By assumption,
#'  \eqn{b(0)=0} and \eqn{b(-i)=-b(i)} and \eqn{s(-i)=s(i)} for
#'  \eqn{i=d/6,...,d}. The values of \eqn{b(x)} and \eqn{s(x)} for any \eqn{x}
#'  in the interval \eqn{[-d,d]} are found using cubic spline interpolation for
#'  the given values of \eqn{b(i)} and \eqn{s(i)} for
#'  \eqn{i=-d,-5d/6,...,0,d/6,...,5d/6,d}. The choices of \eqn{d} for \eqn{m =
#'  1, 2} and \eqn{>2} are \eqn{d=20, 10} and \eqn{6} respectively.
#'
#'
#'  The vector \eqn{(b(d/6), b(2d/6), \dots, b(5d/6), s(0), s(d/6), \dots,
#'  s(5d/6))} that specifies the Kabaila and Giri(2009) confidence interval that
#'  utilizes uncertain prior information (CIUUPI), with minimum coverage
#'  probability \code{1 - alpha}, is obtained using
#'  \code{\link{bsciuupi2}}.
#'
#'  In the examples, we continue with the same 2 x 2 factorial example described
#'  in the documentation for \code{\link{find_rho}}.
#'
#'@seealso \code{\link{find_rho}}, \code{\link{bsciuupi2}}
#'
#'@references Kabaila, P. and Giri, R. (2009).  Confidence intervals in
#'  regression utilizing prior information. Journal of Statistical Planning and
#'  Inference, 139, 3419-3429.
#'
#' @examples
#' alpha <- 0.05
#' m <- 8
#'
#' # Find the vector (b(d/6),...,b(5d/6),s(0),...,s(5d/6)) that specifies the
#' # Kabaila & Giri (2009) CIUUPI for the first definition of the
#' # scaled expected length (default) (takes about 30 mins to run):
#' \donttest{
#' bsvec <- bsciuupi2(alpha, m, rho = -0.7071068)
#' }
#'
#' # The result bsvec is (to 7 decimal places) the following:
#' bsvec <- c(-0.0287487, -0.2151595, -0.3430403, -0.3125889, -0.0852146,
#'             1.9795390,  2.0665414,  2.3984471,  2.6460159,  2.6170066,  2.3925494)
#'
#' # Graph the functions b and s
#' x <- seq(0, 8, by = 0.1)
#' splineval <- bsspline2(x, bsvec, alpha, m)
#'
#' plot(x, splineval[, 2], type = "l", main = "b function",
#' ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue")
#' plot(x, splineval[, 3], type = "l", main = "s function",
#' ylab = " ", las = 1, lwd = 2, xaxs = "i",  col = "blue")
#'
#'@export

bsspline2 <- function(x, bsvec, alpha, m, natural = 1){

  # Set input
  n.ints <- 6
  d <- choice_d(m)$d

  # Specify where the knots for b and s are located
  knots.all <- seq(-d, d, by = d/n.ints)

  # Set t.alpha for m and alpha
  t.alpha <- stats::qt(1 - alpha/2, m)

  # Find b and s functions
  sspl <- spline_s(bsvec, n.ints, knots.all, t.alpha, natural)
  bspl <- spline_b(bsvec, n.ints, knots.all, t.alpha, natural)

  x1 <- x[which(x <= -d)]
  x2 <- x[which(x > -d & x < d)]
  x3 <- x[which(x >= d)]

  bspl.res <- c(rep(0, length(x1)), bspl(x2), rep(0, length(x3)))
  sspl.res <- c(rep(t.alpha, length(x1)), sspl(x2), rep(t.alpha, length(x3)))

  out <- data.frame(x = x, b = bspl.res, s = sspl.res)

}
