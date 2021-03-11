#'Compute the first definition of the scaled expected length of the Kabaila &
#'Giri (2009) CIUUPI
#'
#'Evaluate the first definition of the scaled expected length of the Kabaila &
#'Giri (2009) confidence interval that utilizes uncertain prior information
#'(CIUUPI), with minimum coverage 1 - \code{alpha}, at \code{gam}.
#'
#'@param gam A value of gamma or vector of gamma values at which the first
#'  definition of the scaled expected length function is evaluated
#'@param bsvec The vector (b(d/6),b(2d/6),...,b(5d/6),s(0),s(d/6),...,s(5d/6))
#'  computed using \code{bsciuupi2}
#'@param alpha The minimum coverage probability is 1 - \code{alpha}
#'@param m Degrees of freedom \code{n - p}
#'@param rho A known correlation
#'@param natural Equal to 1 (default) if the b and s functions are obtained by
#'  natural cubic spline interpolation or 0 if obtained by clamped cubic spline
#'  interpolation. This parameter must take the same value as that used in
#'  \code{bsciuupi2}
#'
#'@return The value(s) of the first definition of the scaled expected length of
#'  the Kabaila & Giri (2009) CIUUPI at \code{gam}.
#'
#'@details
#'
#'Suppose that \deqn{y = X \beta + \epsilon} where \eqn{y} is a random
#'\eqn{n}-vector of responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix with
#'linearly independent columns, \eqn{\beta} is an unknown parameter
#'\eqn{p}-vector and \eqn{\epsilon} is a random \eqn{n}-vector with components
#'that are independent and identically normally distributed with zero mean and
#'unknown variance. The parameter of interest is \eqn{\theta = } \code{a}'
#'\eqn{\beta}. The uncertain prior information is that \eqn{\tau = } \code{c}'
#'\eqn{\beta} takes the value \code{t}, where \code{a} and \code{c} are
#'specified linearly independent vectors and \code{t} is a specified number.
#'\code{rho} is the known correlation between the least squares estimators of
#'\eqn{\theta} and \eqn{\tau}. It is determined by the \eqn{n} by \eqn{p} design
#'matrix X and the \eqn{p}-vectors a and c using \code{\link{find_rho}}.
#'
#'The Kabaila & Giri (2009) CIUUPI is specified by the vector
#'(b(d/6),...,b(5d/6),s(0),...,s(5d/6)), \code{alpha}, \code{m} and
#'\code{natural}
#'
#'The first definition of the scaled expected length of the Kabaila and
#'Giri(2009) CIUUPI is the expected length of this confidence interval divided
#'by the expected length of the usual confidence interval with coverage
#'probability \code{1 - alpha}.
#'
#'In the examples, we continue with the same 2 x 2 factorial example described
#'in the documentation for \code{\link{find_rho}}.
#'
#'@seealso \code{\link{find_rho}}, \code{\link{bsciuupi2}}
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
#'
#' # Graph the squared scaled expected length function
#' gam <- seq(0, 10, by = 0.1)
#' sel <- sel1ciuupi2(gam, bsvec, alpha, m, rho = -0.7071068)
#' plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#' main = "Squared Scaled Expected Length", col = "blue",
#' xlab = expression(paste("|", gamma, "|")))
#' abline(h = 1, lty = 2)
#'
#'@references
#'
#'Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
#'prior information.  Journal of Statistical Planning and Inference, 139, 3419 -
#'3429.
#'
#'@export


sel1ciuupi2 <- function(gam, bsvec, alpha, m, rho, natural = 1){

  # Specify the values of the inputs to other functions
  n.ints <- 6
  N <- 33
  eps <- 10^{-10}
  n.nodes <- 20
  d <- choice_d(m)$d

  # Set t.alpha for m and alpha
  t.alpha <- stats::qt(1 - alpha/2, m)

  # Find the nodes and weights of the legendre quadrature
  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
  nodes <- quad.info$nodes
  weights <- quad.info$weights

  # Specify where the knots for b and s are located
  # as inputs to other functions
  knots <- seq(0, d, by = d/n.ints)
  knots.all <- seq(-d, d, by = d/n.ints)

  # Specify the values of the inputs to compute the
  # outer integral in SEL1
  h1 <- OuterPara(m, nu=m+2, N, eps)$h
  wvec <- OuterPara(m, nu=m+2, N, eps)$wvec
  psinu.zvec <- OuterPara(m, nu=m+2, N, eps)$psinu.zvec
  cons <- OuterPara(m, nu=m+2, N, eps)$cons2
  exp.w <- OuterPara(m, nu=m+2, N, eps)$exp.w

  # Specify the function s
  s.spl <- spline_s(bsvec, n.ints, knots.all, t.alpha, natural)

  # Compute the scaled expected length
  res <- rep(0, length(gam))
  for(i in 1:length(gam)){
    res[i] <- compute_sel1_trapez(gam[i], knots, t.alpha, nodes, weights,
                                  s.spl, wvec, psinu.zvec, h1, cons, exp.w)
  }

  out <- res

}
