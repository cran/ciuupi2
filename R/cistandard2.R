#'Compute the usual confidence interval
#'
#'Compute the usual 1 - \code{alpha} confidence interval
#'
#'@param X A known \eqn{n} by \eqn{p} matrix
#'@param a A \eqn{p}-vector used to specify the parameter of interest
#'@param y The \eqn{n}-vector of observed responses
#'@param alpha 1 - \code{alpha} is the coverage probability of the
#'  confidence interval
#'
#'@return The usual 1 - \code{alpha} confidence interval.
#'
#'@details Suppose that \deqn{Y = X \beta + \epsilon} is a random \eqn{n}-vector
#'  of responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix with linearly
#'  independent columns, \eqn{\beta} is an unknown parameter \eqn{p}-vector and
#'  \eqn{\epsilon} is the random error with components that
#'   are independent and identically normally distributed with zero mean and
#'   unknown variance. The parameter of interest is \eqn{\theta = } \code{a}'
#'   \eqn{\beta}, where \code{a} is a specified \eqn{p}-vector.
#'   Then \code{cistandard2}
#'  computes the usual 1 - \code{alpha} confidence interval for \eqn{\theta},
#'  for given \eqn{n}-vector of observed responses \code{y}.
#'
#'
#'  In the examples, we continue with the same 2 x 2 factorial example described
#'  in the documentation for \code{\link{find_rho}}, for the vector of observed
#'  responses \eqn{y} = (-1.3, 0.8, 2.6, 5.8, 0.3, 1.3, 4.3, 5.0, -0.4, 1.0,
#'  5.2, 6.2).
#'
#'   The design matrix \eqn{X} and the vector \code{a} (denoted in
#'   R by a.vec) are entered into R using the commands
#'   in the following example.
#'
#' @examples
#' col1 <- rep(1,4)
#' col2 <- c(-1, 1, -1, 1)
#' col3 <- c(-1, -1, 1, 1)
#' col4 <- c(1, -1, -1, 1)
#' X.single.rep <- cbind(col1, col2, col3, col4)
#' X <- rbind(X.single.rep, X.single.rep, X.single.rep)
#' a.vec <- c(0, 2, 0, -2)
#' y <- c(-1.3, 0.8, 2.6, 5.8, 0.3, 1.3, 4.3, 5.0, -0.4, 1.0, 5.2, 6.2)
#'
#' # Calculate the usual 95% confidence interval
#' res <- cistandard2(X, a=a.vec, y, alpha = 0.05)
#' res
#'
#' # The usual 1 - alpha confidence interval for theta is (-0.08185, 3.08185)
#'
#'@seealso \code{\link{find_rho}}
#'
#'@references Kabaila, P. and Giri, K. (2009) Confidence intervals in regression
#'  utilizing prior information.  Journal of Statistical Planning and Inference,
#'  139, 3419 - 3429.
#'
#'@export
#'
cistandard2 <- function(X, a, y, alpha){

  # Do the QR decomposition of the X matrix and find X transpose X inverse
  qrstr <- qr(X)
  R <- qr.R(qrstr)
  XTXinv <- chol2inv(R)

  # Find beta hat, theta hat and tau hat
  beta.hat <- XTXinv %*% t(X) %*% y
  theta.hat <- as.numeric(t(a) %*% beta.hat)

  # Find sigma hat
  n <- dim(X)[1]
  p <- dim(X)[2]
  sigsq <- (t(y - X %*% beta.hat) %*% (y - X %*% beta.hat)) / (n - p)
  sig <- as.numeric(sqrt(sigsq))

  # Find variance of theta hat on sigma squared
  v.theta <- as.numeric(t(a) %*% XTXinv %*% a)

  # Find the standard confidence interval
  m <- n - p
  t.alpha <- stats::qt(1 - alpha/2, m)
  standard.ci <- theta.hat + c(-1, 1) * sig * sqrt(v.theta) * t.alpha

  data.frame(lower = standard.ci[1], upper = standard.ci[2])

}
