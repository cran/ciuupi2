#' Find rho
#'
#' Find the correlation rho for given \eqn{n} by \eqn{p} design matrix X and
#' given \eqn{p}-vectors a and c
#'
#' @param X The \eqn{n} by \eqn{p} design matrix
#' @param a A vector used to specify the parameter of interest
#' @param c A vector used to specify the parameter about which we have uncertain
#'   prior information
#'
#' @return The value of the correlation rho.
#'
#' @details Suppose that \deqn{y = X \beta + \epsilon} where \eqn{y} is a random
#'   \eqn{n}-vector of responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix
#'   with linearly independent columns, \eqn{\beta} is an unknown parameter
#'   \eqn{p}-vector and \eqn{\epsilon} is a random \eqn{n}-vector with
#'   components that are independent and identically normally distributed with
#'   zero mean and unknown variance. The parameter of interest is \eqn{\theta =
#'   } \code{a}' \eqn{\beta}. The uncertain prior information is that \eqn{\tau
#'   = } \code{c}' \eqn{\beta} takes the value \code{t}, where \code{a} and
#'   \code{c} are specified linearly independent nonzero \eqn{p}-vectors and
#'   \code{t} is a specified number. \code{rho} is the known correlation between
#'   the least squares estimators of \eqn{\theta} and \eqn{\tau}. It is
#'   determined by the \eqn{n} by \eqn{p} design matrix X and the
#'   \eqn{p}-vectors a and c.
#'
#' @section \eqn{X}, \code{a} and \code{c} for a particular example: Consider
#'   the same 2 x 2 factorial example as that described in Section 4 of Kabaila
#'   and Giri (2009), except that the number of replicates is 3 instead of 20.
#'   In this case, \eqn{X} is a 12 x 4 matrix, \eqn{\beta} is an unknown
#'   parameter 4-vector and \eqn{\epsilon} is a random 12-vector with components
#'   that are independent and identically normally distributed with zero mean
#'   and unknown variance. In other words, the length of the response vector
#'   \eqn{y} is \eqn{n} = 12 and the length of the parameter vector \eqn{\beta}
#'   is \eqn{p} = 4, so that \eqn{m = n - p} = 8. The parameter of interest is
#'   \eqn{\theta = } \code{a}' \eqn{\beta}, where the column vector \code{a} =
#'   (0, 2, 0, -2). Also, the parameter \eqn{\tau = } \code{c}' \eqn{\beta},
#'   where the column vector \code{c} = (0, 0, 0, 1). The uncertain prior
#'   information is that \eqn{\tau = } \code{t}, where \code{t} = 0.
#'
#'   The design matrix \eqn{X} and the vectors \code{a} and \code{c} (denoted in
#'   R by a.vec and c.vec, respectively) are entered into R using the commands
#'   in the following example.
#'
#' @references Kabaila, P. and Giri, R. (2009).  Confidence intervals in
#'   regression utilizing prior information. Journal of Statistical Planning and
#'   Inference, 139, 3419-3429.
#'
#' @examples
#' col1 <- rep(1,4)
#' col2 <- c(-1, 1, -1, 1)
#' col3 <- c(-1, -1, 1, 1)
#' col4 <- c(1, -1, -1, 1)
#' X.single.rep <- cbind(col1, col2, col3, col4)
#' X <- rbind(X.single.rep, X.single.rep, X.single.rep)
#' a.vec <- c(0, 2, 0, -2)
#' c.vec <- c(0, 0, 0, 1)
#'
#' # Find the value of rho
#' rho <- find_rho(X, a=a.vec, c=c.vec)
#' rho
#'
#' # The value of rho is -0.7071068
#'
#' @export
find_rho <- function(X, a, c){

# Do the QR decomposition of the matrix X and find X transpose X inverse
qrstr <- qr(X)
R <- qr.R(qrstr)
XTXinv <- chol2inv(R)

# Compute rho
rho <- (t(a) %*% XTXinv %*% c) /
  sqrt( t(a) %*% XTXinv %*% a %*% t(c) %*% XTXinv %*% c)
rho <- as.numeric(rho)

}
