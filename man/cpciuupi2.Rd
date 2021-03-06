% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpciuupi2.R
\name{cpciuupi2}
\alias{cpciuupi2}
\title{Compute the coverage probability of the Kabaila & Giri (2009) CIUUPI}
\usage{
cpciuupi2(gam, bsvec, alpha, m, rho, natural = 1)
}
\arguments{
\item{gam}{A value of gamma or vector of gamma values at which the coverage
probability function is evaluated}

\item{bsvec}{The vector (b(d/6),b(2d/6),...,b(5d/6),s(0),s(d/6),...,s(5d/6))
computed using \code{bsciuupi2}}

\item{alpha}{The minimum coverage probability is 1 - \code{alpha}}

\item{m}{Degrees of freedom \code{n - p}}

\item{rho}{A known correlation}

\item{natural}{Equal to 1 (default) if the b and s functions are obtained by
natural cubic spline interpolation or 0 if obtained by clamped cubic spline
interpolation. This parameter must take the same value as that used in
\code{bsciuupi2}}
}
\value{
The value(s) of the coverage probability of the Kabaila & Giri (2009)
 CIUUPI at \code{gam}.
}
\description{
Evaluate the coverage probability of the Kabaila & Giri (2009) confidence
interval that utilizes uncertain prior information (CIUUPI),
with minimum coverage 1 - \code{alpha}, at \code{gam}.
}
\details{
Suppose that \deqn{y = X \beta + \epsilon} where \eqn{y} is a random
 \eqn{n}-vector of responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix
 with linearly independent columns, \eqn{\beta} is an unknown parameter
 \eqn{p}-vector and \eqn{\epsilon} is a random \eqn{n}-vector with
 components that are independent and identically normally distributed with
 zero mean and unknown variance. The parameter of interest is \eqn{\theta = }
 \code{a}' \eqn{\beta}. The uncertain prior information is that \eqn{\tau = }
 \code{c}' \eqn{\beta} takes the value \code{t}, where \code{a} and \code{c}
 are specified linearly independent vectors and \code{t} is a specified
 number. \code{rho} is the known correlation between the least squares
 estimators of \eqn{\theta} and \eqn{\tau}. It is determined by the \eqn{n}
 by \eqn{p} design matrix X and the \eqn{p}-vectors a and c using
 \code{\link{find_rho}}.

 In the examples, we continue with the same 2 x 2 factorial example described
 in the documentation for \code{\link{find_rho}}.
}
\examples{
alpha <- 0.05
m <- 8

# Find the vector (b(d/6),...,b(5d/6),s(0),...,s(5d/6)) that specifies the
# Kabaila & Giri (2009) CIUUPI for the first definition of the
# scaled expected length (default) (takes about 30 mins to run):
\donttest{
bsvec <- bsciuupi2(alpha, m, rho = -0.7071068)
}

# The result bsvec is (to 7 decimal places) the following:
bsvec <- c(-0.0287487, -0.2151595, -0.3430403, -0.3125889, -0.0852146,
            1.9795390,  2.0665414,  2.3984471,  2.6460159,  2.6170066,  2.3925494)


# Graph the coverage probability function
gam <- seq(0, 10, by = 0.1)
cp <- cpciuupi2(gam, bsvec, alpha, m, rho = -0.7071068)
plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
main = "Coverage Probability", col = "blue",
xlab = expression(paste("|", gamma, "|")), ylim = c(0.9490, 0.9510))
abline(h = 1-alpha, lty = 2)

}
\references{
Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
prior information.  Journal of Statistical Planning and Inference, 139, 3419 -
3429.
}
\seealso{
\code{\link{find_rho}}, \code{\link{bsciuupi2}}
}
