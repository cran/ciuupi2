% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ciuupi2.R
\name{ciuupi2}
\alias{ciuupi2}
\title{Compute the Kabaila & Giri (2009) CIUUPI}
\usage{
ciuupi2(alpha, X, a, c, bsvec, t, y, natural = 1)
}
\arguments{
\item{alpha}{1 - \code{alpha} is the minimum coverage probability of the
confidence interval}

\item{X}{The \eqn{n} by \eqn{p} design matrix}

\item{a}{A vector used to specify the parameter of interest}

\item{c}{A vector used to specify the parameter about which we have uncertain
prior information}

\item{bsvec}{The vector (b(d/6),b(2d/6),...,b(5d/6),s(0),s(d/6),...,s(5d/6))
computed using \code{bsciuupi2}}

\item{t}{A number used to specify the uncertain prior information, which has
the form \eqn{\tau =} t}

\item{y}{The \eqn{n}-vector of observed responses}

\item{natural}{Equal to 1 (default) if the b and s functions are evaluated by
natural cubic spline interpolation or 0 if evaluated by clamped cubic spline
interpolation. This parameter must take the same value as that used in
\code{bsciuupi2}}
}
\value{
The Kabaila & Giri (2009) confidence interval, with minimum coverage 1
 - \code{alpha}, that utilizes the uncertain prior information.
}
\description{
Compute the Kabaila and Giri (2009) confidence interval that utilizes
uncertain prior information (CIUUPI), with minimum coverage 1 - \code{alpha},
for a given vector \code{y} of observed responses.
}
\details{
Suppose that \deqn{y = X \beta + \epsilon} where \eqn{y} is a random
 \eqn{n}-vector of responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix
 with linearly independent columns, \eqn{\beta} is an unknown parameter
 \eqn{p}-vector and \eqn{\epsilon} is a random \eqn{n}-vector with components
 that are independent and identically normally distributed with zero mean and
 unknown variance. The parameter of interest is \eqn{\theta = } \code{a}'
 \eqn{\beta}. The uncertain prior information is that \eqn{\tau = } \code{c}'
 \eqn{\beta} takes the value \code{t}, where \code{a} and \code{c} are
 specified linearly independent vectors nonzero \eqn{p}-vectors and \code{t}
 is a specified number. Given the vector \code{bsvec}, computed using
 \code{\link{bsciuupi2}}, the design matrix \code{X}, the vectors \code{a}
 and \code{c} and the number t, \code{ciuupi2} computes the confidence
 interval for \eqn{\theta} that utilizes the uncertain prior information that
 \eqn{\tau} = \code{t} for given \eqn{n}-vector of observed responses
 \code{y}.

 In the examples, we continue with the same 2 x 2 factorial example described
 in the documentation for \code{\link{find_rho}}, for the vector of observed
 responses \eqn{y} = (-1.3, 0.8, 2.6, 5.8, 0.3, 1.3, 4.3, 5.0, -0.4, 1.0,
 5.2, 6.2).
}
\examples{
# Specify the design matrix X and vectors a and c
# (denoted in R by a.vec and c.vec, respectively)
col1 <- rep(1,4)
col2 <- c(-1, 1, -1, 1)
col3 <- c(-1, -1, 1, 1)
col4 <- c(1, -1, -1, 1)
X.single.rep <- cbind(col1, col2, col3, col4)
X <- rbind(X.single.rep, X.single.rep, X.single.rep)
a.vec <- c(0, 2, 0, -2)
c.vec <- c(0, 0, 0, 1)


# Compute the vector (b(d/6),...,b(5d/6),s(0),...,s(5d/6)) that specifies the
# Kabaila & Giri (2009) CIUUPI, with minimum coverage 1 - alpha, for the
# first definition of the scaled expected length (default)
# for given alpha, m and rho (takes about 30 mins to run):
\donttest{
bsvec <- bsciuupi2(alpha = 0.05, m = 8, rho = -0.7071068)
}

# The result bsvec is (to 7 decimal places) the following:
bsvec <- c(-0.0287487, -0.2151595, -0.3430403, -0.3125889, -0.0852146,
            1.9795390,  2.0665414,  2.3984471,  2.6460159,  2.6170066,  2.3925494)


# Specify t and y
t <- 0
y <- c(-1.3, 0.8, 2.6, 5.8, 0.3, 1.3, 4.3, 5.0, -0.4, 1.0, 5.2, 6.2)

# Find the Kabaila and Giri (2009) CIUUPI, with minimum coverage 1 - alpha,
# for the first definition of the scaled expected length
res <- ciuupi2(alpha=0.05, X, a=a.vec, c=c.vec, bsvec, t, y, natural = 1)
res

# The Kabaila and Giri (2009) CIUUPI, with minimum coverage 1 - alpha,
# is (0.14040, 2.85704).
# The usual 1 - alpha confidence interval for theta is (-0.08185, 3.08185).

}
\references{
Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
prior information.  Journal of Statistical Planning and Inference, 139, 3419 -
3429.
}
\seealso{
\code{\link{find_rho}}, \code{\link{bsciuupi2}}
}
