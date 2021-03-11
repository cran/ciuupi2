choice_d <- function(m, length.out = 60){
  # Compute the value of d and the vector
  # (0:d+constant) of length 60 for a given m.
  #
  # Inputs:
  # m: degrees of freedom n - p
  # length.out: length of gams vector
  #
  # Output:
  # A list which contains the value d and vector gams.
  #
  # Written by P Kabaila in September 2020
  # Modified by N Ranathunga in October 2020

  # Set a cutoff value to the Normal curve
  init.cutoff.d <- 1.545

  # Find a value for d
  init.prob.d <- stats::pnorm(init.cutoff.d, 0, 1)
  multiplier <- 6 / init.cutoff.d
  d <- round(multiplier * stats::qt(init.prob.d, m), 1)

  # Find a maximum possible value for gams
  init.prob.extra <- stats::pnorm(2, 0, 1)
  extra <- round(stats::qt(init.prob.extra, m), 1)
  max.gamma.constr <- d + extra

  # Find the gams vector
  gams <- c(seq(0, d, length.out = length.out),
            seq(d, max.gamma.constr, by = (2 * d / length.out)))

  out <- list(d=d, gams=gams)

}
