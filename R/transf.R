transf <- function(z){
  # This function applies the
  # transformation (2.6) of Mori(1988)
  # to the value z
  #
  # Inputs:
  # z: a given value
  #
  # Written by N. Ranathunga in September 2020

  out <- exp((z/2) - exp(-z))

}
