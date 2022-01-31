## Utility functions for analysis of GRN output

library(e1071)

# calculate_bimodality calculates Sarle's bimodality coefficient for a finite sample (https://en.wikipedia.org/wiki/Multimodal_distribution#Bimodality_coefficient).
# Note that a unimodal distribution will tend toward 0, a uniform distribution will give a value of ~0.555, and a Bernoulli distribution will give a value of 1.
# Accepts:
# - x: a numeric vector. Can be any length, but will return a warning if less than length 4
#
# Returns:
# - bc: a numeric value of Sarle's bimodality coefficient, between 0 and 1. Will return NA if length<4. Will return NaN if skewness or kurtosis is not calculable (e.g., for all values 0).
calculate_bimodality <- function(x) {
  
  if(length(x) < 4) {
    return(NA)
  }
  
  sx <- skewness(x)
  kx <- kurtosis(x) # excess kurtosis = m_4 / m_2^2 - 3
  nx <- length(x)
  
  bc <- (sx^2 + 1) / (kx + (3*(nx-1)^2)/((nx-2)*(nx-3)))
  
  return(bc)
  
}
