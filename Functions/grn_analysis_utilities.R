## Utility functions for analysis of GRN output

library(e1071)

# calculate_bimodality calculates Sarle's bimodality coefficient for a finite sample (https://en.wikipedia.org/wiki/Multimodal_distribution#Bimodality_coefficient).
# Note that a unimodal distribution will tend toward 0, a uniform distribution will give a value of ~0.555, and a Bernoulli distribution will give a value of 1.
# Accepts:
# - x: a numeric vector.
#
# Returns:
# - bc: a numeric value of Sarle's bimodality coefficient, between 0 and 1. Will return 0 for spike-and-slab dist (e.g., when all values 0)
calculate_bimodality <- function(x) {
  
  nx <- length(x)
  
  if(sum(x == mean(x)) == nx) {
    bc = 0
  } else {
    
    sx <- skewness(x)
    kx <- kurtosis(x) # excess kurtosis = m_4 / m_2^2 - 3
    
    bc <- (sx^2 + 1) / (kx + (3*(nx-1)^2)/((nx-2)*(nx-3)))
  }
  return(bc)
  
}
