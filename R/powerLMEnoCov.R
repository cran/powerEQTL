# created on Apr 22, 2020
#  (1) power for simple linear mixed effects model without covariate

#
# We assume the following linear mixed effects model to characterize
#   the association between x and y:
#  d_{ij} = beta_{0i} + epsilon_{ij},
#    beta_{0i} ~ N(beta_0, sigma^2_{\beta})
#    epsilon_{ij} ~ N(0, sigma^2_{epsilon})
#
# delta - slope under alternative hypothesis
# n - number of subjects
# m - number of observations per subject
# sigma.y - standard deviation of the outcome y
# rho - intra-class correlation (i.e., correlation between y_{ij} and y_{ik})
#        rho = sigma^2_{beta} / (sigma^2_{beta}+sigma^2_{epsilon})
# alpha - type I error rate 
# nTests = number of tests
powerLMEnoCov=function(
  delta, 
  n, 
  m, 
  sigma.y, 
  rho=0.8, 
  alpha=0.05,
  nTests=1)
{
  alpha2=alpha/nTests
  za2=qnorm(1-alpha2/2)
  
  part0=delta*sqrt(m*n)/(sigma.y*sqrt(1+(m-1)*rho))
  part1 = za2-part0
  
  part2 = -za2-part0
  
  power = 1- pnorm(part1) + pnorm(part2)
  
  return(power)
  
}
