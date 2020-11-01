# created on June 24, 2020
#  (1) add wrapper function
#  (2) rename 'delta' to 'slope'
#  (3) rename 'alpha' to 'FWER'

# created on Apr 22, 2020
#  (1) power for simple linear mixed effects model without covariate

#
# We assume the following linear mixed effects model to characterize
#   the association between x and y:
#  d_{ij} = beta_{0i} + epsilon_{ij},
#    beta_{0i} ~ N(beta_0, sigma^2_{\beta})
#    epsilon_{ij} ~ N(0, sigma^2)
#
# slope - slope under alternative hypothesis
# n - number of subjects
# m - number of observations per subject
# sigma.y - standard deviation of the outcome y
# rho - intra-class correlation (i.e., correlation between y_{ij} and y_{ik})
#        rho = sigma^2_{beta} / (sigma^2_{beta}+sigma^2)
# FWER - family-wise type I error rate 
# nTests = number of tests
powerLMEnoCov.default=function(
  slope, 
  n, 
  m, 
  sigma.y, 
  rho=0.8, 
  FWER=0.05,
  nTests=1)
{
  alpha2=FWER/nTests
  za2=qnorm(1-alpha2/2)
  
  part0=slope*sqrt(m*n)/(sigma.y*sqrt(1+(m-1)*rho))
  part1 = za2-part0
  
  part2 = -za2-part0
  
  power = 1- pnorm(part1) + pnorm(part2)
  
  return(power)
  
}


powerLMEnoCov=function(
  slope, 
  n, 
  m, 
  sigma.y,
  power = NULL,
  rho=0.8, 
  FWER=0.05,
  nTests=1,
  n.lower = 2.01,
  n.upper = 1e+30)
{
  if (is.null(slope) == TRUE &&
      is.null(n) == FALSE && is.null(power) == FALSE) {
    slope = minSlopeLMEnoCov(n=n, 
                                      m=m, 
                                      sigma.y=sigma.y,
                                      power = power,
                                      rho=rho, 
                                      FWER=FWER,
                                      nTests=nTests,
                                      n.lower = n.lower,
                                      n.upper = n.upper
    )
    
    names(slope) = "slope"
    return(slope)
  } else if (is.null(slope) == FALSE &&
             is.null(n) == TRUE && is.null(power) == FALSE) {
    n = ssLMEnoCov(    slope=slope, 
                                m=m, 
                                sigma.y=sigma.y,
                                power = power,
                                rho=rho, 
                                FWER=FWER,
                                nTests=nTests,
                                n.lower = n.lower,
                                n.upper = n.upper
    )
    
    names(n)="n"
    return(n)
  } else if (is.null(slope) == FALSE &&
             is.null(n) == FALSE && is.null(power) == TRUE) {
    power = powerLMEnoCov.default(
      slope = slope, 
      n = n, 
      m = m, 
      sigma.y = sigma.y, 
      rho=rho, 
      FWER=FWER,
      nTests=nTests)
    
    names(power) = "power"
    return(power)
  } else {
    stop("One and only one of the 3 parameters (slope, n, power) can be NULL!\n")
  }    
}
