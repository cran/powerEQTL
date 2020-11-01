# modified on June 22, 2020
#  (1) rename 'alpha' in powerEQTL.scRNAseq' by 'FWER'
#  (3) rename 'delta' in powerEQTL.scRNAseq' by 'slope'
#
# created on Apr 7, 2020
#  (1) power for simple linear mixed effects model
#
# We assume the following linear mixed effects model to characterize
#   the association between x and y:
#  y_{ij} = beta_{0i} + beta_1 * x_i + epsilon_{ij},
#    beta_{0i} ~ N(beta_0, sigma^2_{\beta})
#    epsilon_{ij} ~ N(0, sigma^2)
#
# slope - slope under alternative hypothesis
# n - number of subjects
# m - number of observations per subject
# sigma.y - standard deviation of the outcome y
# sigma.x - standard deviation of the predictor x
# rho - intra-class correlation (i.e., correlation between y_{ij} and y_{ik})
#        rho = sigma^2_{beta} / (sigma^2_{beta}+sigma^2)
# FWER - family-wise type I error rate 
# nTests = number of tests

powerLME.default=function(
  slope, 
  n, 
  m, 
  sigma.y, 
  sigma.x, 
  rho=0.8, 
  FWER=0.05,
  nTests=1)
{
  alpha2=FWER/nTests
  za2=qnorm(1-alpha2/2)
  
  part0=sigma.x*slope*sqrt(m*(n-1))/(sigma.y*sqrt(1+(m-1)*rho))
  part1 = za2-part0
  
  part2 = -za2-part0
  
  power = 1- pnorm(part1) + pnorm(part2)
  
  return(power)
  
}

powerLME=function(
  slope, 
  n, 
  m, 
  sigma.y, 
  sigma.x,
  power = NULL,
  rho=0.8, 
  FWER=0.05,
  nTests=1,
  n.lower = 2.01,
  n.upper = 1e+30)
{

  if (is.null(slope) == TRUE &&
             is.null(n) == FALSE && is.null(power) == FALSE) {
    slope = minSlope.LME(n = n, 
                         m = m, 
                         sigma.y = sigma.y, 
                         sigma.x = sigma.x,
                         power = power,
                         rho=rho, 
                         FWER=FWER,
                         nTests=nTests
    )
    
    names(slope) = "slope"
    return(slope)
  } else if (is.null(slope) == FALSE &&
             is.null(n) == TRUE && is.null(power) == FALSE) {
    n = ssLME(  slope = slope, 
                         m = m, 
                         sigma.y = sigma.y, 
                         sigma.x = sigma.x,
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
    power = powerLME.default(
      slope = slope, 
      n = n, 
      m = m, 
      sigma.y = sigma.y, 
      sigma.x = sigma.x, 
      rho=rho, 
      FWER=FWER,
      nTests=nTests)
    
    names(power) = "power"
    return(power)
  } else {
    stop("One and only one of the 3 parameters (slope, n, power) can be NULL!\n")
  }  
  
}
