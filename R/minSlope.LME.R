# created on June 24, 2020

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

#            
# n.lower - lower bound for sample size
# n.upper - upper bound for sample size
# verbose - flag indicating if intermedaite results should be output

diffPower4slope.LME=function(  slope,
                            n,
                            m, 
                            sigma.y, 
                            sigma.x,
                            power = 0.8,
                            rho=0.8, 
                            FWER=0.05,
                            nTests=1)
{
  power.est = powerLME.default(
    slope = slope, 
    n = n, 
    m = m, 
    sigma.y = sigma.y, 
    sigma.x = sigma.x, 
    rho=rho, 
    FWER=FWER,
    nTests=nTests)
  
  diff=power.est-power
  return(diff)
  
}


minSlope.LME=function(n, 
                 m, 
                 sigma.y, 
                 sigma.x,
                 power = 0.8,
                 rho=0.8, 
                 FWER=0.05,
                 nTests=1
)
{
  res.uni=uniroot(f=diffPower4slope.LME,
                  interval = c(1.0e-6, 1.0e+30),
                  n = n,
                  m = m, 
                  sigma.y = sigma.y, 
                  sigma.x = sigma.x,
                  power = power,
                  rho=rho, 
                  FWER=FWER,
                  nTests=nTests
  )
  return(res.uni$root)
}

