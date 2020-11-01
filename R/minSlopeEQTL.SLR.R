# modified on June 21, 2020
#  (1) use exact power calculation formula
#
# created on Dec. 8, 2016
#  (1) minimum slope calculation for eQTL analysis based on simple linear regression
#

# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# power - desired power
# sigma.y - standard deviation of the outcome


diffPower4slope.SLR=function(slope,
                          MAF,
                          n,
                          power = 0.8,
                          sigma.y=0.13,
                          FWER=0.05,
                          nTests=200000)
{
  power.est = powerEQTL.SLR.default(MAF = MAF,
                            slope=slope,
                            n=n,
                            sigma.y=sigma.y,
                            FWER=FWER,
                            nTests=nTests)
  diff=power.est-power
  return(diff)
}
  

minSlopeEQTL.SLR=function(MAF,
                   n=200,
                   power=0.8,
                   sigma.y=0.13,
                   FWER=0.05,
                   nTests=200000)
{
  if(n <= 2)
  {
    stop("n must be > 2")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }
  
  sigma.x=sqrt(2*MAF*(1-MAF))
  # Dupont and Plummer (1998) Formula (1): slope = rho * sigma_y / sigma_x
  # where rho is the corr(y, x). So 0 < rho < 1 if we only consider positive correlation.
  # That is, 0 < slope < sigma_y/sigma_x
  res.root=uniroot(f=diffPower4slope.SLR,
                   interval = c(0.000001, sigma.y/sigma.x-0.000001),
                   MAF = MAF,
                   n= n,
                   power = power,
                   sigma.y=sigma.y,
                   FWER=FWER,
                   nTests=nTests)
  
  return(res.root$root)
}

