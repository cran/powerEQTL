# modified on June 21, 2020
#  (1) use exact power calculation formula
#
# created on Jan. 15, 2017
#  (1) minimum MAF calculation for eQTL analysis based on simple linear regression
#


# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# power - desired power
# sigma.y - standard deviation of the outcome


diffPower4MAF.SLR=function(MAF,
                             slope,
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


minMAFeQTL.SLR=function(slope,
                   n=200,
                   power=0.8,
                   sigma.y=0.13,
                   FWER=0.05,
                   nTests=200000)
{
  a = 1/4 - sigma.y^2/(2*slope^2)
  if(a > 0)
  {
    upp.MAF = 0.5 - sqrt(1/4 - sigma.y^2/(2*slope^2))
  } else {
    upp.MAF = 0.5
  }
  upp.MAF = upp.MAF - (1.0e-6)
  
  res.root=uniroot(f=diffPower4MAF.SLR,
    interval=c(0.000001, upp.MAF - 1.0e-6),
    slope=slope,   
    n=n,
    sigma.y=sigma.y,
    FWER=FWER,
    nTests=nTests,
    power=power)
                        
  
  return(res.root$root)
}

