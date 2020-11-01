# modified on June 21, 2020
#  (1) power calculation for eQTL analysis based on simple linear regression
#      using exact formula
#

# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# sigma.y - standard deviation of the outcome


powerEQTL.SLR=function(MAF,
                       slope=0.13,
                       n=200,
                       power = NULL,
                       sigma.y=0.13,
                       FWER=0.05,
                       nTests=200000,
                       n.lower = 2.01,
                       n.upper = 1e+30)
{
  if(is.null(MAF)==TRUE && is.null(slope) == FALSE &&
     is.null(n) == FALSE && is.null(power) == FALSE)
  {
    MAF = minMAFeQTL.SLR(slope=slope,
                            n=n,
                            power=power,
                            sigma.y=sigma.y,
                            FWER=FWER,
                            nTests=nTests)
    names(MAF) = "MAF"
    return(MAF)
  } else if (is.null(MAF)==FALSE && is.null(slope) == TRUE &&
             is.null(n) == FALSE && is.null(power) == FALSE) {
    slope = minSlopeEQTL.SLR(MAF = MAF,
                             n= n,
                             power=power,
                             sigma.y=sigma.y,
                             FWER=FWER,
                             nTests=nTests)
    names(slope) = "slope"
    return(slope)
  } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
             is.null(n) == TRUE && is.null(power) == FALSE) {
    n = ssEQTL.SLR(MAF = MAF,
                        slope=slope,
                        power=power,
                        sigma.y=sigma.y,
                        FWER=FWER,
                        nTests=nTests,
                        n.lower = n.lower,
                        n.upper = n.upper)
    names(n)="n"
    return(n)
  } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
             is.null(n) == FALSE && is.null(power) == TRUE) {
    power = powerEQTL.SLR.default(MAF = MAF,
                                  slope=slope,
                                  n=n,
                                  sigma.y=sigma.y,
                                  FWER=FWER,
                                  nTests=nTests) 
    names(power) = "power"
    return(power)
  } else {
    stop("One and only one of the 4 parameters (MAF, slope, n, power) can be NULL!\n")
  }
  
}
  

powerEQTL.SLR.default=function(MAF,
                          slope=0.13,
                          n=200,
                          sigma.y=0.13,
                          FWER=0.05,
                          nTests=200000)
{
  sigma2.x = 2*MAF*(1-MAF)
  sigma2.y = sigma.y^2
  alpha = FWER/nTests

  delta = slope
  
  bound = sigma.y/sqrt(sigma2.x)

  if(delta >= bound || delta <= - bound)
  {
    stop("slope must be in the interval (-a, a), where a = sigma.y/sqrt(2MAF(1-MAF))!\n")
  }
  
  numer.ncp = delta *sqrt((n-1)*sigma2.x)
  denom.ncp= sqrt(sigma2.y - delta^2*sigma2.x)
  
  lambda = numer.ncp / denom.ncp
  
  mydf = n - 2
  cutoff = qt(1-alpha/2, df=mydf, ncp=0)
  
  power = 1 - pt(cutoff, df=mydf, ncp=lambda)
  power = power + pt(-cutoff, df=mydf, ncp=lambda)
  
  return(power)

}

