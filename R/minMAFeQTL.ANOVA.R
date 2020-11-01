# created on Jan. 15, 2017
# calculate MAF based on simple linear regression
#

# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# sigma - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu1-mu2 and deltaVec[2]=mu3-mu2
#  group 1 is mutation homozygotes
#  group 2 is heterozygotes
#  group 3 is wildtype homozygotes

diffPower4MAF.ANOVA=function(MAF,
                             deltaVec = c(-0.13, 0.13),
                             n = 200,
                             sigma = 0.13,
                             FWER=0.05,
                             nTests=200000,
                             power=0.8)
{
  est.power=powerEQTL.ANOVA.default(MAF = MAF,
                                    deltaVec=deltaVec,
                                    n=n,
                                    sigma=sigma,
                                    FWER=FWER,
                                    nTests=nTests)
  diff=(est.power-power)
  return(diff)

}


# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu2-mu1 and deltaVec[2]=mu3-mu2
# verbose - flag indicating if we should output intermediate results

minMAFeQTL.ANOVA=function(
                          deltaVec=c(-0.13, 0.13),
                          n=200,
                          power = 0.8,
                          sigma=0.13,
                          FWER=0.05,
                          nTests=200000)
{
  res.uni=uniroot(f=diffPower4MAF.ANOVA,
                  interval=c(0.000001, 0.5 - 1.0e-6),
                  deltaVec = deltaVec,
                  n = n,
                  sigma = sigma,
                  FWER=FWER,
                  nTests=nTests,
                  power=power
                  )
  
  return(res.uni$root)
}

