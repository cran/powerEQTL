# created on June 21, 2020
#  (1) write a wrapper function
#  (2) rename input parameters
#  (3) remove powerEQTL.ANOVA2 related functions

# created on Dec. 11, 2016
# (1) allow mu2-mu1 not equal to mu3-mu2,
#   where mu1 is the mean expression level for mutation homozygote,
# mu2 is the mean expression level for heterozygote,
#  and mu3 is the mean expression level for wildtype homozygote
#
# created on Dec. 8, 2016
#  (1) power calculation for eQTL analysis based on ANOVA or simple linear regression
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


powerEQTL.ANOVA=function(MAF,
                         deltaVec=c(-0.13, 0.13),
                       n=200,
                       power = NULL,
                       sigma=0.13,
                       FWER=0.05,
                       nTests=200000,
                       n.lower = 4,
                       n.upper = 1e+30)
{
  if(is.null(MAF)==TRUE && 
     is.null(n) == FALSE && is.null(power) == FALSE)
  {
    MAF = minMAFeQTL.ANOVA(
      deltaVec=deltaVec,
      n=n,
      power = power,
      sigma=sigma,
      FWER=FWER,
      nTests=nTests)
    
    names(MAF) = "MAF"
    return(MAF)
  } else if (is.null(MAF)==FALSE &&
             is.null(n) == TRUE && is.null(power) == FALSE) {
    n = ssEQTL.ANOVA(MAF = MAF,
                     deltaVec=deltaVec,
                     power=power,
                     sigma=sigma,
                     FWER=FWER,
                     nTests=nTests,
                     n.lower = n.lower,
                     n.upper = n.upper)
      
    names(n)="n"
    return(n)
  } else if (is.null(MAF)==FALSE && 
             is.null(n) == FALSE && is.null(power) == TRUE) {
    power = powerEQTL.ANOVA.default(MAF = MAF,
                                    deltaVec=deltaVec,
                                    n=n,
                                    sigma=sigma,
                                    FWER=FWER,
                                    nTests=nTests) 
    names(power) = "power"
    return(power)
  } else {
    stop("One and only one of the 3 parameters (MAF, n, power) can be NULL!\n")
  }
  
}


powerEQTL.ANOVA.default=function(MAF,
                   deltaVec=c(-0.13, 0.13),
                   n=200,
                   sigma=0.13,
                   FWER=0.05,
                   nTests=200000)
{
  if(length(deltaVec)!=2)
  {
    stop("'deltaVec' has 2 and only 2 elements!\n1st element = mu2 - mu1; 2nd element = mu3 - mu2!\n")
  }
  
  gm1 = -deltaVec[1] # mu2 - mu1
  gm2 = 0
  gm3 = deltaVec[2] # mu3 - mu2

  w1=MAF^2 # mutation homozygotes
  w2=2*MAF*(1-MAF) # heterozygotes
  w3=(1-MAF)^2 # wildtype homozygotes


  alpha=FWER/nTests

  k=3
  mydf1=k-1
  mydf2=n-k
  q=qf(p=1-alpha, df1=mydf1, df2=mydf2)

  wVec=c(w1, w2, w3)
  muVec=c(gm1, gm2, gm3)

  mu=sum(wVec*muVec, na.rm=TRUE)

  myncp = n*sum(wVec*(muVec-mu)^2, na.rm=TRUE)
  myncp=myncp/(sigma^2)

  power=1-pf(q=q, df1=mydf1, df2=mydf2,
             ncp=myncp)

  return(power)

}



