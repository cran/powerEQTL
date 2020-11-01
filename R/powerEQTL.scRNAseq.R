# modified on June 22, 2020
#  (1) rename 'alpha' in powerEQTL.scRNAseq' by 'FWER'
#  (2) rename 'theta' in powerEQTL.scRNAseq' by 'MAF'
#  (3) rename 'delta' in powerEQTL.scRNAseq' by 'slope'
#
# modified on Apr 7, 2020
#  (1) changed the function name power.eQTL.scRNAseq to 
#      powerEQTL.scRNAseq
#  (2) simplify R code
#
# created on June 23, 2019
#  (1) define function to calculate power for eQTL based on scRNAseq

# define function
# We assume the following linear mixed effects model to characterize
#   the association between genotype and gene expression:
#  y_{ij} = beta_{0i} + beta_1 * x_i + epsilon_{ij},
#    beta_{0i} ~ N(beta_0, sigma^2_{\beta})
#    epsilon_{ij} ~ N(0, sigma^2)
#
# slope - slope under alternative hypothesis
# n - number of subjects
# m - number of cells per subject
# sigma.y - standard deviation of the gene expression
# MAF - minor allele frequency (between 0 and 0.5)
# rho - intra-class correlation (i.e., correlation between y_{ij} and y_{ik})
#        rho = sigma^2_{beta} / (sigma^2_{beta}+sigma^2)
# FWER - family-wise type I error rate
# nTests = number of genes * number of SNPs


powerEQTL.scRNAseq=function(
  slope, 
  n, 
  m, 
  power = NULL,
  sigma.y, 
  MAF=0.2, 
  rho=0.8, 
  FWER=0.05,
  nTests=1,
  n.lower=2.01,
  n.upper=1e+30)
{
  if(is.null(MAF)==TRUE && is.null(slope) == FALSE &&
     is.null(n) == FALSE && is.null(power) == FALSE)
  {
    MAF = minMAFEQTL.scRNAseq(
      slope = slope,
      n = n, 
      m = m,
      power = power,
      sigma.y = sigma.y, 
      rho=rho, 
      FWER=FWER,
      nTests=nTests)

    names(MAF) = "MAF"
    return(MAF)
  } else if (is.null(MAF)==FALSE && is.null(slope) == TRUE &&
             is.null(n) == FALSE && is.null(power) == FALSE) {
    slope = minSlopeEQTL.scRNAseq(
      n = n, 
      m = m,
      power = power,
      sigma.y = sigma.y, 
      MAF=MAF, 
      rho=rho, 
      FWER=FWER,
      nTests=nTests)

    names(slope) = "slope"
    return(slope)
  } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
             is.null(n) == TRUE && is.null(power) == FALSE) {
    n = ssEQTL.scRNAseq(  
          slope = slope, 
          m = m,
          power = power,
          sigma.y = sigma.y, 
          MAF=MAF, 
          rho=rho, 
          FWER=FWER,
          nTests=nTests,
          n.lower=n.lower,
          n.upper=n.upper)

    names(n)="n"
    return(n)
  } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
             is.null(n) == FALSE && is.null(power) == TRUE) {
    power = powerEQTL.scRNAseq.default(
      slope = slope, 
      n = n, 
      m = m, 
      sigma.y = sigma.y, 
      MAF=MAF, 
      rho=rho, 
      FWER=FWER,
      nTests=nTests)
      
    names(power) = "power"
    return(power)
  } else {
    stop("One and only one of the 4 parameters (MAF, slope, n, power) can be NULL!\n")
  }
  
}


powerEQTL.scRNAseq.default=function(
  slope, 
  n, 
  m, 
  sigma.y, 
  MAF=0.2, 
  rho=0.8, 
  FWER=0.05,
  nTests=1)
{
  alpha2=FWER/nTests
  za2=qnorm(1-alpha2/2)
  
  sigma.x=sqrt(2*MAF*(1-MAF))
  part0=sigma.x*slope*sqrt(m*(n-1))/(sigma.y*sqrt(1+(m-1)*rho))
  part1 = za2-part0
  
  part2 = -za2-part0
  
  power = 1- pnorm(part1) + pnorm(part2)
  
  return(power)
  
}
