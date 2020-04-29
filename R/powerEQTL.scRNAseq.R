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
#    epsilon_{ij} ~ N(0, sigma^2_{epsilon})
#
# delta - slope under alternative hypothesis
# n - number of subjects
# m - number of cells per subject
# sigma.y - standard deviation of the gene expression
# theta - minor allele frequency (between 0 and 0.5)
# rho - intra-class correlation (i.e., correlation between y_{ij} and y_{ik})
#        rho = sigma^2_{beta} / (sigma^2_{beta}+sigma^2_{epsilon})
# alpha - type I error rate
# nTests = number of genes * number of SNPs
powerEQTL.scRNAseq=function(
  delta, 
  n, 
  m, 
  sigma.y, 
  theta=0.2, 
  rho=0.8, 
  alpha=0.05,
  nTests=1)
{
  alpha2=alpha/nTests
  za2=qnorm(1-alpha2/2)
  
  sigma.x=sqrt(2*theta*(1-theta))
  part0=sigma.x*delta*sqrt(m*(n-1))/(sigma.y*sqrt(1+(m-1)*rho))
  part1 = za2-part0
  
  part2 = -za2-part0
  
  power = 1- pnorm(part1) + pnorm(part2)
  
  return(power)
  
}
