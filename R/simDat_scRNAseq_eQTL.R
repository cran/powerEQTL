# modified on Feb. 8, 2021
#  (1) remove un-used input option 'rho'
#
# created on Feb. 7, 2021
#  (1) simulate scRNAseq - genotype data
#  (2) counts follow zero-inflated negative binomial distribution:
#         Y_{ij} = 0 with prob p; Y_{ij} = NB(mu_i, theta) with prob 1-p
#         i = 1, ..., nSubj; j = 1, ..., nCellPersubj
#   mu_i = exp(int_i + slope * geno_i) 
#   int_i ~ N(int, sigma^2)


# one gene and one SNP;
# all cells in a subject has the same SNP genotypes;
# generate counts
simDat.eQTL.scRNAseq = function(nSubj = 50,
  nCellPerSubj = 100,
  # probability that an excess zero occurs
  zero.p = 0.01,
  # mean intercept for genotype in NB distribution: 
  #  mu_i = exp(int_i + slope * geno_i)
  # int_i ~ N(int, sigma^2)
  m.int = 0, # mean of the random intercept
  # standard deviation of the random intercept
  sigma.int = 1,
  # slope for genotype in NB distribution: 
  #  mu_i = exp(int_i + slope * geno_i)
  # int_i ~ N(int, sigma^2)
  slope = 1,
  # dispersion parameter of NB distribution
  # the smaller theta is, the larger variance of NB random variable is
  theta  = 1,
  MAF = 0.45 # SNP MAF
)
{
  # N = nSubj*nCellPerSubj
  N = nSubj * nCellPerSubj

  id = rep(seq_len(nSubj), each = nCellPerSubj)

  # generate genotypes
  # under Hardy-Weinberg equilibrium
  # p(AA) = MAF^2
  # p(AB) = 2*MAF*(1-MAF)
  # p(BB) = (1-MAF)^2
  # A = minor allele, B = major allele
  pAA = MAF^2
  pAB = 2*MAF*(1-MAF)
  pBB = (1-MAF)^2

  geno = rep(sample(x=c(2, 1, 0), size = nSubj, replace = TRUE,
	 prob = c(pAA, pAB, pBB)), each = nCellPerSubj)

  ##
  counts = rep(0, N)
  NBflag=sample(x=c(0,1), size=N, prob = c(zero.p, 1-zero.p), replace=TRUE)
  pos.NB = which(NBflag == 1)
  N.NB = length(pos.NB)


  # generate counts from NB distribution
  intercept = rep(rnorm(n=nSubj, mean = m.int, sd = sigma.int), each = nCellPerSubj)
  mu = exp(intercept + slope * geno)
  # The variance is mu + mu^2/size in this parametrization
  # the smaller size is, the larger variance of NB random variable is
  counts[pos.NB] = rnbinom(n = N.NB, size = theta, mu = mu)

  ##
  frame = data.frame(id = id, geno = geno, counts = counts)

  invisible(frame)

}

