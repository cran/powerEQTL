# modified on July 21, 2021
#  (1) remove dependency on glmmTMB since it is no longer available in CRAN
#
# modified on June 17, 2021
#  (1) drop glmmADMB since it is not in mainstream repositories

# modified on June 13, 2021
#  (1) added glmmADMB
#
# modified on Feb. 15, 2021
#  (1) rename 'nSubj' back to 'n' and 'nCellPerSubj' back to 'm' to be consistent with other R functions in this package
#  (2) use asymptotic distribution of test statistics to calculate power to get stable results and to improve speed. 
#
# modified on Feb. 14, 2021
#  (1) rename 'n' to 'nSubj' and 'm' to 'nCellPerSubj'
#  (2) calculate cutoff based on data from null distribution
#
# modified on Feb. 10, 2021
#  (1) NBZIMM not in CRAN yet. temporarily remove NBZIMM
#
# created on Feb. 8, 2021
#  (1) calculate power for eQTL via scRNAseq data by using simulation approach
#  (1) consider NBZIMM since it is much faster than GLMMadaptive and glmmTMB
#
########################################################

powerEQTL.scRNAseq.sim=function(
  slope, 
  n, 
  m, 
  power = NULL,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim  =  1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1, # number of computer cores used by 'mclapply'
  n.lower = 2.01,
  n.upper = 1e+4,
  slope.lower =  1e-6,
  slope.upper = log(1.0e+6),
  MAF.lower = 0.05,
  MAF.upper = 0.49
)
{
  if(is.null(MAF)==TRUE && is.null(slope) == FALSE &&
     is.null(n) == FALSE && is.null(power) == FALSE)
  {
    MAF = minMAFEQTL.scRNAseq.sim(
      slope = slope,
      n = n, 
      m = m,
      power = power,
  m.int = m.int, # mean of random intercept
  sigma.int = sigma.int, # SD of the random intercept
  zero.p = zero.p, # probability that an excess zero occurs
  theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  FWER = FWER,
  nTests = nTests,
  nSim = nSim, # number of simulations
  estMethod = estMethod, # parameter estimation method for ZINB
  nCores = nCores, # number of computer cores used by 'mclapply'
  MAF.lower = MAF.lower,
  MAF.upper = MAF.upper
)
    names(MAF) = "MAF"
    return(MAF)
  } else if (is.null(MAF)==FALSE && is.null(slope) == TRUE &&
             is.null(n) == FALSE && is.null(power) == FALSE) {
    slope = minSlopeEQTL.scRNAseq.sim(
      n = n, 
      m = m,
      power = power,
  m.int = m.int, # mean of random intercept
  sigma.int = sigma.int, # SD of the random intercept
  zero.p = zero.p, # probability that an excess zero occurs
  theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = MAF, 
  FWER = FWER,
  nTests = nTests,
  nSim = nSim, # number of simulations
  estMethod = estMethod, # parameter estimation method for ZINB
  nCores = nCores, # number of computer cores used by 'mclapply'
  slope.lower = slope.lower,
  slope.upper = slope.upper
      )

    names(slope) = "slope"
    return(slope)
  } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
             is.null(n) == TRUE && is.null(power) == FALSE) {
    n = ssEQTL.scRNAseq.sim(  
          slope = slope, 
          m = m,
          power = power,
  m.int = m.int, # mean of random intercept
  sigma.int = sigma.int, # SD of the random intercept
  zero.p = zero.p, # probability that an excess zero occurs
  theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = MAF, 
  FWER = FWER,
  nTests = nTests,
  nSim  =  nSim, # number of simulations
  estMethod = estMethod, # parameter estimation method for ZINB
  nCores = nCores, # number of computer cores used by 'mclapply'
  n.lower = n.lower,
  n.upper = n.upper)

    names(n)="n"
    return(n)
  } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
             is.null(n) == FALSE && is.null(power) == TRUE) {
    power = powerEQTL.scRNAseq.sim.default(
      slope = slope, 
      n = n, 
      m = m, 
      m.int = m.int, # mean of random intercept
      sigma.int = sigma.int, # SD of the random intercept
      zero.p = zero.p, # probability that an excess zero occurs
      theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
      # the smaller theta is, the larger variance of NB random variable is
      MAF = MAF, 
      FWER = FWER,
      nTests = nTests,
      nSim = nSim, # number of simulations
      estMethod = estMethod, # parameter estimation method for ZINB
      nCores = nCores # number of computer cores used by 'mclapply'
    )$power  
    names(power) = "power"
    return(power)
  } else {
    stop("One and only one of the 4 parameters (MAF, slope, n, power) can be NULL!\n")
  }
  
}

powerEQTL.scRNAseq.sim.default = function(
  slope, # slope
  n, # number of subjects
  m, # number of cells per subjects 
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim = 1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1 # number of computer cores used by 'mclapply'
)
{
  #nSubj = 102
  #nCellPerSubj = 227868

  nSubj = n
  nCellPerSubj = m
  
  if(estMethod == "GLMMadaptive")  {
    pvalVec = unlist(mclapply(1:nSim, function(i) {
      # simulate data under alternative hypothesis: slope not equal to zero
      simDati = simDat.eQTL.scRNAseq(nSubj=nSubj,
        nCellPerSubj = nCellPerSubj,
        # probability that an excess zero occurs
        zero.p = zero.p,
        m.int = m.int,
        # standard deviation of the random intercept
        sigma.int = sigma.int,
        # slope for genotype in NB distribution: mu = exp(beta0 + beta1*SNP)
        slope = slope,
        # dispersion parameter of NB distribution
        # the smaller theta is, the larger variance of NB random variable is
        theta  = theta,
        MAF = MAF # SNP MAF
      )

      #########################
      # test if slope is significantly different from zero 
      res.try = try(f <- GLMMadaptive::mixed_model(fixed = counts ~ geno, 
                    random = ~ 1 | id, 
      	      data = simDati,
                    zi_fixed = ~ 1, 
	      family = zi.negative.binomial()), silent = TRUE)
      aa = attr(res.try, which="class")
      if(aa == "try-error")
      {
        pval = NA 
      } else {
        pval = summary(f)$coef_table[2, 4]
      }

      return(pval)

    }, mc.cores=nCores))
  } else {
    stop("Currently available choice for 'estMethod' is GLMMadaptive!")
  }

  # power
  power = mean( pvalVec < FWER/nTests, na.rm =TRUE)
  res = list(power = power, pvalVec=pvalVec)

  invisible(res)
}


# cutoff of test statistic is based on null distribution of the test statistic
powerEQTL.scRNAseq.sim.default.exact = function(
  slope, # slope
  n, # number of subjects
  m, # number of cells per subjects 
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim = 1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1 # number of computer cores used by 'mclapply'
)
{
  #nSubj = 102
  #nCellPerSubj = 227868

  nSubj = n
  nCellPerSubj = m
  
  if(estMethod == "GLMMadaptive")  {
    res.sim = mclapply(1:nSim, function(i) {
      # simulate data under alternative hypothesis: slope not equal to zero
      simDati = simDat.eQTL.scRNAseq(nSubj=nSubj,
        nCellPerSubj = nCellPerSubj,
        # probability that an excess zero occurs
        zero.p = zero.p,
        m.int = m.int,
        # standard deviation of the random intercept
        sigma.int = sigma.int,
        # slope for genotype in NB distribution: mu = exp(beta0 + beta1*SNP)
        slope = slope,
        # dispersion parameter of NB distribution
        # the smaller theta is, the larger variance of NB random variable is
        theta  = theta,
        MAF = MAF # SNP MAF
      )

      # simulate data under null hypothesis: slope = 0
      simDat0 = simDat.eQTL.scRNAseq(nSubj=nSubj,
        nCellPerSubj = nCellPerSubj,
        # probability that an excess zero occurs
        zero.p = zero.p,
        m.int = m.int,
        # standard deviation of the random intercept
        sigma.int = sigma.int,
        # slope for genotype in NB distribution: mu = exp(beta0 + beta1*SNP)
        slope = 0,
        # dispersion parameter of NB distribution
        # the smaller theta is, the larger variance of NB random variable is
        theta  = theta,
        MAF = MAF # SNP MAF
      )
 
      
      #########################
      # test if slope is significantly different from zero 
      res.try = try(f <- GLMMadaptive::mixed_model(fixed = counts ~ geno, 
                    random = ~ 1 | id, 
      	      data = simDati,
                    zi_fixed = ~ 1, 
	      family = zi.negative.binomial()), silent = TRUE)
      aa = attr(res.try, which="class")
      if(aa == "try-error")
      {
        stat1 = NA 
      } else {
        stat1 = summary(f)$coef_table[2, 3]
      }

      # under H0
      res.try0 = try(f0 <- GLMMadaptive::mixed_model(fixed = counts ~ geno, 
                    random = ~ 1 | id, 
      	      data = simDat0,
                    zi_fixed = ~ 1, 
	      family = zi.negative.binomial()), silent = TRUE)
      
      aa0 = attr(res.try0, which="class")
      if(aa0 == "try-error")
      {
        stat0 = NA 
      } else {
        stat0 = summary(f0)$coef_table[2, 3]
      }

      res = c(stat1, stat0)
      return(res)

    }, mc.cores=nCores)
  } else {
    stop("Currently available choice for 'estMethod' is GLMMadaptive!")
  }
  
  mat = t(sapply(res.sim, function(x) {x}))
  colnames(mat) = c("Ha", "H0")

  stat1 = mat[,1] # test statistic under Ha
  stat0 = mat[,2] # test statistic under H0
  # cutoff
  alpha2 = FWER/nTests
  cutoff.upp = quantile(stat0, prob= 1 - alpha2/2, na.rm=TRUE)
  cutoff.low = quantile(stat0, prob= alpha2/2, na.rm=TRUE)


  # power
  if(is.na(cutoff.upp)==TRUE || is.na(cutoff.low) == TRUE || any(is.na(stat1) == FALSE) == FALSE)
  {
    power = 0
  } else {
    power = mean(stat1 < cutoff.low | stat1  > cutoff.upp, na.rm=TRUE)
  }

  res = list(power = power, stat0=stat0, stat1=stat1, cutoff.low=cutoff.low, cutoff.upp = cutoff.upp)
  invisible(res)
}

###########################################################

# difference between estimated power and desired power
diffPower4ss.scRNAseq.sim=function( n, 
                                slope, 
                                 m,
                                power,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim  =  1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1 # number of computer cores used by 'mclapply'
)
{
  n = ceiling(n)
  est.power=powerEQTL.scRNAseq.sim.default(
    slope = slope, 
    n = n, 
    m = m, 
  m.int = m.int, # mean of random intercept
      sigma.int = sigma.int, # SD of the random intercept
      zero.p = zero.p, # probability that an excess zero occurs
      theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
      # the smaller theta is, the larger variance of NB random variable is
      MAF = MAF, 
      FWER = FWER,
      nTests = nTests,
      nSim = nSim, # number of simulations
      estMethod = estMethod, # parameter estimation method for ZINB
      nCores = nCores # number of computer cores used by 'mclapply'
   )$power  

  diff=est.power - power
  return(diff)

}


ssEQTL.scRNAseq.sim=function(  slope, 
                           m,
                           power = 0.8,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim  =  1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1, # number of computer cores used by 'mclapply'
  n.lower = 2.01,
  n.upper = 1e+4)
{

  loop = 0
  
  while(loop < 10)
  {
    loop = loop + 1


    res.try = try(res.root <- uniroot(f=diffPower4ss.scRNAseq.sim,
                     interval = c(n.lower, n.upper),
                     slope = slope, 
                     m = m,
                     power = power,
    m.int = m.int, # mean of random intercept
    sigma.int = sigma.int, # SD of the random intercept
    zero.p = zero.p, # probability that an excess zero occurs
    theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
    # the smaller theta is, the larger variance of NB random variable is
    MAF = MAF, 
    FWER = FWER,
    nTests = nTests,
    nSim = nSim, # number of simulations
    estMethod = estMethod, # parameter estimation method for ZINB
    nCores = nCores # number of computer cores used by 'mclapply'
                    ), silent = TRUE)

    aa = attr(res.try, which="class")
    if(is.null(aa) == FALSE)
    {
      res.root = list(root = NA)
    } else {
      break
    }   

  }

  return(ceiling(res.root$root))
}

###########################################
# difference between estimated power and desired power
diffPower4slope.scRNAseq.sim=function(
                                slope,
                                n,
                                m,
                                power,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim = 1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1 # number of computer cores used by 'mclapply'
)
{
  est.power=powerEQTL.scRNAseq.sim.default(
    slope = slope, 
    n = n, 
    m = m, 
  m.int = m.int, # mean of random intercept
      sigma.int = sigma.int, # SD of the random intercept
      zero.p = zero.p, # probability that an excess zero occurs
      theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
      # the smaller theta is, the larger variance of NB random variable is
      MAF = MAF, 
      FWER = FWER,
      nTests = nTests,
      nSim = nSim, # number of simulations
      estMethod = estMethod, # parameter estimation method for ZINB
      nCores = nCores # number of computer cores used by 'mclapply'
   )$power  


  diff=est.power - power
  return(diff)
  
}

minSlopeEQTL.scRNAseq.sim=function(n, 
                           m,
                           power = 0.8,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  MAF = 0.2, 
  FWER = 0.05,
  nTests = 1,
  nSim  =  1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1, # number of computer cores used by 'mclapply'
  slope.lower =  1e-6,
  slope.upper = log(1.0e+6)
)
{
  
  loop = 0
  
  while(loop < 10)
  {
    loop = loop + 1

    res.try = try(res.root <- uniroot(f=diffPower4slope.scRNAseq.sim,
                     interval = c(slope.lower, slope.upper),
                     n = n, 
                     m = m,
                     power = power,
    m.int = m.int, # mean of random intercept
    sigma.int = sigma.int, # SD of the random intercept
    zero.p = zero.p, # probability that an excess zero occurs
    theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
    # the smaller theta is, the larger variance of NB random variable is
    MAF = MAF, 
    FWER = FWER,
    nTests = nTests,
    nSim = nSim, # number of simulations
    estMethod = estMethod, # parameter estimation method for ZINB
    nCores = nCores # number of computer cores used by 'mclapply'
    ), silent = TRUE)

    aa = attr(res.try, which="class")
    if(is.null(aa) == FALSE)
    {
      res.root = list(root = NA)
    } else {
      break
    }   
  }
  
  return(res.root$root)
}

##################################################

# difference between estimated power and desired power
diffPower4MAF.scRNAseq.sim=function(
  MAF,
  slope,
  n,
  m,
  power,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  FWER = 0.05,
  nTests = 1,
  nSim = 1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1 # number of computer cores used by 'mclapply'
)
{
  est.power=powerEQTL.scRNAseq.sim.default(
    slope = slope, 
    n = n, 
    m = m, 
  m.int = m.int, # mean of random intercept
      sigma.int = sigma.int, # SD of the random intercept
      zero.p = zero.p, # probability that an excess zero occurs
      theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
      # the smaller theta is, the larger variance of NB random variable is
      MAF = MAF, 
      FWER = FWER,
      nTests = nTests,
      nSim = nSim, # number of simulations
      estMethod = estMethod, # parameter estimation method for ZINB
      nCores = nCores # number of computer cores used by 'mclapply'
   )$power  

  diff=est.power - power
  return(diff)
  
}

minMAFEQTL.scRNAseq.sim=function(slope, 
                             n, 
                               m,
                               power = 0.8,
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.1, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  # the smaller theta is, the larger variance of NB random variable is
  FWER = 0.05,
  nTests = 1,
  nSim  =  1000, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1, # number of computer cores used by 'mclapply'
  MAF.lower = 0.05,
  MAF.upper = 0.049
)
{

  loop = 0
  
  while(loop < 10)
  {
    loop = loop + 1
    res.try = try(res.root <- uniroot(f=diffPower4MAF.scRNAseq.sim,
                     interval = c(MAF.lower, MAF.upper),
                     slope = slope,
                     n = n, 
                     m = m,
                     power = power,
    m.int = m.int, # mean of random intercept
        sigma.int = sigma.int, # SD of the random intercept
        zero.p = zero.p, # probability that an excess zero occurs
        theta = theta, # dispersion parameter of NB distribution NB(mu, theta)
        # the smaller theta is, the larger variance of NB random variable is
        FWER = FWER,
        nTests = nTests,
        nSim = nSim, # number of simulations
        estMethod = estMethod, # parameter estimation method for ZINB
        nCores = nCores # number of computer cores used by 'mclapply'
    ), silent = TRUE)
    aa = attr(res.try, which="class")
    if(is.null(aa) == FALSE)
    {
      res.root = list(root = NA)
    } else {
      break
    }   
  }
  
  return(res.root$root)
}

