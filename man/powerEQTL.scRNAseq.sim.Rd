\name{powerEQTL.scRNAseq.sim}
\alias{powerEQTL.scRNAseq.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Association Between Genotype and Gene Expression Based on Single Cell RNAseq Data via Simulation from ZINB Mixed Effects Regression Model
}
\description{
Power calculation for association between genotype and gene expression based on single cell RNAseq data via simulation from ZINB mixed effects regression model. This function can be used to calculate one of the 4 parameters (power, sample size, minimum detectable slope, and minimum allowable MAF) by setting the corresponding parameter
as NULL and providing values for the other 3 parameters.
}
\usage{
powerEQTL.scRNAseq.sim(slope, 
		       n, 
		       m, 
                       power = NULL,
		       m.int = -1, 
		       sigma.int = 1, 
		       zero.p = 0.1, 
		       theta = 1, 
		       MAF = 0.2, 
		       FWER = 0.05, 
		       nTests = 1, 
		       nSim = 1000, 
		       estMethod = "GLMMadaptive", 
		       nCores = 1,
                       n.lower = 2.01,
                       n.upper = 1e+4,
                       slope.lower =  1e-6,
                       slope.upper = log(1.0e+6),
                       MAF.lower = 0.05,
                       MAF.upper = 0.49
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slope}{
	  numeric. Slope (see details).
}
  \item{n}{
	  integer. Total number of subjects.
}
  \item{m}{
	  integer. Number of cells per subject.
}
  \item{power}{
numeric. Power for testing if the slope is equal to zero.
}
  \item{m.int}{
	  numeric. Mean of random intercept (see details).
}
  \item{sigma.int}{
	  numeric. Standard deviation of random intercept (see details).
}

  \item{zero.p}{
	  numeric. Probability that an excess zero occurs.
}
  \item{theta}{
	  numeric. Dispersion parameter of negative binomial distribution.
  The smaller \code{theta} is, the larger variance of NB random variable is.
}
  \item{MAF}{
	  numeric. Minor allele frequency of the SNP.
}
  \item{FWER}{
numeric. Family-wise Type I error rate.
}
  \item{nTests}{
	  integer. Number of tests (i.e., number of all (SNP, gene) pairs) in eQTL analysis.
}
  \item{nSim}{
	  integer. Number of simulated datasets to be generated.
}
  \item{estMethod}{
	  character. Indicates which method would be used to fit zero inflated 
	  negative binomial mixed effects model.
	  The possible choices are \dQuote{GLMMadaptive}
	  and \dQuote{glmmTMB}.
}
  \item{nCores}{
	  integer. Number of computer cores used by \code{mclapply} for parallel computing. For Windows, \code{nCores=1}.
}
  \item{n.lower}{
numeric. Lower bound of the total number of subjects. Only used when calculating total number of subjects.
} \item{n.upper}{
numeric. Upper bound of the total number of subjects. Only used when calculating total number of subjects.
}
  \item{slope.lower}{
numeric. Lower bound of the slope. Only used when calculating minimum slope.
} \item{slope.upper}{
numeric. Upper bound of the slope Only used when calculating minimum slope.
}
  \item{MAF.lower}{
numeric. Lower bound of the MAF. Only used when calculating minimum MAF.
} \item{MAF.upper}{
numeric. Upper bound of the MAF Only used when calculating minimum MAF.
}

}
\details{
	This function calculates the power for testing if genotypes of a SNP is associated with the expression of a gene via \code{nSim} datasets generated from zero-inflated negative binomial (ZINB) regression model.

	Each dataset is generated from zero-inflated negative binomial mixed effects regression model with only one covariate: genotype. That is, the read counts of a gene follows a mixture of 2-component distributions. One component takes only one value: zero. The other component is negative binomial distribution, which takes non-negative values 0, 1, 2, .... The log mean of the negative binomial distribution is linear function of the genotype.

	For each dataset, the p-value for testing if the slope for genotype is equal to zero will be calculated.

	The proportion of p-values \eqn{< \alpha} is the estimated power, where \eqn{\alpha = FWER/nTests}.
	
	Each simulated dataset contains gene expression levels of one gene and genotypes of one SNP for subjects with multiple cells.
	The gene expression levels (read counts) follow zero-inflated negative binomial distribution. Denote \eqn{Y_{ij}} as the read counts for the \eqn{j}-th cell of
	the \eqn{i}-th subject, \eqn{i=1,\ldots, n}, \eqn{j=1,\ldots, m},
	\eqn{n} is the number of subjects, and \eqn{m} is the number of cells per subject.
	Denote \eqn{p} as the probability that \eqn{Y_{ij}=0} is an excess zero.
	With probability \eqn{1-p}, \eqn{Y_{ij}} follows a negative binomial distribution \eqn{NB(\mu, \theta)}, where \eqn{\mu} is the mean (i.e., \eqn{\mu=E(Y_{ij})}) and \eqn{\theta} is the dispersion parameter.
The variance of the NB distribution is \eqn{\mu + \mu^2/\theta}.
	The relationship between gene expression and genotype for the \eqn{i}-th subject is characterized by the equation 
	\deqn{\mu_i = \exp(\beta_{0i} + \beta_{1} x_i),}
	where \eqn{\beta_{0i}} is the random intercept following a normal
	distribution \eqn{N(\beta_0, \sigma^2)} to account for within-subject correlation of gene expression, \eqn{\beta_0} is the mean of the random intercept, \eqn{\sigma} is the standard deviation of the random intercept, \eqn{\beta_1} is the slope, and \eqn{x_i} is the additive-coded genotype for the SNP with minor allele frequency \eqn{MAF}.
	
	We assume that the SNP satisfies the Hardy-Weinberg Equilibrium. That is, the
	probabilities of the 3 genotypes \eqn{(0, 1, 2)} are \eqn{(1-MAF)^2, 2 MAF (1-MAF), MAF^2}, respectively.


For simplicity, we assume that excess zeros are caused by technical issues, hence are not related to genotypes.

}
\value{
power if the input parameter \code{power = NULL}.

sample size (total number of subjects) if the input parameter \code{n = NULL};

minimum detectable slope if the input parameter \code{slope = NULL};

minimum allowable MAF if the input parameter \code{MAF = NULL}.

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Xianjun Dong <XDONG@rics.bwh.harvard.edu>,
Xiaoqi Li<xli85@bwh.harvard.edu>,
Tzuu-Wang Chang <Chang.Tzuu-Wang@mgh.harvard.edu>,
Scott T. Weiss <restw@channing.harvard.edu>,
Weiliang Qiu <weiliang.qiu@gmail.com>
}
\note{
The speed of simulation approach is slow. It is recommended to run the function \code{powerEQTL.scRNAseq.sim} in parallel computing environment (e.g., Unix/Linux) and set \code{nCores > 1}.
Also, it is important to set appropriate ranges to search sample size, minimum detectable slope, or minimum allowable MAF. If the function (e.g., glmmTMB) to fit data returns an error message 
for a simulated dataset, then the function \code{powerEQTL.scRNAseq.sim} will try a few more runs with different simulated datasets. If all the runs failed, zero (for power estimation) or \code{NA} or \code{NaN} (for estimating slope, sample size, or MAF) would be output.  It is recommended to set \code{nSim >= 1000} to get stable results.
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

nSubj = 10
nCellPerSubj = 10

# calculate power
power = powerEQTL.scRNAseq.sim(
  slope = 1.62, # slope
  n = nSubj, # total number of subjects
  m = nCellPerSubj, # number of cells per subject 
  power = NULL, # power to be estimated
  m.int = -1, # mean of random intercept
  sigma.int = 1, # SD of the random intercept
  zero.p = 0.01, # probability that an excess zero occurs
  theta = 1, # dispersion parameter of NB distribution NB(mu, theta)
  MAF = 0.45, 
  FWER = 0.05,
  nTests = 1,
  nSim  =  5, # number of simulations
  estMethod = "GLMMadaptive", # parameter estimation method for ZINB
  nCores = 1 # number of computer cores used by 'mclapply'
)

print(power)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ method }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
