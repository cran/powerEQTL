\name{simDat.eQTL.scRNAseq}
\alias{simDat.eQTL.scRNAseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Generate Gene Expression Levels Of One Gene And Genotypes Of One SNP For Subjects With Multiple Cells Based On ZINB Mixed Effects Regression Model
}
\description{
	Generate gene expression levels of one gene and genotypes of one SNP for subjects with multiple cells based on ZINB mixed effects regression model.
}
\usage{
simDat.eQTL.scRNAseq(nSubj = 50, 
		     nCellPerSubj = 100, 
		     zero.p = 0.01, 
		     m.int = 0, 
		     sigma.int = 1, 
		     slope = 1, 
		     theta = 1, 
		     MAF = 0.45)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSubj}{
	  integer. Total number of subjects.
}
  \item{nCellPerSubj}{
	  integer. Number of cells per subject.
}
  \item{zero.p}{
	  numeric. Probability that an excess zero occurs.
}
  \item{m.int}{
	  numeric. Mean of random intercept (see details).
}
  \item{sigma.int}{
	  numeric. Standard deviation of random intercept (see details).
}
  \item{slope}{
	  numeric. Slope (see details).
}
  \item{theta}{
	  numeric. dispersion parameter of negative binomial distribution.
  The smaller \code{theta} is, the larger variance of NB random variable is.
}
  \item{MAF}{
	  numeric. Minor allele frequency of the SNP.
}
}
\details{
	This function simulates gene expression levels of one gene and genotypes of one SNP for subjects with multiple cells based on zero-inflated negative binomial (ZINB) regression model with only one 
covariate: genotype.  
That is, the read counts of a gene follows a mixture of 2-component distributions. One component takes only one value: zero. The other component is negative binomial distribution, which takes non-negative values 0, 1, 2, .... The log mean of the negative binomial distribution is linear function of the genotype.

	Denote \eqn{Y_{ij}} as the read counts for the \eqn{j}-th cell of
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
	probabilities of the 3 genotypes \eqn{(0, 1, 2)} are \eqn{(1-MAF)^2,
		2 MAF (1-MAF), MAF^2}, respectively.

For simplicity, we assume that excess zeros are caused by technical issues, hence are not related to genotypes.

}
\value{
	A data frame with 3 columns:
\item{id }{subject id}
\item{geno }{additive-coded genotype of the SNP}
\item{counts }{gene expression of the gene}
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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
frame = simDat.eQTL.scRNAseq(nSubj = 5, 
		     nCellPerSubj = 3, 
		     zero.p = 0.01, 
		     m.int = 0, 
		     sigma.int = 1, 
		     slope = 1, 
		     theta = 1, 
		     MAF = 0.45)
print(dim(frame))
print(frame[1:10,])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ method }
