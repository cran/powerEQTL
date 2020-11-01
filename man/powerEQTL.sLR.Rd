\name{powerEQTL.SLR}
\alias{powerEQTL.SLR}

\title{
Power Calculation for EQTL Analysis Based on Simple Linear Regression
}
\description{
Power calculation for eQTL analysis that tests if a SNP is associated to a gene probe by using simple linear regression. This function can be used to calculate one of the 4 parameters (power, sample size, minimum detectable slope, and minimum allowable MAF) by setting the corresponding parameter
as NULL and providing values for the other 3 parameters.
}
\usage{
powerEQTL.SLR(
  MAF,
  slope = 0.13,
  n = 200,
  power = NULL,
  sigma.y = 0.13,
  FWER = 0.05,
  nTests = 2e+05,
  n.lower = 2.01,
  n.upper = 1e+30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAF}{
numeric. Minor allele frequency.
}
  \item{slope}{
numeric. Slope of the simple linear regression.
}
  \item{n}{
integer. Number of subjects.
}
  \item{power}{
numeric. Power for testing if the slope is equal to zero.
  }
  \item{sigma.y}{
numeric. Standard deviation of the outcome \eqn{y_i} in simple linear regression.
}
  \item{FWER}{
numeric. Family-wise Type I error rate.
}
  \item{nTests}{
integer. Number of tests in eQTL analysis.
} 
  \item{n.lower}{
numeric. Lower bound of the sample size. Only used when calculating sample size.
} \item{n.upper}{
numeric. Upper bound of the sample size. Only used when calculating sample size.
}
}
\details{
To test if a SNP is associated with a gene probe, we use the simple linear regression
\deqn{y_i = \beta_0+\beta_1 x_i + \epsilon_i,}{y_i = beta_0 + beta_1 * x_i + epsilon_i,}
where \eqn{y_i}{y_i} is the gene expression level of the \eqn{i}{i}-th subject,
\eqn{x_i}{x_i} is the genotype of the \eqn{i}{i}-th subject, and
\eqn{\epsilon_i} is the random error term with mean zero and standard deviation \eqn{\sigma}. Additive coding for genotype is used. That is,
\eqn{x_i=0} indicates wildtype homozygotes;
\eqn{x_i=1} indicates heterozygotes; and \eqn{x_i=2} indicates mutation heterozygotes.

To test if the SNP is associated with the gene probe, we test the null hypothesis \eqn{H_0: \beta_1=0}{H_0: beta_1 = 0} versus the alternative hypotehsis \eqn{H_1: \beta_1 = \delta}, where \eqn{\delta\neq 0}.


Denote \eqn{\theta}{theta} as the minor allele frequency (MAF) of the SNP. Under Hardy-Weinberg equilibrium, we can calculate the variance of genotype of the SNP: 
\eqn{\sigma^2_x=2 \theta (1-\theta)}{sigma^2_x = 2 * theta * (1 - theta)}, 
where \eqn{\sigma^2_x}{sigma^2_x} is the variance of the
predictor (i.e. the SNP) \eqn{x_i}{x_i}.

The exact power calculation formula can be derived as
\deqn{1-T_{n-2, \lambda}(t_{n-2}(\alpha/2)) + T_{n-2, \lambda}(-t_{n-2}(\alpha/2)),}
where \eqn{T_{n-2, \lambda}(a)} is the value at \eqn{a} of cumulative distribution function of non-central t distribution with \eqn{n-2} degrees of freedom
and non-centrality parameter \eqn{\lambda=\delta/\sqrt{\sigma^2/[(n-1)\tilde{\sigma}^2_{x}]}}. And \eqn{\tilde{\sigma}^2_{x}=\sum_{i=1}^n(x_i - \bar{x})^2/(n-1)}. 

Dupont and Plummer (1998) mentioned the following relationship:
\deqn{\sigma^2 = \sigma^2_y - \beta_1^2 \sigma^2_x.}
So we can plug in the above equation to the power calculation formula.

Under Hardy-Weinberg equilibrium, we have \eqn{\sigma_x^2=2\theta(1-\theta)},
where \eqn{\theta} is the minor allele frequency (MAF).

Hence, the non-centrality parameter can be rewritten as
\deqn{\lambda=\frac{\delta}{\sqrt{
\left(\sigma_y^2 - \delta^2 2\left(1-\hat{\theta}\right)\hat{\theta}\right)/
\left[(n-1)2\left(1-\hat{\theta}\right)\hat{\theta}\right]
}}}
}
\value{
power if the input parameter \code{power = NULL}.

sample size if the input parameter \code{n = NULL};

minimum detectable slope if the input parameter \code{slope = NULL};

minimum allowable MAF if the input parameter \code{MAF = NULL}.

}
\references{
  Dupont, W.D. and Plummer, W.D.. 
  Power and Sample Size Calculations for Studies Involving Linear Regression.
  Controlled Clinical Trials. 1998;19:589-601.

Dong X, Li X, Chang T-W, Weiss ST, and Qiu W. 
powerEQTL: an R package and R shiny application for calculating sample size and power of bulk tissue and single-cell eQTL analysis.
manuscript. (2020)
}
\author{
Xianjun Dong <XDONG@rics.bwh.harvard.edu>,
Tzuu-Wang Chang <Chang.Tzuu-Wang@mgh.harvard.edu>,
Scott T. Weiss <restw@channing.harvard.edu>,
Weiliang Qiu <weiliang.qiu@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{
# calculate power
powerEQTL.SLR(
  MAF = 0.1,
  slope = 0.13,
  n = 179,
  power = NULL,
  sigma.y = 0.13,
  FWER = 0.05,
  nTests = 2e+05)
  
# calculate sample size
powerEQTL.SLR(
  MAF = 0.1,
  slope = 0.13,
  n = NULL,
  power = 0.8,
  sigma.y = 0.13,
  FWER = 0.05,
  nTests = 2e+05)
  
# calculate minimum detectable slope
powerEQTL.SLR(
  MAF = 0.1,
  slope = NULL,
  n = 179,
  power = 0.8,
  sigma.y = 0.13,
  FWER = 0.05,
  nTests = 2e+05)  
  
# calculate minimum allowable MAF
powerEQTL.SLR(
  MAF = NULL,
  slope = 0.13,
  n = 179,
  power = 0.8,
  sigma.y = 0.13,
  FWER = 0.05,
  nTests = 2e+05) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

