\name{powerEQTL.scRNAseq}
\alias{powerEQTL.scRNAseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Association Between Genotype and Gene Expression Based on Single Cell RNAseq Data
}
\description{
Power calculation for association between genotype and gene expression based on single cell RNAseq data. This function can be used to calculate one of the 4 parameters (power, sample size, minimum detectable slope, and minimum allowable MAF) by setting the corresponding parameter
as NULL and providing values for the other 3 parameters.
}
\usage{
powerEQTL.scRNAseq(
  slope, 
  n, 
  m, 
  power = NULL,
  sigma.y, 
  MAF = 0.2, 
  rho = 0.8, 
  FWER = 0.05,
  nTests = 1,
  n.lower = 2.01,
  n.upper = 1e+30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slope}{
numeric. Slope under alternative hypothesis.
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
  \item{sigma.y}{
numeric. Standard deviation of the gene expression.
}
  \item{MAF}{
numeric. Minor allele frequency (between 0 and 0.5).
}
  \item{rho}{
numeric. Intra-class correlation (i.e., correlation between 
\eqn{y_{ij}} and \eqn{y_{ik}} for the \eqn{j}-th  and \eqn{k}-th cells of the \eqn{i}-th subject).
}
  \item{FWER}{
numeric. Family-wise Type I error rate for one pair (SNP, gene).
}
  \item{nTests}{
integer. Number of tests (i.e., number of all (SNP, gene) pairs) in eQTL analysis.
}
  \item{n.lower}{
numeric. Lower bound of the total number of subjects. Only used when calculating total number of subjects.
} \item{n.upper}{
numeric. Upper bound of the total number of subjects. Only used when calculating total number of subjects.
}

}
\details{
We assume the following simple linear mixed effects model for each
(SNP, gene) pair to characterize the association between
genotype and gene expression:
\deqn{y_{ij} = \beta_{0i} + \beta_1 * x_i + \epsilon_{ij},}
where 
\deqn{\beta_{0i} \sim N\left(\beta_0, \sigma^2_{\beta}\right),}
and
\deqn{\epsilon_{ij} \sim N\left(0, \sigma^2\right),}
\eqn{i=1,\ldots, n}, \eqn{j=1,\ldots, m}, \eqn{n} is the number of
subjects, \eqn{m} is the number of cells per subject,
\eqn{y_{ij}} is the gene expression level for the \eqn{j}-th cell
of the \eqn{i}-th subject, \eqn{x_i} is the genotype for the
\eqn{i}-th subject using additive coding. That is, \eqn{x_i=0} indicates the \eqn{i}-th subject is a wildtype homozygote,
\eqn{x_i=1} indicates the \eqn{i}-th subject is a heterozygote,
and \eqn{x_i=2} indicates the \eqn{i}-th subject is a mutation
homozygote.

We would like to test the following hypotheses:
\deqn{H_0: \beta_1=0,}
and
\deqn{H_1: \beta_1 = \delta,}
where \eqn{\delta\neq 0}.

For a given SNP, we assume Hardy-Weinberg Equilibrium and denote the minor allele
frequency of the SNP as \eqn{\theta}.

We can derive the power calculation formula is
\deqn{power=1- \Phi\left(z_{\alpha^{*}/2}-a\times b\right)
+\Phi\left(-z_{\alpha^{*}/2} - a\times b\right),}
where
\deqn{a=
\frac{\sqrt{2\theta\left(1-\theta\right)}}{\sigma_y}}
and
\deqn{
b=\frac{\delta\sqrt{m(n-1)}}{\sqrt{1+(m-1)\rho}}
}
and \eqn{z_{\alpha^{*}/2}} is the upper \eqn{100\alpha^{*}/2}
percentile of the standard normal distribution,
\eqn{\alpha^{*}=FWER/nTests}, nTests is the number of
(SNP, gene) pairs,
\eqn{\sigma_y=\sqrt{\sigma^2_{\beta}+\sigma^2}},
and \eqn{\rho=\sigma^2_{\beta}/\left(\sigma^2_{\beta}+\sigma^2\right)} is the intra-class correlation.
}
\value{
power if the input parameter \code{power = NULL}.

sample size (total number of subjects) if the input parameter \code{n = NULL};

minimum detectable slope if the input parameter \code{slope = NULL};

minimum allowable MAF if the input parameter \code{MAF = NULL}.
}
\references{
Dong X, Li X, Chang T-W, Weiss ST, and Qiu W. 
powerEQTL: an R package and R shiny application for calculating sample size and power of bulk tissue and single-cell eQTL analysis.
manuscript. (2020)
}
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
  n = 102
  m = 227868
  
  # calculate power
  power = powerEQTL.scRNAseq(
    slope = 0.6, 
    n = n, 
    m = m,
    power = NULL,
    sigma.y = 0.29, 
    MAF = 0.05, 
    rho = 0.8, 
    nTests = 1e+6)

  print(power)
  
  # calculate sample size (total number of subjects)
  n = powerEQTL.scRNAseq(
    slope = 0.6, 
    n = NULL, 
    m = m,
    power = 0.9567288,
    sigma.y = 0.29, 
    MAF = 0.05, 
    rho = 0.8, 
    nTests = 1e+6)

  print(n)

  # calculate slope
  slope = powerEQTL.scRNAseq(
    slope = NULL, 
    n = n, 
    m = m,
    power = 0.9567288,
    sigma.y = 0.29, 
    MAF = 0.05, 
    rho = 0.8, 
    nTests = 1e+6)

  print(slope)
  
  # calculate MAF
  MAF = powerEQTL.scRNAseq(
    slope = 0.6, 
    n = n, 
    m = m,
    power = 0.9567288,
    sigma.y = 0.29, 
    MAF = NULL, 
    rho = 0.8, 
    nTests = 1e+6)
  print(MAF)

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
