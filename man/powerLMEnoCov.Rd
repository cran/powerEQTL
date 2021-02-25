\name{powerLMEnoCov}
\alias{powerLMEnoCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Simple Linear Mixed Effects Model Without Covariate
}
\description{
Power calculation for simple linear mixed effects model without covariate.
This function can be used to calculate one of the 3 parameters (power, sample size, and minimum detectable slope) by setting the corresponding parameter
as NULL and providing values for the other 2 parameters.
}
\usage{
powerLMEnoCov(
  slope, 
  n, 
  m, 
  sigma.y,
  power = NULL,
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
integer. Number of subjects.
}
  \item{m}{
integer. Number of pairs of replicates per subject.
}
  \item{sigma.y}{
numeric. Standard deviation of the outcome y.
}
  \item{power}{
  numeric. Desired power.
  }
  \item{rho}{
numeric. Intra-class correlation (i.e., correlation between 
\eqn{y_{ij}} and \eqn{y_{ik}} for the \eqn{j}-th and \eqn{k}-th observations of the \eqn{i}-th subject).
}
  \item{FWER}{
numeric. Family-wise Type I error rate. 
}
  \item{nTests}{
integer. Number of tests (e.g., number of genes in differential expression analysis based on scRNAseq to compare gene expression before and after treatment).
}
  \item{n.lower}{
numeric. Lower bound of the sample size. Only used when calculating sample size.
} \item{n.upper}{
numeric. Upper bound of the sample size. Only used when calculating sample size.
}

}
\details{
In an experiment, there are \eqn{n} samples. For each sample, we get 
\eqn{m} pairs of replicates. For each pair, one replicate will receive no treatment. The other
replicate will receive treatment. The outcome is the expression of a gene. The overall goal of the experiment is to check if the treatment affects gene expression level or not. Or equivalently, the overall goal of the experiment is to test if the mean within-pair difference of gene expression is equal to zero or not. In the design stage, we would like to calculate the power/sample size of the experiment for testing if the mean within-pair difference of gene expression is equal to zero or not. 


We assume the following linear mixed effects model to characterize the 
relationship between the within-pair difference of gene expression \eqn{y_{ij}} and
the mean of the within-pair difference \eqn{\beta_{0i}}:
\deqn{y_{ij} = \beta_{0i} + \epsilon_{ij},}
where
\deqn{\beta_{0i} \sim N\left(\beta_0, \sigma^2_{\beta}\right),}
and
\deqn{\epsilon_{ij} \sim N\left(0, \sigma^2\right),}
\eqn{i=1,\ldots, n}, \eqn{j=1,\ldots, m}, \eqn{n} is the number of
subjects, \eqn{m} is the number of pairs of replicates per subject,
\eqn{y_{ij}} is the within-pair difference of outcome value for the \eqn{j}-th pair of the \eqn{i}-th subject. 

We would like to test the following hypotheses:
\deqn{H_0: \beta_0=0,}
and
\deqn{H_1: \beta_0 = \delta,}
where \eqn{\delta\neq 0}. If we reject the null hypothesis \eqn{H_0}
based on a sample, we then get evidence that the treatment affects
the gene expression level.

We can derive the power calculation formula:
\deqn{power=1- \Phi\left(z_{\alpha^{*}/2}-a\times b\right)
+\Phi\left(-z_{\alpha^{*}/2} - a\times b\right),}
where
\deqn{a= 
\frac{1
}{\sigma_y}}
and
\deqn{
b=\frac{\delta\sqrt{mn}}{\sqrt{1+(m-1)\rho}}
}
and \eqn{z_{\alpha^{*}/2}} is the upper \eqn{100\alpha^{*}/2}
percentile of the standard normal distribution,
\eqn{\alpha^{*}=\alpha/nTests}, nTests is the number of 
tests,
\eqn{\sigma_y=\sqrt{\sigma^2_{\beta}+\sigma^2}}
and \eqn{\rho=\sigma^2_{\beta}/\left(\sigma^2_{\beta}+\sigma^2\right)} is the intra-class correlation.
}
\value{
power if the input parameter \code{power = NULL}.

sample size if the input parameter \code{n = NULL};

minimum detectable slope if the input parameter \code{slope = NULL}.
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
  n = 17
  m = 5
  sigma.y = 0.68
  slope = 1.3*sigma.y
  print(slope)
  
  # estimate power
  power = powerLMEnoCov(
    slope = slope, 
    n = n, 
    m = m,
    sigma.y = sigma.y,
    power = NULL,
    rho = 0.8, 
    FWER = 0.05,
    nTests = 20345)

  print(power)

  # estimate sample size
  n = powerLMEnoCov(
    slope = slope, 
    n = NULL, 
    m = m,
    sigma.y = sigma.y,
    power = 0.8721607,
    rho = 0.8, 
    FWER = 0.05,
    nTests = 20345)

  print(n)
  
  # estimate slope
  slope = powerLMEnoCov(
    slope = NULL, 
    n = n, 
    m = m,
    sigma.y = sigma.y,
    power = 0.8721607,
    rho = 0.8, 
    FWER = 0.05,
    nTests = 20345)

  print(slope)
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
