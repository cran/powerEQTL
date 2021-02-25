\name{powerLME}
\alias{powerLME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Simple Linear Mixed Effects Model
}
\description{
Power calculation for simple linear mixed effects model.
This function can be used to calculate one of the 3 parameters (power, sample size, and minimum detectable slope) by setting the corresponding parameter
as NULL and providing values for the other 2 parameters.
}
\usage{
powerLME(
  slope, 
  n, 
  m, 
  sigma.y, 
  sigma.x,
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
integer. Number of observations per subject.
}
  \item{sigma.y}{
numeric. Standard deviation of the outcome y.
}
  \item{sigma.x}{
numeric. Standard deviation of the predictor x.
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
integer. Number of tests (e.g., number of genes in differential expression analysis based on scRNAseq to compare gene expression between diseased subjects and healthy subjects).
}
  \item{n.lower}{
numeric. Lower bound of the sample size. Only used when calculating sample size.
} \item{n.upper}{
numeric. Upper bound of the sample size. Only used when calculating sample size.
}

}
\details{
We assume the following simple linear mixed effects model to characterize the association between
the predictor x and the outcome y:
\deqn{y_{ij} = \beta_{0i} + \beta_1 * x_i + \epsilon_{ij},}
where
\deqn{\beta_{0i} \sim N\left(\beta_0, \sigma^2_{\beta}\right),}
and
\deqn{\epsilon_{ij} \sim N\left(0, \sigma^2\right),}
\eqn{i=1,\ldots, n}, \eqn{j=1,\ldots, m}, \eqn{n} is the number of
subjects, \eqn{m} is the number of observations per subject,
\eqn{y_{ij}} is the outcome value for the \eqn{j}-th observation
of the \eqn{i}-th subject, \eqn{x_i} is the predictor value for the
\eqn{i}-th subject. For example, \eqn{x_i} is the binary variable indicating if the \eqn{i}-th subject is a diseased subject or not.

We would like to test the following hypotheses:
\deqn{H_0: \beta_1=0,}
and
\deqn{H_1: \beta_1 = \delta,}
where \eqn{\delta\neq 0}.

We can derive the power calculation formula is
\deqn{power=1- \Phi\left(z_{\alpha^{*}/2}-a\times b\right)
+\Phi\left(-z_{\alpha^{*}/2} - a\times b\right),}
where
\deqn{a= 
\frac{\hat{\sigma}_x
}{\sigma_y}}
and
\deqn{
b=\frac{\delta\sqrt{m(n-1)}}{\sqrt{1+(m-1)\rho}}
}
and \eqn{z_{\alpha^{*}/2}} is the upper \eqn{100\alpha^{*}/2}
percentile of the standard normal distribution,
\eqn{\alpha^{*}=\alpha/nTests}, nTests is the number of 
tests,
\eqn{\sigma_y=\sqrt{\sigma^2_{\beta}+\sigma^2}},
\eqn{\hat{\sigma}_x=\sqrt{\sum_{i=1}^n\left(x_i-\bar{x}\right)^2/(n-1)}},
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
  n = 102
  m = 227868
  
  # calculate power
  power = powerLME(
    slope = 0.6, 
    n = n, 
    m = m,
    sigma.y = 0.29, 
    sigma.x = 0.308,
    power = NULL,
    rho = 0.8, 
    FWER = 0.05,
    nTests = 1e+6)

  print(power)
  
  # calculate sample size
  n = powerLME(
    slope = 0.6, 
    n = NULL, 
    m = m,
    sigma.y = 0.29, 
    sigma.x = 0.308,
    power = 0.9562555,
    rho = 0.8, 
    FWER = 0.05,
    nTests = 1e+6)

  print(n)
  
  # calculate slope
  slope = powerLME(
    slope = NULL, 
    n = n, 
    m = m,
    sigma.y = 0.29, 
    sigma.x = 0.308,
    power = 0.9562555,
    rho = 0.8, 
    FWER = 0.05,
    nTests = 1e+6)

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
