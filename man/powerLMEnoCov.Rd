\name{powerLMEnoCov}
\alias{powerLMEnoCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Simple Linear Mixed Effects Model Without Covariate
}
\description{
Power calculation for simple linear mixed effects model without covariate.
}
\usage{
powerLMEnoCov(
  delta, 
  n, 
  m, 
  sigma.y, 
  rho = 0.8, 
  alpha = 0.05,
  nTests = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
numeric. slope under alternative hypothesis, i.e., \eqn{\beta_1=\delta} under \eqn{H_a: \beta_1\neq 0}.
}
  \item{n}{
integer. number of subjects
}
  \item{m}{
integer. number of observations per subject
}
  \item{sigma.y}{
numeric. conditional standard deviation of the outcome y given the predictor x
}
  \item{rho}{
numeric. intra-class correlation (i.e., correlation between 
\eqn{y_{ij}} and \eqn{y_{ik}} for the \eqn{j}-th and \eqn{k}-th observations of the \eqn{i}-th subject)
}
  \item{alpha}{
numeric. Type I error rate for one test 
}
  \item{nTests}{
integer. Number of tests
}

}
\details{
We assume the following simple linear mixed effects model to characterize the association between
the predictor x and the outcome y:
\deqn{y_{ij} = \beta_{0i} + \epsilon_{ij},}
where
\deqn{\beta_{0i} \sim N\left(\beta_0, \sigma^2_{\beta}\right),}
and
\deqn{\epsilon_{ij} \sim N\left(0, \sigma^2_{\epsilon}\right),}
\eqn{i=1,\ldots, n}, \eqn{j=1,\ldots, m}, \eqn{n} is the number of
subjects, \eqn{m} is the number of observations per subject,
\eqn{y_{ij}} is the outcome value for the \eqn{j}-th observation
of the \eqn{i}-th subject. 

We would like to test the following hypotheses:
\deqn{H_0: \beta_0=0,}
and
\deqn{H_a: \beta_0 = \delta,}
where \eqn{\delta\neq 0}.

We can derive the power calculation formula is
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
\eqn{\sigma_y=\sqrt{\sigma^2_{\beta}+\sigma^2_{\epsilon}}}
and \eqn{\rho=\sigma^2_{\beta}/\left(\sigma^2_{\beta}+\sigma^2_{\epsilon}\right)} is the intra-class correlation.
}
\value{
the power to be calculated.
}
\references{
  Dong X and Qiu W. 
  Power Calculation for Association Between Genotype and Gene Expression Based on Single Cell RNAseq Data.
  manuscript. (2020)
}
\author{
Xianjun Dong <XDONG@rics.bwh.harvard.edu>,
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
  power = powerLMEnoCov(
    delta = 1.3*sigma.y, 
    n = n, 
    m = m,
    sigma.y = sigma.y, 
    rho = 0.8, 
    alpha = 0.05,
    nTests = 20345)

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
