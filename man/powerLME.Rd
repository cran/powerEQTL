\name{powerLME}
\alias{powerLME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Simple Linear Mixed Effects Model
}
\description{
Power calculation for simple linear mixed effects model.
}
\usage{
powerLME(
  delta, 
  n, 
  m, 
  sigma.y, 
  sigma.x, 
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
  \item{sigma.x}{
numeric. standard deviation of the predictor x
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
\deqn{y_{ij} = \beta_{0i} + \beta_1 * x_i + \epsilon_{ij},}
where
\deqn{\beta_{0i} \sim N\left(\beta_0, \sigma^2_{\beta}\right),}
and
\deqn{\epsilon_{ij} \sim N\left(0, \sigma^2_{\epsilon}\right),}
\eqn{i=1,\ldots, n}, \eqn{j=1,\ldots, m}, \eqn{n} is the number of
subjects, \eqn{m} is the number of observations per subject,
\eqn{y_{ij}} is the outcome value for the \eqn{j}-th observation
of the \eqn{i}-th subject, \eqn{x_i} is the predictor value for the
\eqn{i}-th subject. 

We would like to test the following hypotheses:
\deqn{H_0: \beta_1=0,}
and
\deqn{H_a: \beta_1 = \delta,}
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
\eqn{\sigma_y=\sqrt{\sigma^2_{\beta}+\sigma^2_{\epsilon}}},
\eqn{\hat{\sigma}_x=\sqrt{\sum_{i=1}^n\left(x_i-\bar{x}\right)^2/(n-1)}},
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
  n = 102
  m = 227868
  power = powerLME(
    delta = 0.6, 
    n = n, 
    m = m,
    sigma.y = 0.29, 
    sigma.x = 0.308, 
    rho = 0.8, 
    alpha = 0.05,
    nTests = 1e+6)

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
