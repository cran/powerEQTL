\name{powerEQTL.scRNAseq}
\alias{powerEQTL.scRNAseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Association Between Genotype and Gene Expression Based on Single Cell RNAseq Data
}
\description{
Power calculation for association between genotype and gene expression based on single cell RNAseq data.
}
\usage{
powerEQTL.scRNAseq(
  delta, 
  n, 
  m, 
  sigma.y, 
  theta = 0.2, 
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
integer. number of cells per subject
}
  \item{sigma.y}{
numeric. conditional standard deviation of the gene expression given SNP
}
  \item{theta}{
numeric. minor allele frequency (between 0 and 0.5)
}
  \item{rho}{
numeric. intra-class correlation (i.e., correlation between 
\eqn{y_{ij}} and \eqn{y_{ik}} for the \eqn{j}-th  and \eqn{k}-th cells of the \eqn{i}-th subject)
}
  \item{alpha}{
numeric. Type I error rate for one pair (SNP, gene)
}
  \item{nTests}{
integer. Number of tests, i.e., number of (SNP, gene) pairs, in eQTL analysis
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
\deqn{\epsilon_{ij} \sim N\left(0, \sigma^2_{\epsilon}\right),}
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
\deqn{H_a: \beta_1 = \delta,}
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
\eqn{\alpha^{*}=\alpha/nTests}, nTests is the number of
(SNP, gene) pairs,
\eqn{\sigma_y=\sqrt{\sigma^2_{\beta}+\sigma^2_{\epsilon}}},
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
  power = powerEQTL.scRNAseq(
    delta = 0.6, 
    n = n, 
    m = m,
    sigma.y = 0.29, 
    theta = 0.05, 
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
