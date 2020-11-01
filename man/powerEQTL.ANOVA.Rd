\name{powerEQTL.ANOVA}
\alias{powerEQTL.ANOVA}
\title{
Power Calculation for EQTL Analysis Based on Un-Balanced One-Way ANOVA
}
\description{
Power calculation for eQTL analysis that tests if a SNP is associated to a gene probe by using un-balanced one-way ANOVA. This function can be used to calculate one of the 3 parameters (power, sample size, and minimum allowable MAF) by setting the corresponding parameter
as NULL and providing values for the other 2 parameters.
}
\usage{
powerEQTL.ANOVA(MAF,
                deltaVec=c(-0.13, 0.13),
                n=200,
                power = NULL,
                sigma = 0.13,
                FWER = 0.05,
                nTests = 200000,
                n.lower = 4,
                n.upper = 1e+30)
}
\arguments{
  \item{MAF}{
numeric. Minor allele frequency.
}
  \item{deltaVec}{
numeric. A vector having 2 elements. The first element is equal to
\eqn{\mu_2-\mu_1}{mu_2 - mu_1} and the second elementis equalt to
\eqn{\mu_3-\mu_2}{mu_3 - mu_2}, where \eqn{\mu_1}{mu_1} is
the mean gene expression level for the mutation homozygotes,
\eqn{\mu_2}{mu_2} is the mean gene expression level for the heterozygotes,
and \eqn{\mu_3}{mu_3} is the mean gene expression level for the 
wild-type gene expression level.
}
  \item{n}{
integer. Total number of subjects.
  }
  \item{power}{
numeric. Power for testing if 3 genotypes have the same mean 
gene expression levels.
  }
  \item{sigma}{
numeric. Standard deviation of the random error.
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
If we would like to test potential non-linear relationship between genotype of a SNP and expression of a gene, we can use un-balanced one-way ANOVA. Actually, an article published by the GTEx Consortium in 2013 used this approach. 

Suppose there are \eqn{k=3} groups of subjects: (1) mutation homozygotes; (2) heterozygotes; and (3) wildtype homozygotes. We would like to test if the mean expression \eqn{\mu_i}, \eqn{i=1, \ldots, k}, of the gene is the same among the \eqn{k} groups of subjects. We can use the following one-way ANOVA model to characterize the relationship between observed gene expression level \eqn{y_{ij}} and the population mean expression level \eqn{\mu_i}:
\deqn{y_{ij} = \mu_i + \epsilon_{ij}, \quad \epsilon_{ij} \sim N\left(0, \sigma^2\right),}
where \eqn{i=1,\ldots, k}, \eqn{j= 1, \ldots, n_i},
\eqn{y_{ij}} is the observed gene expression level for the \eqn{j}-th subject in the \eqn{i}-th group, \eqn{\mu_i} is the mean gene expression level of the \eqn{i}-th group, \eqn{\epsilon_{ij}} is the random error, \eqn{k} is the number of groups, \eqn{n_i} is the number of subjects in the \eqn{i}-th group. 
Denote the total number of subjects as \eqn{N = \sum_{i=1}^{k} n_i}.
That is, we have \eqn{n_1} mutation homozygotes, \eqn{n_2} heterozygotes, and \eqn{n_3} wildtype homozygotes.

We would like to test the null hypothesis \eqn{H_0} and alternative hypothesis \eqn{H_1}:
\deqn{H_0:  \mu_1 = \mu_2 = \mu_3,}
\deqn{H_1:  \mbox{not all means are the same}.}



According to O'Brien and Muller (1993), the power calculation formula for unbalanced one-way ANOVA is
\deqn{power=Pr\left(\left.F\geq F_{1-\alpha}\left(k-1, N-k\right)\right| F\sim
F_{k-1, N-k, \lambda}\right),}{power = Pr(F >= F(1 - alpha, k - 1, N - k)| F ~ F(k - 1, N - k, lambda)),}
where \eqn{k=3}{k = 3} is the number of groups of subjects, \eqn{N}{N} is the total number
of subjects, \eqn{F_{1-\alpha}\left(k-1, N-k\right)}{F_{1 - alpha}(k - 1, N - k)} is the 
\eqn{100(1-\alpha)}{100 * (1 - alpha)}-th percentile of central F distribution with degrees of freedoms \eqn{k-1}{k - 1} and \eqn{N-k}{N - k}, and \eqn{F_{k-1, N-k, \lambda}}{F_{k - 1, N - k, lambda}} is the non-central F distribution
with degrees of freedoms \eqn{k-1}{k - 1} and \eqn{N-k}{N - k} and non-central parameter (ncp)
\eqn{\lambda}{lambda}. The ncp \eqn{\lambda}{lambda} is equal to
\deqn{\lambda=\frac{N}{\sigma^2}\sum_{i=1}^{k} w_i \left(\mu_i-\mu\right)^2,}{
lambda = N * sum(wi * (mu_i - mu)^2, i = 1,.., k)/sigma^2,
} where \eqn{\mu_i}{mu_i} is the mean gene expression level
for the \eqn{i}{i}-th group of subjects, \eqn{w_i}{w_i} is the weight for the \eqn{i}{i}-th group of subjects, \eqn{\sigma^2}{sigma^2} is the variance of the random errors in ANOVA (assuming each group has equal variance), and
\eqn{\mu}{mu} is the weighted mean gene expression level
\deqn{\mu=\sum_{i=1}^{k}w_i \mu_i.}{mu = sum(w_i * mu_i, i = 1, ..., k).}
The weights \eqn{w_i=n_i/N}{w_i=n_i/N} are the sample proportions for the 3 groups of subjects, where \eqn{N=n_1+n_2+n_3} is the total number of subjects. Hence, 
\eqn{\sum_{i=1}^{3}w_i = 1}{sum(w_i, i = 1, 2, 3) = 1}.
Based on Hardy-Weinberg Equilibrium, we have
\eqn{w_1 = \theta^2}, \eqn{w_2 = 2\theta(1-\theta)}, and \eqn{w_3 = (1-\theta)^2}, 
where \eqn{\theta} is MAF.

Without loss of generality, we set \eqn{\mu_1 = -\delta_1}, \eqn{\mu_2=0},
and \eqn{\mu_3 = \delta_2}.

}
\value{
power if the input parameter \code{power = NULL}.

sample size if the input parameter \code{n = NULL};

minimum detectable slope if the input parameter \code{slope = NULL}.
}
\references{
The GTEx Consortium. 
The Genotype-Tissue Expression (GTEx) project. 
Nature Genetics, 45:580-585, 2013.

O'Brien, RG and Muller, KE. 
Unified power analysis for t-tests through multivariate hypotheses. 
In LK Edwards, editor, Applied Analysis of Variance in Behavioral Science,
pages 297-344. New York: Dekker, 1993.

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
%\seealso{
%
%}
\examples{

# calculate power
powerEQTL.ANOVA(MAF = 0.1,
                deltaVec = c(-0.13, 0.13),
                n = 282,
                power = NULL,
                sigma = 0.13,
                FWER = 0.05,
                nTests = 200000)
                
# calculate sample size
powerEQTL.ANOVA(MAF = 0.1,
                deltaVec = c(-0.13, 0.13),
                n = NULL,
                power = 0.8,
                sigma = 0.13,
                FWER = 0.05,
                nTests = 200000) 
                
# calculate minimum allowable MAF
powerEQTL.ANOVA(MAF = NULL,
                deltaVec = c(-0.13, 0.13),
                n = 282,
                power = 0.8,
                sigma = 0.13,
                FWER = 0.05,
                nTests = 200000)                 
                       
                       
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

