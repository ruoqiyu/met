\name{gfks2_stat}
\alias{gfks2_stat}
\title{Bivariate GFKS statistics}
\usage{
gfks2_stat(x,z)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the discretized covariates. }
\item{z}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
}
\value{
\item{gfks}{A vector of GFKS statistics for all marginal and bivariate distributions.}
\item{where}{A vector of locations v that the GFKS statistics occur, i.e., the location of the maximum difference of the empirical distribution function in the treated and control groups.}
\item{direction}{A vector of directions indicating which quadrant the GFKS statistics occur -- 1 for '<=' and '<=', 2 for '<=' and '>', 3 for '>' and '<=', 4 for '>' and '>'.}
}
\description{
Compute the GFKS statistics and identify its location for all marginal distributions and bivariate distributions for x. This function and its use are discussed in Yu (2020).
}

\examples{
library(optmatch)
data("SSRI")
attach(SSRI)
X<-cbind(female,black,education)
dist<-DiPs::maha_dense(z,X)
o<-DiPs::match(z, dist, SSRI)
M0<-o$data
Xm<-cbind(M0$female,M0$black,M0$education)
gfks2_stat(Xm,M0$z)
detach(SSRI)
}

\references{
Yu, R. (2020) Evaluating and Improving a Matched Comparison of Antidepressants and Bone Density. Under revision.
}
