\name{gfks1_stat}
\alias{gfks1_stat}
\title{Univariate GFKS statistics}
\usage{
gfks1_stat(x,z)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the discretized covariates. }
\item{z}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
}
\value{
\item{gfks}{A vector of GFKS statistics for each discrete covariate.}
\item{where}{A vector of locations that the GFKS statistics occur, i.e., the location of the maximum difference of the empirical distribution function of each discrete covariate in the treated and control groups.}
}
\description{
Compute the GFKS statistics and identify its location for each column of x. This function and its use are discussed in Yu (2020).
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
gfks1_stat(Xm,M0$z)
detach(SSRI)
}

\references{
Yu, R. (2020) Evaluating and Improving a Matched Comparison of Antidepressants and Bone Density. Under revision.
}
