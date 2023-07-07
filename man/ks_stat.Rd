\name{ks_stat}
\alias{ks_stat}
\title{Kolmogorov-Smirnov statistics}
\usage{
ks_stat(x,z)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the covariates. }
\item{z}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
}
\value{
A vector of Kolmogorov-Smirnov statistics for each column of x.
}
\description{
Compute the Kolmogorov-Smirnov statistics for each individual covariate in x.
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
ks_stat(Xm,M0$z)
detach(SSRI)
}

