\name{t_stat}
\alias{t_stat}
\title{Two sample t statistics}
\usage{
t_stat(x,z)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the covariates. }
\item{z}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
}
\value{
A vector of two-sample t statistics for each column of x.
}
\description{
Compute the two-sample t statistics for each individual covariate in x.
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
t_stat(Xm,M0$z)
detach(SSRI)
}

