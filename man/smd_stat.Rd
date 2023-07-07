\name{smd_stat}
\alias{smd_stat}
\title{Standardized mean differences}
\usage{
smd_stat(x,z,xf,zf)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the covariates after matching. }
\item{z}{A vector whose ith coordinate is 1 for a matched treated unit and is 0 for a matched control.}
\item{xf}{A matrix with length(zf) rows giving the covariates before matching. }
\item{zf}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
}
\value{
A vector of standardized mean differences for each column of x.
}
\description{
Compute the standardized mean differences for each individual covariate in x, i.e., the mean difference of each column of x is standarzied by the pooled standard deviation calculated based on the original data before matching (xf and zf).
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
smd_stat(Xm,M0$z,X,z)
detach(SSRI)
}

