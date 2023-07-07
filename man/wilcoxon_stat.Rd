\name{wilcoxon_stat}
\alias{wilcoxon_stat}
\title{Wilcoxon rank statistics}
\usage{
wilcoxon_stat(x,z)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the covariates. }
\item{z}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
}
\value{
A vector of Wilcoxon rank statistics for each column of x.
}
\description{
Compute the Wilcoxon rank statistics for each individual covariate in x.
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
wilcoxon_stat(Xm,M0$z)
detach(SSRI)
}

