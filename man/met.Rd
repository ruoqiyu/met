\name{met}
\alias{met}
\title{Balance evaluation of matched samples}
\usage{
met(x,z,method='GFKS-1',prob=c(0,0.11,0.35,0.65,0.89,1),
continuous=rep(FALSE,dim(x)[2]),nperm=1000,xf=NULL,zf=NULL)
}
\arguments{
\item{x}{A matrix with length(z) rows giving the discretized covariates. }
\item{z}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{method}{Balance evaluation methods, can be chosen from 'GFKS-1', 'GFKS-2','SMD','KS','t','wilcoxon'.}
\item{prob}{If method is chosen as 'GFKS-1'or 'GFKS-2', continuous variables in x are discretized based on prob, i.e., only the quantiles at prob are considered.}
\item{continuous}{A vector of length dim(x)[2] indicating whether the covariates are discrete (F) or continuous (T). If continuous, GFKS statistics are evaluted at the quantiles at prob.}
\item{nperm}{Number of simulated randomized experiments to approximate the p-values.}
\item{xf}{A matrix with length(zf) rows giving the covariates before matching. xf should be specified if method is chosen as 'GFKS-1', 'GFKS-2', 'SMD'.}
\item{zf}{A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. zf should be specified if method is chosen as 'GFKS-1', 'GFKS-2', 'SMD'.}
}
\value{
\item{statistic}{A summary statistic to combine all statistics -- the minimum p-value.}
\item{pvalue}{A vector of locations v that the GFKS statistics occur, i.e., the location of the maximum difference of the empirical distribution function in the treated and control groups.}
\item{imbalance}{A vector indicating which marginal distribution or joint distribution is most imbalanced. Not NULL if method is chosen as 'GFKS-1' or 'GFKS-2'.}
\item{location}{Location that statistic occur. Not NULL if method is chosen as 'GFKS-1' or 'GFKS-2'.}
\item{direction}{Direction indicating which quadrant statistic occur. 1 for '<=' and '<=', 2 for '<=' and '>', 3 for '>' and '<=', 4 for '>' and '>'. Not NULL if method is chosen as 'GFKS-1' or 'GFKS-2'.}
\item{ind.pvalue}{P-values for the observational match.}
\item{null.pvalues}{P-values for the simulated randomized experiments.}
\item{variable}{A new binary variable that the next iteration match should balance. Not NULL if method is chosen as 'GFKS-1' or 'GFKS-2'.}
}
\description{
Evaluate a matched sample with minimum p-value of multiple measures of imbalance. This function and its use are discussed in Yu (2020).
}

\examples{
library(optmatch)
data("SSRI")
attach(SSRI)
X<-cbind(female,black,education,age,bmi)
dist<-DiPs::maha_dense(z,X)
o<-DiPs::match(z, dist, SSRI)
M0<-o$data
Xm<-cbind(M0$female,M0$black,M0$education,M0$age,M0$bmi)
met(Xm,M0$z,'GFKS-1',continuous=c(FALSE,FALSE,FALSE,TRUE,TRUE),nperm=100,xf=X,zf=z)
detach(SSRI)

\donttest{
#real data application in Yu(2020)
library(DiPs)
data("SSRI")
attach(SSRI)
X<-cbind(age,female,black,hispanic,povertyNA,povertyFill,education,height,weight,
bmi,cotinine,hd,diabetes,insurance,weighmore,weighless,weightchange,physicalact,dietsup)
XX<-cbind(age,female,black,hispanic,povertyNA,povertyFill,education,height,weight,
bmi,cotinine,hd,diabetes,insurance,weighmore,weighless,weightchange,physicalact,dietsup,pr)
detach(SSRI)

dat=SSRI
#basic match
dist<-maha_dense(dat$z,X)
dist<-addcaliper(dist, dat$z, dat$pr, c(-.2,.2), stdev = TRUE, penalty = 1000)
o<-match(dat$z,dist, dat,fine=factor(dat$education), ncontrol = 4)
M0<-o$data
Xm<-subset(M0, select=c('age','female','black','hispanic','povertyNA','povertyFill',
'education','height','weight','bmi','cotinine','hd','diabetes','insurance','weighmore','weighless',
'weightchange','physicalact','dietsup','pr'))
btb0<-check(XX,Xm,dat$z,M0$z)
round(btb0,3)

Result0_t=met(Xm,M0$z,method='t',nperm=2000)
Result0_t$pvalue

Result0_smd=met(Xm,M0$z,method='SMD',nperm=2000,xf=XX,zf=dat$z)
Result0_smd$pvalue

Result0_w=met(Xm,M0$z,method='wilcoxon',nperm=2000)
Result0_w$pvalue

Result0_ks=met(Xm,M0$z,method='KS',nperm=2000)
Result0_ks$pvalue

Result0=met(Xm,M0$z,method='GFKS-2',continuous=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),nperm=2000,xf=XX,zf=dat$z)
Result0$statistic
Result0$pvalue
Result0$imbalance
Result0$location
Result0$direction

#iteration 1
dist<-maha_dense(dat$z,X)
dist<-addcaliper(dist, dat$z, dat$pr, c(-.2,.2), stdev = TRUE, penalty = 1000)
i1=as.factor(dat$diabetes<=0 & dat$weighmore<=0)
o<-match(dat$z, dist, dat, fine=as.numeric(factor(dat$education):i1), ncontrol = 4)
M1<-o$data
Xm<-subset(M1, select=c('age','female','black','hispanic','povertyNA','povertyFill',
'education','height','weight','bmi','cotinine','hd','diabetes','insurance','weighmore','weighless',
'weightchange','physicalact','dietsup','pr'))
btb1<-check(XX,Xm,dat$z,M1$z)
round(btb1,3)

Result1_t=met(Xm,M1$z,method='t',nperm=2000)
Result1_t$pvalue

Result1_smd=met(Xm,M1$z,method='SMD',nperm=2000,xf=XX,zf=dat$z)
Result1_smd$pvalue

Result1_w=met(Xm,M1$z,method='wilcoxon',nperm=2000)
Result1_w$pvalue

Result1_ks=met(Xm,M1$z,method='KS',nperm=2000)
Result1_ks$pvalue

Result1=met(Xm,M1$z,method='GFKS-2',continuous=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),nperm=2000,xf=XX,zf=dat$z)
Result1$statistic
Result1$pvalue
Result1$imbalance
Result1$location
Result1$direction

#iteration 2
dist<-maha_dense(dat$z,X)
dist<-addcaliper(dist, dat$z, dat$pr, c(-.2,.2), stdev = TRUE, penalty = 1000)
i2=as.factor(dat$weightchange<=quantile(dat$weightchange[dat$z==1],probs=0.65) &
as.numeric(dat$pr)<=quantile(dat$pr[dat$z==1],probs=0.89))
o<-match(dat$z, dist, dat, fine=as.numeric(factor(dat$education):i1:i2), ncontrol = 4)
M2<-o$data
Xm<-subset(M2, select=c('age','female','black','hispanic','povertyNA','povertyFill',
'education','height','weight','bmi','cotinine','hd','diabetes','insurance','weighmore','weighless',
'weightchange','physicalact','dietsup','pr'))
btb2<-check(XX,Xm,dat$z,M2$z)
round(btb2,3)

Result2_t=met(Xm,M2$z,method='t',nperm=2000)
Result2_t$pvalue

Result2_smd=met(Xm,M2$z,method='SMD',nperm=2000,xf=XX,zf=dat$z)
Result2_smd$pvalue

Result2_w=met(Xm,M2$z,method='wilcoxon',nperm=2000)
Result2_w$pvalue

Result2_ks=met(Xm,M2$z,method='KS',nperm=2000)
Result2_ks$pvalue

Result2=met(Xm,M2$z,method='GFKS-2',continuous=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),nperm=2000,xf=XX,zf=dat$z)
Result2$statistic
Result2$pvalue
Result2$imbalance
Result2$location
Result2$direction

#iteration 3
dist<-maha_dense(dat$z,X)
dist<-addcaliper(dist, dat$z, dat$pr, c(-.2,.2), stdev = TRUE, penalty = 1000)
i3=as.factor(dat$weightchange<=quantile(dat$weightchange[dat$z==1],probs=0.35))
o<-match(dat$z, dist, dat, fine=as.numeric(factor(dat$education):i1:i2:i3), ncontrol = 4)
M3<-o$data
Xm<-subset(M3, select=c('age','female','black','hispanic','povertyNA','povertyFill',
'education','height','weight','bmi','cotinine','hd','diabetes','insurance','weighmore','weighless',
'weightchange','physicalact','dietsup','pr'))
btb3<-check(XX,Xm,dat$z,M3$z)
round(btb3,3)

Result3_t=met(Xm,M3$z,method='t',nperm=2000)
Result3_t$pvalue

Result3_smd=met(Xm,M3$z,method='SMD',nperm=2000,xf=XX,zf=dat$z)
Result3_smd$pvalue

Result3_w=met(Xm,M3$z,method='wilcoxon',nperm=2000)
Result3_w$pvalue

Result3_ks=met(Xm,M3$z,method='KS',nperm=2000)
Result3_ks$pvalue

Result3=met(Xm,M3$z,method='GFKS-2',continuous=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),nperm=2000,xf=XX,zf=dat$z)
Result3$statistic
Result3$pvalue
Result3$imbalance
Result3$location
Result3$direction
}
}

\references{
Yu, R. (2020) Evaluating and Improving a Matched Comparison of Antidepressants and Bone Density. Under revision.
}
