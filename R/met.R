met<-function(x,z,method='GFKS-1',prob=c(0,0.11,0.35,0.65,0.89,1),continuous=rep(FALSE,dim(x)[2]),nperm=1000,xf=NULL,zf=NULL){
  nvar=dim(x)[2]
  if (method%in%c('GFKS-1','GFKS-2')){
    num=length(prob)-1
    xx=x
    xxf=xf
    for (i in 1:nvar){
      if (continuous[i]){
        cutoff=stats::quantile(x[z==1,i],probs=prob)
        tb=table(cutoff)
        if (any(tb!=1)){
          cutoff=sort(c(unique(cutoff),unique(cutoff)[tb!=1]*0.999999))
        }
        c=cut(x[,i],breaks=cutoff,labels = 1:num)
        if (sum(is.na(c))>0){
          c[x[,i]>max(x[z==1,i])]<-num
          c[x[,i]<=min(x[z==1,i])]<-1
        }
        xx[,i]=c
        cf=cut(xf[,i],breaks=cutoff,labels = 1:num)
        if (sum(is.na(cf))>0){
          cf[xf[,i]>max(x[z==1,i])]<-num
          cf[xf[,i]<=min(x[z==1,i])]<-1
        }
        xxf[,i]=cf
      }
    }
  }

  if (method=='GFKS-1') {
    res=gfks1_stat(xx,z)
    stat=res$gfks
    res$direction=rep(1,length(stat))
  }
  else if (method=='GFKS-2'){
    res=gfks2_stat(xx,z)
    stat=res$gfks
  }
  else if (method=='SMD') stat=smd_stat(x,z,xf,zf)
  else if (method=='t') stat=t_stat(x,z)
  else if (method=='wilcoxon') stat=wilcoxon_stat(x,z)
  else if (method=='KS') stat=ks_stat(x,z)

  statlist=matrix(0,nrow=nperm,ncol=length(stat))

  #  set.seed(1)
  n=length(z)
  for (i in 1:nperm){
    znew=numeric(n)
    ix=sample(1:n,sum(z))
    znew[ix]=1
    if (method=='GFKS-1') statlist[i,]=gfks1_stat(xx,znew)$gfks
    if (method=='GFKS-2') statlist[i,]=gfks2_stat(xx,znew)$gfks
    if (method=='SMD') statlist[i,]=smd_stat(x,znew,xf,zf)
    if (method=='t') statlist[i,]=t_stat(x,znew)
    if (method=='wilcoxon') statlist[i,]=wilcoxon_stat(x,znew)
    if (method=='KS') statlist[i,]=ks_stat(x,znew)
  }

  pvalues0=matrix(ncol=length(stat),nrow=nperm)
  for (i in 1:nperm){
    znew=numeric(n)
    ix=sample(1:n,sum(z))
    znew[ix]=1
    if (method=='GFKS-1') sta=gfks1_stat(xx,znew)$gfks
    else if (method=='GFKS-2') sta=gfks2_stat(xx,znew)$gfks
    else if (method=='SMD') sta=smd_stat(x,znew,xf,zf)
    else if (method=='t') sta=t_stat(x,znew)
    else if (method=='wilcoxon') sta=wilcoxon_stat(x,znew)
    else if (method=='KS') sta=ks_stat(x,znew)
    for (kk in 1:length(stat)){
      if (method=='wilcoxon') pvalues0[i,kk]=ecdf(-abs(statlist[,kk]-sum(z)*sum(1-z)/2))(-abs(sta[kk]-sum(z)*sum(1-z)/2))
      else if (method%in%c('SMD','t')) pvalues0[i,kk]=ecdf(-abs(statlist[,kk]))(-abs(sta[kk]))
      else pvalues0[i,kk]=ecdf(-statlist[,kk])(-sta[kk])
    }
  }

  pvalue0=numeric(length(stat))
  for (kk in 1:length(stat)){
    if (method=='wilcoxon') pvalue0[kk]=ecdf(-abs(statlist[,kk]-sum(z)*sum(1-z)/2))(-abs(stat[kk]-sum(z)*sum(1-z)/2))
    else if (method%in%c('SMD','t')) pvalue0[kk]=ecdf(-abs(statlist[,kk]))(-abs(stat[kk]))
    else pvalue0[kk]=ecdf(-statlist[,kk])(-stat[kk])
  }

  minp=apply(pvalues0,1,min)
  minpp=ecdf(minp)(min(pvalue0))
  ix=which(pvalue0==min(pvalue0))

  if (method%in%c('GFKS-1','GFKS-2')){
    location=res$where[ix]
    direction=res$direction[ix]
    if (method=='GFKS1') variable=as.numeric(as.numeric(xxf[,ix])<=location)
    else{
      ulocation=location[[1]]
      udirection=direction[[1]]
      uix=ix[1]
      if (uix<=nvar) variable=as.numeric(as.numeric(xxf[,uix])<=as.numeric(ulocation))
      else{
        nnn=cumsum(nvar:1)
        nnnix=max(which(uix>nnn))
        ix1=nnnix
        ix2=uix-nnn[nnnix]+ix1
        if (udirection==1) variable= as.numeric(as.numeric(xxf[,ix1])<=ulocation[1] & as.numeric(xxf[,ix2])<=ulocation[2])
        else if (udirection==2) variable= as.numeric(as.numeric(xxf[,ix1])<=ulocation[1] & as.numeric(xxf[,ix2])>ulocation[2])
        else if (udirection==3) variable= as.numeric(as.numeric(xxf[,ix1])>ulocation[1] & as.numeric(xxf[,ix2])<=ulocation[2])
        else if (udirection==4) variable= as.numeric(as.numeric(xxf[,ix1])>ulocation[1] & as.numeric(xxf[,ix2])>ulocation[2])
      }
    }
    return (list(statistic=stat,pvalue=minpp,imbalance=ix,location=location,direction=direction,ind.pvalue=pvalue0,null.pvalues=pvalues0,variable=variable))
  }
  else return (list(statistic=stat,pvalue=minpp,imbalance=ix,location=NULL,direction=NULL,ind.pvalue=pvalue0,null.pvalues=pvalues0,variable=NULL))
}
