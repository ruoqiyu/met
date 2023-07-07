gfks2_stat<-function(x,z){
  nvar=dim(x)[2]
  x1=x[z==1,]
  x0=x[z==0,]
  n1<-sum(z)
  n0=sum(1-z)

  #marginal
  ks=numeric(nvar+nvar*(nvar-1)/2)
  loc=vector(mode='list',length=nvar+nvar*(nvar-1)/2)
  direction=vector(mode='list',length=nvar+nvar*(nvar-1)/2)
  for (i in 1:nvar){
    xi=x[,i]
    nlevel=length(unique(xi))
    xxi=as.numeric(factor(xi,labels=1:nlevel))
    tb1=table(xxi[z==1])
    tb0=table(xxi[z==0])
    if (length(tb1)==nlevel) xx1i=tb1
    else{
      xx1i=numeric(nlevel)
      xx1i[as.numeric(names(tb1))]=tb1
    }
    if (length(tb0)==nlevel) xx0i=tb0
    else{
      xx0i=numeric(nlevel)
      xx0i[as.numeric(names(tb0))]=tb0
    }
    di=abs(cumsum(xx1i)/n1-cumsum(xx0i)/n0)
    ks[i]=max(di)
    loc[[i]]=sort(unique(xi))[which.max(di)]
    direction[[i]]=1
  }

  edf_compare2<-function(dd,zz){
    if (is.data.frame(dd)) dd<-data.matrix(dd)
    nlevel1=length(unique(dd[,1]))
    nlevel2=length(unique(dd[,2]))
    d=dd
    d[,1]=as.numeric(factor(d[,1],labels=1:nlevel1))
    d[,2]=as.numeric(factor(d[,2],labels=1:nlevel2))
    dd1<-d[which(zz==1),]
    dd0<-d[which(zz==0),]
    n1=sum(zz)
    n0=sum(1-zz)

    tb1<-table(dd1[,1],dd1[,2])/n1
    if (nrow(tb1)<nlevel1){
      A1=matrix(0,nrow=nlevel1,ncol=ncol(tb1))
      A1[as.numeric(unique(dd1[,1])),]=tb1
    }
    if (ncol(tb1)<nlevel2){
      tb1=matrix(0,nrow=nlevel1,ncol=nlevel2)
      tb1[,as.numeric(unique(dd1[,2]))]=A1
    }
    A1<-tb1
    B1<-A1
    for (i in 1:nlevel1) B1[i,]=cumsum(A1[i,])
    C11<-B1
    for (i in 1:nlevel2) C11[,i]=cumsum(B1[,i])
    C12=matrix(C11[,nlevel2],nrow=nlevel1,ncol=nlevel2)-C11
    C13=matrix(C11[nlevel1,],nrow=nlevel1,ncol=nlevel2,byrow=T)-C11
    C14=matrix(C11[nlevel1,nlevel2],nrow=nlevel1,ncol=nlevel2)-C12-C13-C11

    tb0<-table(dd0[,1],dd0[,2])/n0
    if (nrow(tb0)<nlevel1){
      A0=matrix(0,nrow=nlevel1,ncol=ncol(tb0))
      A0[as.numeric(unique(dd0[,1])),]=tb0
    }
    if (ncol(tb0)<nlevel2){
      tb0=matrix(0,nrow=nlevel1,ncol=nlevel2)
      tb0[,as.numeric(unique(dd0[,2]))]=A0
    }
    A0<-tb0
    B0<-A0
    for (i in 1:nlevel1) B0[i,]=cumsum(A0[i,])
    C01<-B0
    for (i in 1:nlevel2) C01[,i]=cumsum(B0[,i])
    C02=matrix(C01[,nlevel2],nrow=nlevel1,ncol=nlevel2)-C01
    C03=matrix(C01[nlevel1,],nrow=nlevel1,ncol=nlevel2,byrow=T)-C01
    C04=matrix(C01[nlevel1,nlevel2],nrow=nlevel1,ncol=nlevel2)-C02-C03-C01

    abs1=abs(C11-C01)
    abs2=abs(C12-C02)
    abs3=abs(C13-C03)
    abs4=abs(C14-C04)
    mm=pmax(abs1,abs2,abs3,abs4)
    r=max(mm)
    loc=which(mm==r, arr.ind = TRUE)
    direction=which(c(abs1[loc],abs2[loc],abs3[loc],abs4[loc])==r)
    list(dif=r, loc=c(sort(unique(dd[,1]))[loc[1,1]],sort(unique(dd[,2]))[loc[1,2]]),direction=direction)
  }

  pin=nvar
  #joint distribution
  for(k1 in 1:(nvar-1)){
    for (k2 in (k1+1):nvar){
      pin=pin+1
      r12=edf_compare2(x[,c(k1,k2)],z)
      ks[pin]<-r12$dif
      loc[[pin]]<-r12$loc
      direction[[pin]]=r12$direction
    }
  }
  list(gfks=ks, where=loc, direction=direction)
}
