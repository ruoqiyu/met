gfks1_stat<-function(x,z){
  
  nvar=dim(x)[2]
  n1=sum(z)
  n0=sum(1-z)
  ks=numeric(nvar)
  where=numeric(nvar)
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
    dif=abs(cumsum(xx1i)/n1-cumsum(xx0i)/n0)
    ks[i]=max(dif)
    where[i]=sort(unique(xi))[which.max(dif)]
  }
  list(gfks=ks,where=where)
}