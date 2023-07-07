ks_stat<-function(x,z){
  nvar=dim(x)[2]
  x1=x[z==1,]
  x0=x[z==0,]
  ks=numeric(nvar)
  where=numeric(nvar)
  for (i in 1:nvar){
    ks[i]=stats::ks.test(x1[,i],x0[,i])$statistic
  }
  ks
}
