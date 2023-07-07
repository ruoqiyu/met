smd_stat<-function(x,z,xf,zf){

  nvar=dim(x)[2]
  x1=x[z==1,]
  x0=x[z==0,]
  smd=numeric(nvar)
  where=numeric(nvar)
  smd=DiPs::check(xf,x,zf,z)[,4]
  smd
}
