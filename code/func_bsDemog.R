#Bootstrap demography function
#AC 191213

bsDemog <- function(dat, nb, bs.var, germ.var, surv.var, growth.var){
  #Write empty data frame for output
  nrows<-length(unique(dat[, bs.var]))
  estimates<-data.frame(treat.combo=unique(dat[, bs.var]),
                        germ.boot.mean=rep(0,nrows),
                        germ.boot.95CIlo=rep(0,nrows),
                        germ.boot.95CIhi=rep(0,nrows),
                        surv.boot.mean=rep(0,nrows),
                        surv.boot.95CIlo=rep(0,nrows),
                        surv.boot.95CIhi=rep(0,nrows),
                        mass.boot.mean=rep(0,nrows),
                        mass.boot.95CIlo=rep(0,nrows),
                        mass.boot.95CIhi=rep(0,nrows),
                        massContrib.mean=rep(0,nrows),
                        massContrib.95CIlo=rep(0,nrows),
                        massContrib.95CIhi=rep(0,nrows))
  bsVectors<-vector(mode = "list", length = nrows)
  #Bootstrap
  for (i in 1:nrows){
    subdata<-dat[which(dat[, bs.var]==estimates$treat.combo[i]),]
    germ.boot.mean.vector<-NULL
    surv.boot.mean.vector<-NULL
    mass.boot.mean.vector<-NULL
    massContrib.boot.mean.vector<-NULL
    for (j in 1:nb){
      germ.boot.mean.vector[j]<-mean(sample(na.omit(subdata[, germ.var]),replace=T))/10 #untested germination not included; /10 for germ rate
      surv.boot.mean.vector[j]<-mean(sample(subdata[, surv.var],replace=T))
      mass.boot.mean.vector[j]<-mean(sample(na.omit(subdata[, growth.var]),replace=T))#Dead plants not included in growth calculation
      massContrib.boot.mean.vector[j]<-100*germ.boot.mean.vector[j]*surv.boot.mean.vector[j]*mass.boot.mean.vector[j] #Directly estimated contributed growth shows strong kurtosis
    }
    estimates$germ.boot.mean[i]<-mean(germ.boot.mean.vector)
    estimates$surv.boot.mean[i]<-mean(surv.boot.mean.vector)
    estimates$mass.boot.mean[i]<-mean(mass.boot.mean.vector)
    estimates$massContrib.mean[i]<-mean(massContrib.boot.mean.vector)
    estimates$germ.boot.95CIlo[i]<-quantile(germ.boot.mean.vector,probs=0.025)
    estimates$germ.boot.95CIhi[i]<-quantile(germ.boot.mean.vector,probs=0.975)
    estimates$surv.boot.95CIlo[i]<-quantile(surv.boot.mean.vector,probs=0.025)
    estimates$surv.boot.95CIhi[i]<-quantile(surv.boot.mean.vector,probs=0.975)
    estimates$mass.boot.95CIlo[i]<-quantile(mass.boot.mean.vector,probs=0.025)
    estimates$mass.boot.95CIhi[i]<-quantile(mass.boot.mean.vector,probs=0.975)
    estimates$massContrib.95CIlo[i]<-quantile(massContrib.boot.mean.vector,probs=0.025)
    estimates$massContrib.95CIhi[i]<-quantile(massContrib.boot.mean.vector,probs=0.975)
    bsVectors[[i]]<-massContrib.boot.mean.vector
  }
  return(list(estimates,bsVectors))
}
