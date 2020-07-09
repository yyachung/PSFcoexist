#Non-bootstrap demography function
#AC 200225

nobsDemog <- function(dat, levels.var, germ.var, surv.var, growth.var){
  #Write empty data frame for output
  nrows<-length(unique(dat[, levels.var]))
  estimates_noBS<-data.frame(treat.combo=unique(dat[, levels.var]),
                        germ.mean=rep(0,nrows),
                        germ.95CIlo=rep(0,nrows),
                        germ.95CIhi=rep(0,nrows),
                        surv.mean=rep(0,nrows),
                        surv.95CIlo=rep(0,nrows),
                        surv.95CIhi=rep(0,nrows),
                        mass.mean=rep(0,nrows),
                        mass.95CIlo=rep(0,nrows),
                        mass.95CIhi=rep(0,nrows),
                        massContrib.mean=rep(0,nrows),
                        massContrib.95CIlo=rep(0,nrows),
                        massContrib.95CIhi=rep(0,nrows),
                        massContrib.SD=rep(0,nrows),
                        massContrib.SE=rep(0,nrows))
  #Bootstrap
  for (i in 1:nrows){
    subdata<-dat[which(dat[, levels.var]==estimates_noBS$treat.combo[i]),]
    masscontrib.vector<-na.omit(subdata[, germ.var]/10) %o% subdata[, surv.var] %o% na.omit(subdata[, growth.var])
    estimates_noBS$germ.mean[i]<-mean(na.omit(subdata[, germ.var]/10))
    estimates_noBS$surv.mean[i]<-mean(subdata[, surv.var])
    estimates_noBS$mass.mean[i]<-mean(na.omit(subdata[, growth.var]))
    estimates_noBS$massContrib.mean[i]<-mean(masscontrib.vector*100)
    estimates_noBS$germ.95CIlo[i]<-quantile(na.omit(subdata[, germ.var]/10),probs=0.025)
    estimates_noBS$germ.95CIhi[i]<-quantile(na.omit(subdata[, germ.var]/10),probs=0.975)
    estimates_noBS$surv.95CIlo[i]<-quantile(subdata[, surv.var],probs=0.025)
    estimates_noBS$surv.95CIhi[i]<-quantile(subdata[, surv.var],probs=0.975)
    estimates_noBS$mass.95CIlo[i]<-quantile(na.omit(subdata[, growth.var]),probs=0.025)
    estimates_noBS$mass.95CIhi[i]<-quantile(na.omit(subdata[, growth.var]),probs=0.975)
    estimates_noBS$massContrib.95CIlo[i]<-quantile(masscontrib.vector*100,probs=0.025)
    estimates_noBS$massContrib.95CIhi[i]<-quantile(masscontrib.vector*100,probs=0.975)
    estimates_noBS$massContrib.SD[i]<-sd(masscontrib.vector*100)
    estimates_noBS$massContrib.SE[i]<-sd(masscontrib.vector*100)/sqrt(length(masscontrib.vector))
  }
  return(estimates_noBS)
}
