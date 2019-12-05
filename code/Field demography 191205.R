#Code to try and combine all life stage transitions in field PSF experiment
#AC 5 Dec 2019

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")

#Data import------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Try estimating params for all conditions in 2018--------
#Write empty data frame
nrows<-length(unique(all2018$Trans.Don.Treat))
estimates<-data.frame(Trans.Don.Treat=unique(all2018$Trans.Don.Treat),
                      germ.mean=rep(0,nrows),
                      germ.boot.mean=rep(0,nrows),
                      surv.mean=rep(0,nrows),
                      surv.boot.mean=rep(0,nrows),
                      mass.mean=rep(0,nrows),
                      mass.boot.mean=rep(0,nrows))
#Calculate bootstrap means and regular means
nb<-5000
for (i in 1:nrows){
  subdata<-subset(all2018,Trans.Don.Treat==estimates$Trans.Don.Treat[i])
  germ.boot.mean<-NULL
  surv.boot.mean<-NULL
  mass.boot.mean<-NULL
  for (j in 1:nb){
    germ.boot.mean[j]<-mean(sample(subdata$Germ,replace=T))
    surv.boot.mean[j]<-mean(sample(subdata$Surv,replace=T))
    mass.boot.mean[j]<-mean(sample(na.omit(subdata$Abv),replace=T))#Dead plants not included in mass calculation
  }
  estimates$germ.boot.mean[i]<-mean(germ.boot.mean)
  estimates$surv.boot.mean[i]<-mean(surv.boot.mean)
  estimates$mass.boot.mean[i]<-mean(mass.boot.mean)
  estimates$germ.mean[i]<-mean(subdata$Germ)
  estimates$surv.mean[i]<-mean(subdata$Surv)
  estimates$mass.mean[i]<-mean(na.omit(subdata$Abv))
}
#Bootstrap means pretty much approach the sample means, which I guess is expected. Not sure what to do with that just yet.
