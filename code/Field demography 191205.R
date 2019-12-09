#Code to try and combine all life stage transitions in field PSF experiment
#AC 5 Dec 2019

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
library(ggpubr)
library(viridis)

#Data import------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Estimating params for all conditions in 2018--------
all2018$Trans.Don.Treat<-paste(all2018$Transplant,all2018$DonorSpp,all2018$Treatment,sep=".")
#Write empty data frame
nrows<-length(unique(all2018$Trans.Don.Treat))
estimates<-data.frame(Trans.Don.Treat=unique(all2018$Trans.Don.Treat),
                      germ.mean=rep(0,nrows),
                      germ.boot.mean=rep(0,nrows),
                      germ.boot.95CIlo=rep(0,nrows),
                      germ.boot.95CIhi=rep(0,nrows),
                      surv.mean=rep(0,nrows),
                      surv.boot.mean=rep(0,nrows),
                      surv.boot.95CIlo=rep(0,nrows),
                      surv.boot.95CIhi=rep(0,nrows),
                      mass.mean=rep(0,nrows),
                      mass.boot.mean=rep(0,nrows),
                      mass.boot.95CIlo=rep(0,nrows),
                      mass.boot.95CIhi=rep(0,nrows))
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
  estimates$germ.boot.95CIlo[i]<-quantile(germ.boot.mean,probs=0.025)
  estimates$germ.boot.95CIhi[i]<-quantile(germ.boot.mean,probs=0.925)
  estimates$surv.boot.95CIlo[i]<-quantile(surv.boot.mean,probs=0.025)
  estimates$surv.boot.95CIhi[i]<-quantile(surv.boot.mean,probs=0.925)
  estimates$mass.boot.95CIlo[i]<-quantile(mass.boot.mean,probs=0.025)
  estimates$mass.boot.95CIhi[i]<-quantile(mass.boot.mean,probs=0.925)
}#Takes a few seconds
#Bootstrap means pretty much approach the sample means, which I guess is expected. Not sure what to do with that just yet.

#2018 estimated contribution to biomass-----------------
estimates$massContrib.mean<-100*estimates$germ.boot.mean*estimates$surv.boot.mean*estimates$mass.boot.mean

#Propagate errors
estimates$massContrib.95CIlo<-sqrt(((estimates$germ.boot.mean-estimates$germ.boot.95CIlo)/estimates$germ.boot.mean)^2
                                   +((estimates$surv.boot.mean-estimates$surv.boot.95CIlo)/estimates$surv.boot.mean)^2
                                   +((estimates$mass.boot.mean-estimates$mass.boot.95CIlo)/estimates$mass.boot.mean)^2
)
estimates$massContrib.95CIhi<-sqrt(((estimates$germ.boot.mean-estimates$germ.boot.95CIhi)/estimates$germ.boot.mean)^2
                                   +((estimates$surv.boot.mean-estimates$surv.boot.95CIhi)/estimates$surv.boot.mean)^2
                                   +((estimates$mass.boot.mean-estimates$mass.boot.95CIhi)/estimates$mass.boot.mean)^2
)

#Visualize 2018 projections-----------------
#Add labels back to data
estimates$Transplant<-substr(estimates$Trans.Don.Treat,1,4)
estimates$Donorspp<-substr(estimates$Trans.Don.Treat,6,9)
estimates$Treatment<-substr(estimates$Trans.Don.Treat,11,100)
#Delete final row from data since it has no treatment label
estimates<-estimates[1:60,]
#Make Intra/Inter label
estimates$Type[which(estimates$Donorspp!=estimates$Transplant)]<-"Inter"
estimates$Type[which(estimates$Donorspp==estimates$Transplant)]<-"Intra"
estimates$Type <- factor(estimates$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
estimates$Treatment <- factor(estimates$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization

#Make figure
ggplot(data=estimates, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
  facet_grid(~Transplant)+
  geom_pointrange(aes(ymin=massContrib.mean-massContrib.95CIlo, ymax=massContrib.95CIhi+massContrib.mean)
                  , position = position_dodge(width = 1))+  #Most errors are too small to see
  scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
  theme_bw()

#Estimating params for all conditions in 2019 [need to make this a function later]--------
all2019$Trans.Don.Treat<-paste(all2019$Transplant,all2019$DonorSpp,all2019$Treatment,sep=".")
#Write empty data frame
nrows<-length(unique(all2019$Trans.Don.Treat))
estimates<-data.frame(Trans.Don.Treat=unique(all2019$Trans.Don.Treat),
                      germ.mean=rep(0,nrows),
                      germ.boot.mean=rep(0,nrows),
                      germ.boot.95CIlo=rep(0,nrows),
                      germ.boot.95CIhi=rep(0,nrows),
                      surv.mean=rep(0,nrows),
                      surv.boot.mean=rep(0,nrows),
                      surv.boot.95CIlo=rep(0,nrows),
                      surv.boot.95CIhi=rep(0,nrows),
                      mass.mean=rep(0,nrows),
                      mass.boot.mean=rep(0,nrows),
                      mass.boot.95CIlo=rep(0,nrows),
                      mass.boot.95CIhi=rep(0,nrows))
#Calculate bootstrap means and regular means
nb<-5000
for (i in 1:nrows){
  subdata<-subset(all2019,Trans.Don.Treat==estimates$Trans.Don.Treat[i])
  germ.boot.mean<-NULL
  surv.boot.mean<-NULL
  mass.boot.mean<-NULL
  for (j in 1:nb){
    germ.boot.mean[j]<-mean(sample(na.omit(subdata$Germ),replace=T))
    surv.boot.mean[j]<-mean(sample(subdata$Surv,replace=T))
    mass.boot.mean[j]<-mean(sample(na.omit(subdata$Abv),replace=T))#Dead plants not included in mass calculation
  }
  estimates$germ.boot.mean[i]<-mean(germ.boot.mean)
  estimates$surv.boot.mean[i]<-mean(surv.boot.mean)
  estimates$mass.boot.mean[i]<-mean(mass.boot.mean)
  estimates$germ.mean[i]<-mean(na.omit(subdata$Germ))
  estimates$surv.mean[i]<-mean(subdata$Surv)
  estimates$mass.mean[i]<-mean(na.omit(subdata$Abv))
  estimates$germ.boot.95CIlo[i]<-quantile(germ.boot.mean,probs=0.025)
  estimates$germ.boot.95CIhi[i]<-quantile(germ.boot.mean,probs=0.925)
  estimates$surv.boot.95CIlo[i]<-quantile(surv.boot.mean,probs=0.025)
  estimates$surv.boot.95CIhi[i]<-quantile(surv.boot.mean,probs=0.925)
  estimates$mass.boot.95CIlo[i]<-quantile(mass.boot.mean,probs=0.025)
  estimates$mass.boot.95CIhi[i]<-quantile(mass.boot.mean,probs=0.925)
}#Takes a few seconds
#Bootstrap means pretty much approach the sample means

#2019 estimated contribution to biomass-----------------
estimates$massContrib.mean<-100*estimates$germ.boot.mean*estimates$surv.boot.mean*estimates$mass.boot.mean

#Propagate errors
estimates$massContrib.95CIlo<-sqrt(((estimates$germ.boot.mean-estimates$germ.boot.95CIlo)/estimates$germ.boot.mean)^2
                                   +((estimates$surv.boot.mean-estimates$surv.boot.95CIlo)/estimates$surv.boot.mean)^2
                                   +((estimates$mass.boot.mean-estimates$mass.boot.95CIlo)/estimates$mass.boot.mean)^2
)
estimates$massContrib.95CIhi<-sqrt(((estimates$germ.boot.mean-estimates$germ.boot.95CIhi)/estimates$germ.boot.mean)^2
                                   +((estimates$surv.boot.mean-estimates$surv.boot.95CIhi)/estimates$surv.boot.mean)^2
                                   +((estimates$mass.boot.mean-estimates$mass.boot.95CIhi)/estimates$mass.boot.mean)^2
)

#Visualize 2019 projections-----------------
#Add labels back to data
estimates$Transplant<-substr(estimates$Trans.Don.Treat,1,4)
estimates$Donorspp<-substr(estimates$Trans.Don.Treat,6,9)
estimates$Treatment<-substr(estimates$Trans.Don.Treat,11,100)
#Delete final row from data since it has no treatment label
estimates<-estimates[1:60,]
#Make Intra/Inter label
estimates$Type[which(estimates$Donorspp!=estimates$Transplant)]<-"Inter"
estimates$Type[which(estimates$Donorspp==estimates$Transplant)]<-"Intra"
estimates$Type <- factor(estimates$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
estimates$Treatment <- factor(estimates$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization

#Make figure
ggplot(data=estimates, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
  facet_grid(~Transplant)+
  geom_pointrange(aes(ymin=massContrib.mean-massContrib.95CIlo, ymax=massContrib.95CIhi+massContrib.mean)
                  , position = position_dodge(width = 1))+  #All errors are too small to see
  scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
  theme_bw()

#Alternative visualization
ggplot(data=estimates, aes(x=Donorspp, y=massContrib.mean, color=Transplant, shape=Transplant))+
  facet_grid(~Treatment)+
  geom_pointrange(aes(ymin=massContrib.mean-massContrib.95CIlo, ymax=massContrib.95CIhi+massContrib.mean))+  #All errors are too small to see
  geom_line(aes(group=Transplant))+
  scale_color_manual(values=c("#E2D200", "#46ACC8", "#E58601", "#00A08A"))+
  xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
  theme_bw()
