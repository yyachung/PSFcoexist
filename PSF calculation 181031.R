setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(visreg)
library(car)
library(lme4)
library(ggpubr)

#Using Aboveground Biomass as a response

#Below is using area data-------------------------------------
#Data
all<-read.csv("Transplant data 180712.csv")

#Calculate area
all$Area<-all$Length180709/2*all$Width180709/2*pi #calculate as ellipse

all$levels<-paste(all$DonorSpp,all$Treatment,sep=".")

ARTR<-subset(all,Transplant=="ARTR")
POSE<-subset(all,Transplant=="POSE")
PSSP<-subset(all,Transplant=="PSSP")
HECO<-subset(all,Transplant=="HECO")

#ARTR Home vs. Away -------------------------------------------
PSF.ARTR<-data.frame(DonorSpp=character(),
                         Rep=integer(), 
                         Treatment=character(), 
                         lnRatio=double(),
                         stringsAsFactors=FALSE) 
for (i in 1:10){
  subdata<-subset(ARTR,Rep==i)
  home<-subset(subdata,DonorSpp=="ARTR")
  away<-subset(subdata,DonorSpp!="ARTR")
  for (j in 1:4){
    Donor<-unique(away$DonorSpp)[j]
    feedback<-log(home$Area/away$Area[ARTR$DonorSpp==Donor])#should be a vector of 3 values
    PSF.ARTR<-rbind(PSF.ARTR,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=as.character(home$Treatment),
                                  lnRatio=feedback[1:3]),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
boxplot(PSF.ARTR$lnRatio~PSF.ARTR$DonorSpp*PSF.ARTR$Treatment)
PSF.ARTR$levels<-paste(PSF.ARTR$Treatment,PSF.ARTR$DonorSpp,sep=".")
ggerrorplot(PSF.ARTR, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")

#Analysis
hist(PSF.ARTR$lnRatio)#may be a bit spiky
#Effects of Treatment and Donor Spp
ARTR1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.ARTR)
qqPlot(ARTR1)#fine
summary(ARTR1)#not sig diff from each other
drop1(ARTR1,~.,test="F")#no effects
#Differences from zero?
#Use means parameterization to figure out difference from zero
ARTR1m<-lm(lnRatio~levels-1,data=PSF.ARTR)
summary(ARTR1m)#only HECO exlusion p=0.066

#PSSP Home vs. Away -------------------------------------------
PSF.PSSP<-data.frame(DonorSpp=character(),
                     Rep=integer(), 
                     Treatment=character(), 
                     lnRatio=double(),
                     stringsAsFactors=FALSE) 
for (i in 1:10){
  subdata<-subset(PSSP,Rep==i)
  home<-subset(subdata,DonorSpp=="PSSP")
  away<-subset(subdata,DonorSpp!="PSSP")
  for (j in 1:4){
    Donor<-unique(away$DonorSpp)[j]
    feedback<-log(home$Area/away$Area[PSSP$DonorSpp==Donor])#should be a vector of 3 values
    PSF.PSSP<-rbind(PSF.PSSP,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=as.character(home$Treatment),
                                  lnRatio=feedback[1:3]),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
boxplot(PSF.PSSP$lnRatio~PSF.PSSP$DonorSpp*PSF.PSSP$Treatment)
PSF.PSSP$levels<-paste(PSF.PSSP$Treatment,PSF.PSSP$DonorSpp,sep=".")
ggerrorplot(PSF.PSSP, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")#really large error bars

#Analysis
hist(PSF.PSSP$lnRatio)#may be a bit spiky
#Effects of Treatment and Donor Spp
PSSP1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.PSSP)
qqPlot(PSSP1)#fine
summary(PSSP1)#not sig diff from each other
drop1(PSSP1,~.,test="F")#no effects
#Differences from zero?
#Use means parameterization to figure out difference from zero
PSSP1m<-lm(lnRatio~levels-1,data=PSF.PSSP)
summary(PSSP1m)#nothing

#POSE Home vs. Away [For some reason Donor==PSSP isn't working]-------------------------------------------
PSF.POSE<-data.frame(DonorSpp=character(),
                     Rep=integer(), 
                     Treatment=character(), 
                     lnRatio=double(),
                     stringsAsFactors=FALSE) 
for (i in 1:10){
  subdata<-subset(POSE,Rep==i)
  home<-subset(subdata,DonorSpp=="POSE")
  away<-subset(subdata,DonorSpp!="POSE")
  for (j in 1:4){
    Donor<-as.character(unique(away$DonorSpp)[j])
    feedback<-log(home$Area/away$Area[POSE$DonorSpp==Donor])#should be a vector of 3 values
    PSF.POSE<-rbind(PSF.POSE,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=as.character(home$Treatment),
                                  lnRatio=feedback[1:3]),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
boxplot(PSF.POSE$lnRatio~PSF.POSE$DonorSpp*PSF.POSE$Treatment)
PSF.POSE$levels<-paste(PSF.POSE$Treatment,PSF.POSE$DonorSpp,sep=".")
ggerrorplot(PSF.POSE, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")#really large error bars

#Analysis
hist(PSF.POSE$lnRatio)#may be a bit spiky
#Effects of Treatment and Donor Spp
POSE1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.POSE)
qqPlot(POSE1)#fine
summary(POSE1)#not sig diff from each other
drop1(POSE1,~.,test="F")#no effects
#Differences from zero?
#Use means parameterization to figure out difference from zero
POSE1m<-lm(lnRatio~levels-1,data=PSF.POSE)
summary(POSE1m)#nothing
