###########################################################
#This calculates ln-ratio PSF, Is, and conducts associated#
#analyses from the cleaned data                           #
###########################################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(visreg)
library(car)
library(lme4)
library(ggpubr)

#Using Aboveground Biomass as a response---------------------------
all<-read.csv("./data/Transplant data clean 190102.csv")
ARTR<-subset(all,Transplant=="ARTR")
POSE<-subset(all,Transplant=="POSE")
PSSP<-subset(all,Transplant=="PSSP")
HECO<-subset(all,Transplant=="HECO")

control<-subset(all,Treatment=="Control")
exclusion<-subset(all,Treatment=="Exclusion")#only 199 bc of the messed up label
feedback<-subset(all,Treatment=="Feedback")#only 199 bc of the messed up label

#ARTR Home vs. Away [no sig donor/treatment effects, a few sig neg feedbacks]-------------------
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
    feedback.ratio<-log(home$Abv[home$Treatment=="Feedback"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Feedback"])
    exclusion.ratio<-log(home$Abv[home$Treatment=="Exclusion"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Exclusion"])
    control.ratio<-log(home$Abv[home$Treatment=="Control"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Control"])
    PSF.ARTR<-rbind(PSF.ARTR,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=c("Feedback","Exclusion","Control"),
                                  lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
PSF.ARTR$levels<-paste(PSF.ARTR$Treatment,PSF.ARTR$DonorSpp,sep=".")
ggerrorplot(PSF.ARTR, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")#warning for removing 56 NA values (checks out)
#weirdly, "Control" looks like the least negative feedbacks

#Analysis
hist(PSF.ARTR$lnRatio)#fine
#Effects of Treatment and Donor Spp
ARTR1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.ARTR)
qqPlot(ARTR1)#fine
summary(ARTR1)#not sig diff from each other
drop1(ARTR1,~.,test="F")#I'm surprised, actually
#Differences from zero?
#Use means parameterization to figure out difference from zero
ARTR1m<-lm(lnRatio~levels-1,data=PSF.ARTR)
summary(ARTR1m)#Exclusion PSSP, Feedback HECO, and Feedback PSSP are significantly less than zero

#HECO Home vs. Away [no sig effects, one pos feedback]--------
PSF.HECO<-data.frame(DonorSpp=character(),
                     Rep=integer(), 
                     Treatment=character(), 
                     lnRatio=double(),
                     stringsAsFactors=FALSE) 
for (i in 1:10){
  subdata<-subset(HECO,Rep==i)
  home<-subset(subdata,DonorSpp=="HECO")
  away<-subset(subdata,DonorSpp!="HECO")
  for (j in 1:4){
    Donor<-unique(away$DonorSpp)[j]
    feedback.ratio<-log(home$Abv[home$Treatment=="Feedback"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Feedback"])
    exclusion.ratio<-log(home$Abv[home$Treatment=="Exclusion"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Exclusion"])
    control.ratio<-log(home$Abv[home$Treatment=="Control"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Control"])
    PSF.HECO<-rbind(PSF.HECO,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=c("Feedback","Exclusion","Control"),
                                  lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
PSF.HECO$levels<-paste(PSF.HECO$Treatment,PSF.HECO$DonorSpp,sep=".")
ggerrorplot(PSF.HECO, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")
#Woah weird. Mostly trending positive feedbacks

#Analysis
hist(PSF.HECO$lnRatio)#slight poisson
#Effects of Treatment and Donor Spp
HECO1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.HECO)
qqPlot(HECO1)#fine
summary(HECO1)#not sig diff from each other
drop1(HECO1,~.,test="F")
#Differences from zero?
#Use means parameterization to figure out difference from zero
HECO1m<-lm(lnRatio~levels-1,data=PSF.HECO)
summary(HECO1m)#Only Exclusion BARE shows up as sig pos

#POSE Home vs. Away [no sig effects, one neg feedback]--------
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
    Donor<-unique(away$DonorSpp)[j]
    feedback.ratio<-log(home$Abv[home$Treatment=="Feedback"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Feedback"])
    exclusion.ratio<-log(home$Abv[home$Treatment=="Exclusion"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Exclusion"])
    control.ratio<-log(home$Abv[home$Treatment=="Control"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Control"])
    PSF.POSE<-rbind(PSF.POSE,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=c("Feedback","Exclusion","Control"),
                                  lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
PSF.POSE$levels<-paste(PSF.POSE$Treatment,PSF.POSE$DonorSpp,sep=".")
ggerrorplot(PSF.POSE, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")
#Mostly negative feedbacks

#Analysis
hist(PSF.POSE$lnRatio)#fine
#Effects of Treatment and Donor Spp
POSE1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.POSE)
qqPlot(POSE1)#fine
summary(POSE1)#not sig diff from each other
drop1(POSE1,~.,test="F")
#Differences from zero?
#Use means parameterization to figure out difference from zero
POSE1m<-lm(lnRatio~levels-1,data=PSF.POSE)
summary(POSE1m)#Only Exclusion BARE shows up as sig neg

#PSSP Home vs. Away [no sig effects or feedback]--------
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
    feedback.ratio<-log(home$Abv[home$Treatment=="Feedback"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Feedback"])
    exclusion.ratio<-log(home$Abv[home$Treatment=="Exclusion"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Exclusion"])
    control.ratio<-log(home$Abv[home$Treatment=="Control"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Control"])
    PSF.PSSP<-rbind(PSF.PSSP,list(DonorSpp=rep(as.character(Donor),3),
                                  Rep=rep(i,3),
                                  Treatment=c("Feedback","Exclusion","Control"),
                                  lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                    stringsAsFactors=FALSE)
  }
  
}

#Visualization
PSF.PSSP$levels<-paste(PSF.PSSP$Treatment,PSF.PSSP$DonorSpp,sep=".")
ggerrorplot(PSF.PSSP, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")
#weird contrast between positive feedbacks in exclusion and negative feedbacks in Feedback; 
#Control has large SEs likely due to missing one rep (home plant died)

#Analysis
hist(PSF.PSSP$lnRatio)#fine
#Effects of Treatment and Donor Spp
PSSP1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.PSSP)
qqPlot(PSSP1)#fine
summary(PSSP1)#not sig diff from each other
drop1(PSSP1,~.,test="F")
#Differences from zero?
#Use means parameterization to figure out difference from zero
PSSP1m<-lm(lnRatio~levels-1,data=PSF.PSSP)
summary(PSSP1m)#Actually no sig diffs from zero



#Is (only ARTR-PSSP "feedback" sig neg)-----------------------------------
PSF.ARTR$Transplant<-rep("ARTR",120)
PSF.HECO$Transplant<-rep("HECO",120)
PSF.POSE$Transplant<-rep("POSE",120)
PSF.PSSP$Transplant<-rep("PSSP",120)
PSF.all<-rbind(PSF.ARTR,PSF.HECO,PSF.POSE,PSF.PSSP)
PSF.all<-subset(PSF.all,DonorSpp!="BARE")#we only care about plants
#Make spp pair label
PSF.all$Pair[(PSF.all$DonorSpp=="ARTR"&PSF.all$Transplant=="POSE")|(PSF.all$Transplant=="ARTR"&PSF.all$DonorSpp=="POSE")]<-"ARTRPOSE"
PSF.all$Pair[(PSF.all$DonorSpp=="ARTR"&PSF.all$Transplant=="HECO")|(PSF.all$Transplant=="ARTR"&PSF.all$DonorSpp=="HECO")]<-"ARTRHECO"
PSF.all$Pair[(PSF.all$DonorSpp=="ARTR"&PSF.all$Transplant=="PSSP")|(PSF.all$Transplant=="ARTR"&PSF.all$DonorSpp=="PSSP")]<-"ARTRPSSP"
PSF.all$Pair[(PSF.all$DonorSpp=="HECO"&PSF.all$Transplant=="PSSP")|(PSF.all$Transplant=="HECO"&PSF.all$DonorSpp=="PSSP")]<-"HECOPSSP"
PSF.all$Pair[(PSF.all$DonorSpp=="HECO"&PSF.all$Transplant=="POSE")|(PSF.all$Transplant=="HECO"&PSF.all$DonorSpp=="POSE")]<-"HECOPOSE"
PSF.all$Pair[(PSF.all$DonorSpp=="PSSP"&PSF.all$Transplant=="POSE")|(PSF.all$Transplant=="PSSP"&PSF.all$DonorSpp=="POSE")]<-"PSSPPOSE"

#Calculate Is for all species pairs
Is<-data.frame(Pair=character(),
                     Rep=integer(), 
                     Treatment=character(), 
                     lnRatio=double(),
                     stringsAsFactors=FALSE) 
for (i in 1:10){
  subdata<-subset(PSF.all,Rep==i)
  for (j in 1:6){
    pair<-unique(subdata$Pair)[j]
    pairdata<-subset(subdata,Pair==pair)
    feedback.ratio<-sum(pairdata$lnRatio[pairdata$Treatment=="Feedback"])
    exclusion.ratio<-sum(pairdata$lnRatio[pairdata$Treatment=="Exclusion"])
    control.ratio<-sum(pairdata$lnRatio[pairdata$Treatment=="Control"])
    Is<-rbind(Is,list(Pair=rep(as.character(pair),3),
                                  Rep=rep(i,3),
                                  Treatment=c("Feedback","Exclusion","Control"),
                                  lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                    stringsAsFactors=FALSE)
  }

}

#Visualize
ggerrorplot(Is, x = "Pair", y = "lnRatio",
            facet.by = "Treatment")
#Use means parameterization to figure out difference from zero
Is$levels<-paste(Is$Treatment,Is$Pair,sep=".")
m.IS<-lm(lnRatio~levels-1,data=Is)
summary(m.IS)#Only Feedback ARTR-PSSP sig neg


