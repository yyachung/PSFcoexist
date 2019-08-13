#Examining relationships between shoot and root biomass#
#From 2018 field harvest data                          #
#AC 19 Jun 2019                                        #

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")

library(ggplot2)
library(car)

all<-read.csv("./data/Transplant data clean 190102.csv")

#Data visualization------------------------------------
#Everything together
ggplot(data=all)+
  geom_point(aes(x=Abv,y=Blw.all,color=Transplant,shape=Treatment))+
  geom_abline(slope=1,intercept=0)
#A few outliers with high shoot low root, but not any particular species or treatment

#ARTR only
ggplot(data=subset(all,Transplant=="ARTR"))+
  geom_point(aes(x=Abv,y=Blw.all,color=Transplant,shape=Treatment))+
  geom_abline(slope=1,intercept=0)

#HECO only
ggplot(data=subset(all,Transplant=="HECO"))+
  geom_point(aes(x=Abv,y=Blw.all,color=Transplant,shape=Treatment))+
  geom_abline(slope=1,intercept=0)

#POSE only
ggplot(data=subset(all,Transplant=="POSE"))+
  geom_point(aes(x=Abv,y=Blw.all,color=Transplant,shape=Treatment))+
  geom_abline(slope=1,intercept=0)

#PSSP only
ggplot(data=subset(all,Transplant=="PSSP"))+
  geom_point(aes(x=Abv,y=Blw.all,color=Transplant,shape=Treatment))+
  geom_abline(slope=1,intercept=0)

#Root:shoot
all$RSratio<-all$Blw.all/all$Abv
boxplot(all$RSratio~all$Treatment*all$Transplant)#for POSE and PSSP Feedback has the highest root:shoot

#Just root
boxplot(all$Blw.all~all$Treatment*all$Transplant)

#Prelim stats-----------------------------
clean<-subset(all,Treatment!="")

mR<-lm(log(Blw.all)~Abv*Treatment*Transplant,data=clean)
qqPlot(mR$resid)#better
drop1(mR,~.,test="F")#all sig

mRS<-lm(log(RSratio)~Treatment*Transplant,data=clean)
qqPlot(mRS$resid)#better
drop1(mRS,~.,test="F")#all sig

#Just ARTR
ARTR<-subset(clean,Transplant=="ARTR")
m.artr1<-lm(Blw.all~Treatment*DonorSpp,data=ARTR)
qqPlot(m.artr1$resid)
drop1(m.artr1,~.,test="F")#no sig effects

#Just HECO
HECO<-subset(clean,Transplant=="HECO")
m.heco1<-lm(Blw.all~Treatment*DonorSpp,data=HECO)
qqPlot(m.heco1$resid)#not great
drop1(m.heco1,~.,test="F")#treatment effects

#Just PSSP
PSSP<-subset(clean,Transplant=="PSSP")
m.pssp1<-lm(Blw.all~Treatment*DonorSpp,data=PSSP)
qqPlot(m.pssp1$resid)#fine
drop1(m.pssp1,~.,test="F")#treatment effects

#Just POSE
POSE<-subset(clean,Transplant=="POSE")
m.pose1<-lm(Blw.all~Treatment*DonorSpp,data=POSE)
qqPlot(m.pose1$resid)#not great
drop1(m.pose1,~.,test="F")#treatment effects
