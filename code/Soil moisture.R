################################################################
#Code to analyze soil moisture data during the experiment
################################################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(visreg)

#Data 180612 (still no consistent plant or spatial effects)---------------
#Data
all<-read.csv("./data/Transplant data 180612.csv")
dat180418<-subset(all,!is.na(VWC180418))
dat180531<-subset(all,!is.na(VWC180531))

#180418
hist(dat180418$VWC180418)#pretty normal
boxplot(dat180418$VWC180418~dat180418$Transect)#no differences
boxplot(dat180418$VWC180418~dat180418$DonorSpp)#POSE seems a bit low
boxplot(dat180418$VWC180418~dat180418$Rep)#quite a bit of variation
plot(dat180418$Position,dat180418$VWC180418,pch=dat180418$Transect,col=dat180418$Transect)#no super strong slope patterns
m1<-lm(VWC180418~DonorSpp*factor(Transect),data=dat180418)#plant effects
summary(m1)#nothing
m2<-lm(VWC180418~Position*factor(Transect),data=dat180418)#spatial effects
summary(m2)#position:transect 3 p=0.0449 (positive slope)

#180531
hist(dat180531$VWC180531)#pretty normal
boxplot(dat180531$VWC180531~dat180531$Transect)#no differences
boxplot(dat180531$VWC180531~dat180531$DonorSpp)#ARTR seems a bit low
boxplot(dat180531$VWC180531~dat180531$Rep)#not as much variation
plot(dat180531$Position,dat180531$VWC180531,pch=dat180531$Transect,col=dat180531$Transect)#no super strong slope patterns
m1<-lm(VWC180531~DonorSpp*factor(Transect),data=dat180531)#plant effects
summary(m1)#BARE has highest VWC (p=0.0509 from ARTR); HECO:Transect2 (p=0.04)
visreg(m1,"DonorSpp","Transect",type="conditional")
m2<-lm(VWC180531~Position*factor(Transect),data=dat180531)#spatial effects
summary(m2)#no relationship


#Data 171107 (no spatial or plant effects)---------------
#Data
all<-read.csv("./data/Transplant data 171107.csv")
dat<-subset(all,!is.na(VWC171110))

#Prelim examination/dataviz
hist(dat$VWC171110)#pretty normal
boxplot(dat$VWC171110~dat$Transect)#not too different
boxplot(dat$VWC171110~dat$DonorSpp)#HECO seems a bit low, ARTR largest range
boxplot(dat$VWC171110~dat$Rep)#not too different
plot(dat$Position,dat$VWC171110,pch=dat$Transect,col=dat$Transect)#no super strong slope patterns

#Stats
m1<-lm(VWC171110~DonorSpp*factor(Transect),data=dat)#plant effects
summary(m1)#nothing
m2<-lm(VWC171110~Position*factor(Transect),data=dat)#spatial effects
summary(m2)#nothing

# Data 190509 ------------------------------
dat<-read.csv("./data/Soil moisture 190509.csv")
dat$locality<-paste(dat$Transect,dat$Position,dat$Offset,dat$Side,sep=".")

#Prelim examination/dataviz
hist(dat$Moist_190509)#pretty normal
boxplot(dat$Moist_190509~dat$Transect)#not too different
boxplot(dat$Moist_190509~dat$DonorSpp)#HECO seems higher
boxplot(dat$Moist_190509~dat$Transplant)#Same
boxplot(dat$Moist_190509~dat$Treatment)#Same
boxplot(dat$Moist_190509~dat$Rep)#reps 1 and 2 noticeably higher
plot(dat$Position,dat$Moist_190509,pch=dat$Transect,col=dat$Transect)
#Maybe slope pattern for transect 1, definitely microsite patterns
boxplot(dat$Moist_190509~dat$locality)#lots of variation

#Stats
m1<-lm(Moist_190509~DonorSpp*Transplant*Treatment,data=dat)#Manipulation effects
Anova(m1,type = 3)#No sig effects 

#Include spatial random effects
library(lme4)
m2<-lmer(Moist_190509~DonorSpp*Transplant*Treatment
         +(1|locality),data=dat)#spatial effects
summary(m2)#locality accounts for >50% of total variance
Anova(m2,type = 3)#DonorSpp p=0.06

m3<-lmer(Moist_190509~DonorSpp*Transplant*Treatment
         +(1|Rep),data=dat)#spatial effects
summary(m3)#locality accounts for ~25% of total variance
Anova(m3,type = 3)#DonorSpp p=0.06
