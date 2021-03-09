################################################################
#Code to analyze soil moisture data during the experiment
################################################################

#setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(emmeans)
library(car)

#Data 180612 (only small spatial effect in Apr)---------------
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

m1<-lm(VWC180418~DonorSpp+Position*factor(Transect),data=dat180418)
qqPlot(m1)#good
summary(m1)#POSE may be lower than ARTR; #Transect 3 has potential spatial effect
Anova(m1)#position:transect only based on ANOVA hypotheses

#180531
hist(dat180531$VWC180531)#pretty normal
boxplot(dat180531$VWC180531~dat180531$Transect)#no differences
boxplot(dat180531$VWC180531~dat180531$DonorSpp)#ARTR seems a bit low
boxplot(dat180531$VWC180531~dat180531$Rep)#not as much variation
plot(dat180531$Position,dat180531$VWC180531,pch=dat180531$Transect,col=dat180531$Transect)#no super strong slope patterns

m2<-lm(VWC180531~DonorSpp+Position*factor(Transect),data=dat180531)
qqPlot(m2)#good
summary(m2)#BARE may be higher than ARTR
Anova(m2)#no differences


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
m3<-lm(VWC171110~DonorSpp+Position*factor(Transect),data=dat)
qqPlot(m3)#good
summary(m3)#no effects
Anova(m3)#no differences


# Data 190605 ------------------------------
dat<-read.csv("./data/Soil moisture 190605.csv")
dat$locality<-paste(dat$Transect,dat$Position,dat$Offset,dat$Side,sep=".")

#VWC 190509
#Prelim examination/dataviz
hist(dat$Moist_190509)#pretty normal
boxplot(dat$Moist_190509~dat$Transect)#not too different
boxplot(dat$Moist_190509~dat$DonorSpp)#HECO seems higher
boxplot(dat$Moist_190509~dat$Transplant)#Same
boxplot(dat$Moist_190509~dat$Treatment)#Same
boxplot(dat$Moist_190509~dat$Rep)#reps 1 and 2 noticeably higher
plot(dat$Position,dat$Moist_190509,pch=dat$Transect,col=dat$Transect)
boxplot(dat$Moist_190509~dat$locality)#lots of variation
#Stats
m4<-lm(Moist_190509~DonorSpp+Position*factor(Transect),data=dat)
qqPlot(m4)#good
summary(m4)#most effects are significant, negative position slope for Transect 1
Anova(m4)#all effects significant
pairs(emmeans(m4,~DonorSpp))#ARTR=POSE=BARE<PSSP<HECO
pairs(emmeans(m4,~factor(Transect)*Position))#2=3<1

#VWC 190605
#Prelim examination/dataviz
hist(dat$Moist_190605)#pretty normal
boxplot(dat$Moist_190605~dat$Transect)#not too different
boxplot(dat$Moist_190605~dat$DonorSpp)#not too different
boxplot(dat$Moist_190605~dat$Transplant)#Same
boxplot(dat$Moist_190605~dat$Treatment)#Same
boxplot(dat$Moist_190605~dat$Rep)#reps 10 and 2 noticeably higher
plot(dat$Position,dat$Moist_190605,pch=dat$Transect,col=dat$Transect)
boxplot(dat$Moist_190605~dat$locality)#lots of variation
#Stats
m5<-lm(Moist_190605~DonorSpp+Position*factor(Transect),data=dat)
qqPlot(m5)#good
summary(m5)#a few significant effects in each category, negative slopes for Transect 1 and 2, positive for 3
Anova(m5)#Donor and position:transect sig
pairs(emmeans(m5,~DonorSpp))#POSE=PSSP=ARTR=BARE=HECO: a,ab,b,bc,c
pairs(emmeans(m5,~factor(Transect)*Position))#2<3<1

#Include spatial random effects (not used)
library(lme4)
m2<-lmer(Moist_190509~DonorSpp*Transplant*Treatment
         +(1|locality),data=dat)#spatial effects
summary(m2)#locality accounts for >50% of total variance
Anova(m2,type = 3)#DonorSpp p=0.06

m3<-lmer(Moist_190509~DonorSpp*Transplant*Treatment
         +(1|Rep),data=dat)#spatial effects
summary(m3)#locality accounts for ~25% of total variance
Anova(m3,type = 3)#DonorSpp p=0.06
