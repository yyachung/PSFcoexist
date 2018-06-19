setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(visreg)
library(car)
library(lme4)

#Data
all<-read.csv("Transplant data 180612.csv")
ARTR<-subset(all,Transplant=="ARTR")
POSE<-subset(all,Transplant=="POSE")
PSSP<-subset(all,Transplant=="PSSP")
HECO<-subset(all,Transplant=="HECO")

#180514 survival [no treatment effects]---------------------
#Replace NAs with ones
all$Surv180514[is.na(all$Surv180514)] <- 1

#Data visualization
hist(all$Surv180514)
aggregate(all$Surv180514,by=list(all$Treatment,all$DonorSpp,all$Transplant),FUN="mean")
#very little variation in survival except for ARTR

#Global Model
s1<-glm(Surv180514~Transplant*DonorSpp*Treatment,data=all,family=binomial)
hist(resid(s1))
summary(s1)#Three way interactions with so many levels difficult to interpret. Nothing pops out.

#ARTR only
s2<-glm(Surv180514~DonorSpp*Treatment,data=ARTR,family=binomial)
hist(resid(s2))
summary(s2)

#HECO only
s3<-glm(Surv180514~DonorSpp*Treatment,data=HECO,family=binomial)
hist(resid(s3))
summary(s3)

#POSE only
s4<-glm(Surv180514~DonorSpp*Treatment,data=POSE,family=binomial)#did not converge
#Probably all survived

#PSSP only
s5<-glm(Surv180514~DonorSpp*Treatment,data=PSSP,family=binomial)
hist(resid(s5))#way skewed
summary(s5)

#180514 height (all SPP) [height Excl>Feed>Ctrl for all spp; Donor effect for all except HECO; no interactions]------------------------------------------
hist(all$Height180514)#right skewed
boxplot(all$Height180514~all$Transplant)#should analyze spp separately
#Species separate
hist(ARTR$Height180514)#not bad
hist(POSE$Height180514)#not bad
hist(PSSP$Height180514)#not bad
hist(HECO$Height180514)#not bad

#ARTR
h1<-lmer(Height180514~DonorSpp*Treatment+(1|Rep),data=ARTR)
plot(h1)#not bad
hist(resid(h1))#fine
summary(h1)#rep about 13% of total variance
Anova(h1,test.statistic = "F")#Donor (p=0.013) Treat (p<0.0001), interaction (p=0.2)
visreg(h1,"DonorSpp","Treatment")

#HECO
h2<-lmer(Height180514~DonorSpp*Treatment+(1|Rep),data=HECO)
plot(h2)#a bit funnel-y
hist(resid(h2))#fine
summary(h2)#rep <5% total variance
Anova(h2,test.statistic = "F")#Donor (p=0.3) Treat (p<0.0001), interaction (p=0.06)
visreg(h2,"DonorSpp","Treatment")#not sure where the interaction trend is coming from

#POSE
h3<-lmer(Height180514~DonorSpp*Treatment+(1|Rep),data=POSE)
plot(h3)#not bad
hist(resid(h3))#fine
summary(h3)#rep about <5% of total variance
Anova(h3,test.statistic = "F")#Donor (p<0.0001) Treat (p<0.0001), interaction (p=0.4)
visreg(h3,"DonorSpp","Treatment")
#does best on ARTR soil and worst on BARE for all treatments, but difference biggest in exclusion

#PSSP
h4<-lmer(Height180514~DonorSpp*Treatment+(1|Rep),data=PSSP)
plot(h4)#not bad
hist(resid(h4))#fine
summary(h4)#rep about 16% of total variance
Anova(h4,test.statistic = "F")#Donor (p=0.013) Treat (p<0.0001), interaction (p=0.6)
visreg(h4,"DonorSpp","Treatment")
#does best on ARTR soil all treatments

#180514 Tillers (grasses only) [treatment effects for POSE and PSSP; Donor effects for POSE]---------------------------
hist(all$Tillers180514)#right skewed
boxplot(all$Tillers180514~all$Transplant)#POSE has most tillers, and most variation
#Species separate
hist(POSE$Tillers180514)#right skewed...poisson for count data?
hist(PSSP$Tillers180514)#right skewed...poisson for count data?
hist(HECO$Tillers180514)#right skewed...poisson for count data?

#HECO
t1<-glmer(Tillers180514~DonorSpp*Treatment+(1|Rep),data=HECO,family="poisson")
plot(t1)#fine
hist(resid(t1))#fine; warnings for NaNs? weird because default=na.omit
summary(t1)#rep variance 0?
Anova(t1)#no sig effects
visreg(t1,"DonorSpp","Treatment")#lots of NaN-based warnings...

#POSE
t2<-glmer(Tillers180514~DonorSpp*Treatment+(1|Rep),data=POSE,family="poisson")
plot(t2)#fine
hist(resid(t2))#fine
summary(t2)#Sig effects: TrtmntExcl, TrtmntFdbck
Anova(t2)#DonorSpp p=0.0015, Treatment p<0.0001
visreg(t2,"DonorSpp","Treatment")#DonorSPP effect comes from fewer tillers in ARTR soil than all others

#PSSP
t3<-glmer(Tillers180514~DonorSpp*Treatment+(1|Rep),data=PSSP,family="poisson")
plot(t3)#fine
hist(resid(t3))#fine; NaN warnings
summary(t3)#Trend for TrtmntExcl vs. control p=0.087
Anova(t3)#Treatment p=0.0009
visreg(t3,"DonorSpp","Treatment")#NaN warnings
#Treatment effect comes from fewer tillers in Control vs. others

#180606 Canopy (ARTR only) [Treatment effects]-------------------------
ARTR$CanopyArea<-ARTR$CanopyL180606/2*ARTR$CanopyW180606/2*pi #calculate as ellipse
hist(ARTR$CanopyArea)#right-skewed
hist(log(ARTR$CanopyArea))#better

a1<-lmer(log(CanopyArea)~DonorSpp*Treatment+(1|Rep),data=ARTR)
plot(a1)#fine
hist(resid(a1))#fine
summary(a1)#rep <0.1% of variance
Anova(a1,test="F")#Donor p=0.25, treatment p<0.0001, interaction p=0.9
visreg(a1,"DonorSpp","Treatment")
