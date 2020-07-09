###########################################################
#This calculates ln-ratio PSF, Is, and conducts associated#
#analyses from the cleaned data                           #
###########################################################

#setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
setwd("C:/Users/yyach/Desktop/Box Sync/PSFcoexist/R/PSFcoexist")
library(emmeans)
library(car)
library(lme4)
library(ggpubr)

#Using Aboveground Biomass as a response---------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Clean 2018 data
all2018<-all2018[-which(all2018$Treatment==""),]
all2018$Treatment <- factor(all2018$Treatment)#drop blank factor level

#Subset data
ARTR2018<-subset(all2018,Transplant=="ARTR")
POSE2018<-subset(all2018,Transplant=="POSE")
PSSP2018<-subset(all2018,Transplant=="PSSP")
HECO2018<-subset(all2018,Transplant=="HECO")

ARTR2019<-subset(all2019,Transplant=="ARTR")
POSE2019<-subset(all2019,Transplant=="POSE")
PSSP2019<-subset(all2019,Transplant=="PSSP")
HECO2019<-subset(all2019,Transplant=="HECO")

# control<-subset(all,Treatment=="Control")
# exclusion<-subset(all,Treatment=="Exclusion")#only 199 in 2018 bc of the messed up label
# feedback<-subset(all,Treatment=="Feedback")#only 199 in 2018 bc of the messed up label

#Calculate home vs. away PSF ln ratios for all species--------
source("./code/func_PSFcalc.R")
PSF.ARTR18<-PSFcalc(ARTR2018,"ARTR")
PSF.HECO18<-PSFcalc(HECO2018,"HECO")
PSF.POSE18<-PSFcalc(POSE2018,"POSE")
PSF.PSSP18<-PSFcalc(PSSP2018,"PSSP")
PSF.ARTR19<-PSFcalc(ARTR2019,"ARTR")
PSF.HECO19<-PSFcalc(HECO2019,"HECO")
PSF.POSE19<-PSFcalc(POSE2019,"POSE")
PSF.PSSP19<-PSFcalc(PSSP2019,"PSSP")

#Prelim viz
par(mfrow = c(4, 2))
hist(PSF.ARTR18$lnRatio) 
hist(PSF.HECO18$lnRatio) 
hist(PSF.POSE18$lnRatio)  
hist(PSF.PSSP18$lnRatio) 
hist(PSF.ARTR19$lnRatio) 
hist(PSF.HECO19$lnRatio) 
hist(PSF.POSE19$lnRatio)  
hist(PSF.PSSP19$lnRatio) 
dev.off()

#ARTR Home vs. Away [no sig donor/treatment effects, a few sig neg feedbacks]-------------------
#Check sample sizes
table(na.omit(PSF.ARTR18)$DonorSpp,na.omit(PSF.ARTR18)$Treatment)#range 3 to 9
table(na.omit(PSF.ARTR19)$DonorSpp,na.omit(PSF.ARTR19)$Treatment)#range 5 to 9

#Visualization
PSF.ARTR18$levels<-paste(PSF.ARTR18$Treatment,PSF.ARTR18$DonorSpp,sep=".")
ggerrorplot(PSF.ARTR18, x = "DonorSpp", y = "lnRatio",
            facet.by = "Treatment")#warning for removing 56 NA values (checks out)
#weirdly, "Control" looks like the least negative feedbacks

#Analysis
hist(PSF.ARTR18$lnRatio)#fine
#Effects of Treatment and Donor Spp
ARTR1<-lm(lnRatio~Treatment*DonorSpp,data=PSF.ARTR18)
qqPlot(ARTR1)#fine
summary(ARTR1)#not sig diff from each other
drop1(ARTR1,~.,test="F")#I'm surprised, actually
#Differences from zero?
#Use means parameterization to figure out difference from zero
ARTR1m<-lm(lnRatio~levels-1,data=PSF.ARTR18)
summary(ARTR1m)#Exclusion PSSP, Feedback HECO, and Feedback PSSP are significantly less than zero
#None of this has been corrected for multiple comparisons, but maybe doesn't need to?

#Alternative analysis to take care of unbalanced design
hist(log(ARTR2018$Abv))#well that's nice
m.ARTR18 <- lmer(log(Abv)~DonorSpp*Treatment + (1|Rep),data=ARTR2018)
#boundary (singular) fit warning, probably due to some low sample sizes
summary(m.ARTR18)#Also random effect variance is zero
#Get contrasts
emDonor.ARTR18<-emmeans(m.ARTR18,~DonorSpp|Treatment)
pairs(emDonor.ARTR18)#includes too many comparisons
#subset and readjust for only the contrasts that are relevant for PSF
summary(contrast(emDonor.ARTR18, "pairwise")[c(1:4,11:14,21:24)],adjust="mvt")
#Feedback HECO, and Feedback PSSP are significantly less than zero

#HECO Home vs. Away [no sig effects, one pos feedback]--------
#Check sample sizes
table(na.omit(PSF.HECO18)$DonorSpp,na.omit(PSF.HECO18)$Treatment)#range 6 to 10
table(na.omit(PSF.HECO19)$DonorSpp,na.omit(PSF.HECO19)$Treatment)#range 6 to 9

#Visualization
PSF.HECO19$levels<-paste(PSF.HECO19$Treatment,PSF.HECO19$DonorSpp,sep=".")
ggerrorplot(PSF.HECO19, x = "DonorSpp", y = "lnRatio",
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

#POSE Home vs. Away [no sig effects, one neg feedback]-----------------------------
table(na.omit(PSF.POSE18)$DonorSpp,na.omit(PSF.POSE18)$Treatment)#all 10
table(na.omit(PSF.POSE19)$DonorSpp,na.omit(PSF.POSE19)$Treatment)#one 9 others 10

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

#PSSP Home vs. Away [no sig effects or feedback]-----------------------------------
table(na.omit(PSF.PSSP18)$DonorSpp,na.omit(PSF.PSSP18)$Treatment)#9 and 10
table(na.omit(PSF.PSSP19)$DonorSpp,na.omit(PSF.PSSP19)$Treatment)#all 10

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
#Combine all 2018 PSF data
PSF.ARTR18$Transplant<-rep("ARTR",120)
PSF.HECO18$Transplant<-rep("HECO",120)
PSF.POSE18$Transplant<-rep("POSE",120)
PSF.PSSP18$Transplant<-rep("PSSP",120)
PSF.all18<-rbind(PSF.ARTR18,PSF.HECO18,PSF.POSE18,PSF.PSSP18)
PSF.all18<-subset(PSF.all18,DonorSpp!="BARE")#we only care about plants
#Make spp pair label
PSF.all18$Pair[(PSF.all18$DonorSpp=="ARTR"&PSF.all18$Transplant=="POSE")|(PSF.all18$Transplant=="ARTR"&PSF.all18$DonorSpp=="POSE")]<-"ARTRPOSE"
PSF.all18$Pair[(PSF.all18$DonorSpp=="ARTR"&PSF.all18$Transplant=="HECO")|(PSF.all18$Transplant=="ARTR"&PSF.all18$DonorSpp=="HECO")]<-"ARTRHECO"
PSF.all18$Pair[(PSF.all18$DonorSpp=="ARTR"&PSF.all18$Transplant=="PSSP")|(PSF.all18$Transplant=="ARTR"&PSF.all18$DonorSpp=="PSSP")]<-"ARTRPSSP"
PSF.all18$Pair[(PSF.all18$DonorSpp=="HECO"&PSF.all18$Transplant=="PSSP")|(PSF.all18$Transplant=="HECO"&PSF.all18$DonorSpp=="PSSP")]<-"HECOPSSP"
PSF.all18$Pair[(PSF.all18$DonorSpp=="HECO"&PSF.all18$Transplant=="POSE")|(PSF.all18$Transplant=="HECO"&PSF.all18$DonorSpp=="POSE")]<-"HECOPOSE"
PSF.all18$Pair[(PSF.all18$DonorSpp=="PSSP"&PSF.all18$Transplant=="POSE")|(PSF.all18$Transplant=="PSSP"&PSF.all18$DonorSpp=="POSE")]<-"PSSPPOSE"

#Combine all 2019 PSF data
PSF.ARTR19$Transplant<-rep("ARTR",120)
PSF.HECO19$Transplant<-rep("HECO",120)
PSF.POSE19$Transplant<-rep("POSE",120)
PSF.PSSP19$Transplant<-rep("PSSP",120)
PSF.all19<-rbind(PSF.ARTR19,PSF.HECO19,PSF.POSE19,PSF.PSSP19)
PSF.all19<-subset(PSF.all19,DonorSpp!="BARE")#we only care about plants
#Make spp pair label
PSF.all19$Pair[(PSF.all19$DonorSpp=="ARTR"&PSF.all19$Transplant=="POSE")|(PSF.all19$Transplant=="ARTR"&PSF.all19$DonorSpp=="POSE")]<-"ARTRPOSE"
PSF.all19$Pair[(PSF.all19$DonorSpp=="ARTR"&PSF.all19$Transplant=="HECO")|(PSF.all19$Transplant=="ARTR"&PSF.all19$DonorSpp=="HECO")]<-"ARTRHECO"
PSF.all19$Pair[(PSF.all19$DonorSpp=="ARTR"&PSF.all19$Transplant=="PSSP")|(PSF.all19$Transplant=="ARTR"&PSF.all19$DonorSpp=="PSSP")]<-"ARTRPSSP"
PSF.all19$Pair[(PSF.all19$DonorSpp=="HECO"&PSF.all19$Transplant=="PSSP")|(PSF.all19$Transplant=="HECO"&PSF.all19$DonorSpp=="PSSP")]<-"HECOPSSP"
PSF.all19$Pair[(PSF.all19$DonorSpp=="HECO"&PSF.all19$Transplant=="POSE")|(PSF.all19$Transplant=="HECO"&PSF.all19$DonorSpp=="POSE")]<-"HECOPOSE"
PSF.all19$Pair[(PSF.all19$DonorSpp=="PSSP"&PSF.all19$Transplant=="POSE")|(PSF.all19$Transplant=="PSSP"&PSF.all19$DonorSpp=="POSE")]<-"PSSPPOSE"

#Calculate Is for all species pairs
source("./code/func_IScalc.R")
IS2018<-IScalc(PSF.all18)
IS2019<-IScalc(PSF.all19)

#Check sample sizes
table(na.omit(IS2018)$Pair,na.omit(IS2018)$Treatment)#range 3 to 10
table(na.omit(IS2019)$Pair,na.omit(IS2019)$Treatment)#range 4 to 10

#Visualize
ggerrorplot(IS2018, x = "Pair", y = "lnRatio",
            orientation = "horizontal",
            facet.by = "Treatment")
ggerrorplot(IS2019, x = "Pair", y = "lnRatio",
            orientation = "horizontal",
            facet.by = "Treatment")

#Use means parameterization to figure out difference from zero
Is$levels<-paste(Is$Treatment,Is$Pair,sep=".")
m.IS<-lm(lnRatio~levels-1,data=Is)
summary(m.IS)#Only Feedback ARTR-PSSP sig neg

#Locality random effects using Abv Biomass----------------------------------------
all<-read.csv("./data/Transplant data clean 190102.csv")
raw<-read.csv("./data/Transplant data_labels corrected 190102.csv")
Locality<-paste(raw$Transect,raw$Position,raw$Offset,raw$Side,sep=".")
all$Locality<-Locality #Add locality label from the raw data

#Model 
library(lme4)
m1<-lmer(Abv~DonorSpp*Transplant*Treatment
         +(1|Locality),data=all)#fixed effect model rank deficient? Why?
#Problem:trying to estimate 5*4*3 +1 params with n=10?
length(colnames(model.matrix(m1)))#61
#Try fix
m2<-lmer(Abv~DonorSpp+Transplant+Treatment
         +(1|Locality),data=all)#This is fine apparently
length(colnames(model.matrix(m2)))#11
summary(m2)# 95 unique localities <10% total variance
Anova(m2,type=3)# Transplant and treatment are very significant

#Try directly looking at locality as factor?
m3<-lm(Abv~factor(Locality),data=all)
Anova(m3,type=3)#Locality not significant
