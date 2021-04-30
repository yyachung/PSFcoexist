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

#Input data---------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Remove biomass outlier
all2019$Abv[which(all2019$Abv<0.0002)]<-NA

#Subset data
ARTR2018<-subset(all2018,Transplant=="ARTR")
POSE2018<-subset(all2018,Transplant=="POSE")
PSSP2018<-subset(all2018,Transplant=="PSSP")
HECO2018<-subset(all2018,Transplant=="HECO")

ARTR2019<-subset(all2019,Transplant=="ARTR")
POSE2019<-subset(all2019,Transplant=="POSE")
PSSP2019<-subset(all2019,Transplant=="PSSP")
HECO2019<-subset(all2019,Transplant=="HECO")

#Calculate home vs. away PSF ln ratios for all species--------
source("./code/func_PSFcalc.R")

PSF.ARTR18<-PSFcalc(ARTR2018,"ARTR","Germ","Abv")
PSF.HECO18<-PSFcalc(HECO2018,"HECO","Germ","Abv")
PSF.POSE18<-PSFcalc(POSE2018,"POSE","Germ","Abv")
PSF.PSSP18<-PSFcalc(PSSP2018,"PSSP","Germ","Abv")
PSF.ARTR19<-PSFcalc(ARTR2019,"ARTR","Germ","Abv")
PSF.HECO19<-PSFcalc(HECO2019,"HECO","Germ","Abv")
PSF.POSE19<-PSFcalc(POSE2019,"POSE","Germ","Abv")
PSF.PSSP19<-PSFcalc(PSSP2019,"PSSP","Germ","Abv")

#Prelim viz
par(mfrow = c(4, 2))
hist(PSF.ARTR18$lnRatio.germ) 
hist(PSF.HECO18$lnRatio.germ) 
hist(PSF.POSE18$lnRatio.germ)  
hist(PSF.PSSP18$lnRatio.germ) 
hist(PSF.ARTR19$lnRatio.germ) 
hist(PSF.HECO19$lnRatio.germ) 
hist(PSF.POSE19$lnRatio.germ)  
hist(PSF.PSSP19$lnRatio.germ) 
dev.off()

#ARTR Home vs. Away [no sig donor/treatment effects, a few sig neg feedbacks]-------------------
PSF.ARTR<-rbind.data.frame(PSF.ARTR18,PSF.ARTR19)
PSF.ARTR$Year<-rep(c(2018,2019),each=120)
PSF.ARTR$levels<-paste(PSF.ARTR$Treatment,PSF.ARTR$DonorSpp,PSF.ARTR$Year,sep=".")

#ABOVEGROUND BIOMASS
#Check sample sizes 
table(na.omit(PSF.ARTR18)$DonorSpp,na.omit(PSF.ARTR18)$Treatment)#range 3 to 9
table(na.omit(PSF.ARTR19)$DonorSpp,na.omit(PSF.ARTR19)$Treatment)#range 5 to 9

#ANOVA both years together
ARTR.m1<-lm(lnRatio.mass~Treatment*DonorSpp*factor(Year),data=PSF.ARTR)
qqPlot(ARTR.m1)#fine
Anova(ARTR.m1)#treatment:year p=0.08

#Significant difference from zero
ARTR.m2<-lm(lnRatio.mass~levels-1,data=PSF.ARTR)
summary(ARTR.m2)
#2018: Feedback HECO, and Feedback PSSP are significantly less than zero
#2019: Control PSSP, Exclusion PSSP are significantly less than zero
#P<0.1: 2019 Control POSE, 2019 Exclusion BARE, 2018 Exclusion PSSP 

#Visualization
# #need to make a damn dataframe with all the annotation params
# ann_text<-data.frame(DonorSpp=factor("HECO","PSSP","PSSP","PSSP",levels=unique(PSF.ARTR$DonorSpp)),
#                      lmRatio.mass=c(-0.5,-0.3,-0.5,-0.2),
#                      Treatment=factor("Feedback","Feedback","Control","Exclusion",levels=c("Control","Exclusion","Feedback")),
#                      Year=c("2018","2018","2019","2019"),
#                      lab=rep("*",4))
p.mARTR<-ggerrorplot(PSF.ARTR, x = "DonorSpp", y = "lnRatio.mass",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="ARTR feedback log ratio")+
    #  geom_text(data = ann_text,mapping = aes(x = DonorSpp, y = lnRatio.mass, label = lab))+
       geom_hline(yintercept=0,linetype=2)
  

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

#GERMINATION
#ANOVA both years together
ARTR.g1<-lm(lnRatio.germ~Treatment*DonorSpp*factor(Year),data=PSF.ARTR)
qqPlot(ARTR.g1)#fine
Anova(ARTR.g1)#nothing

#Significant difference from zero
ARTR.g2<-lm(lnRatio.germ~levels-1,data=PSF.ARTR)
summary(ARTR.g2)#nothing

#Visualization
p.gARTR<-ggerrorplot(PSF.ARTR, x = "DonorSpp", y = "lnRatio.germ",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="ARTR feedback log ratio")+
         geom_hline(yintercept=0,linetype=2)


#HECO Home vs. Away [no sig effects, one pos feedback]--------
PSF.HECO<-rbind.data.frame(PSF.HECO18,PSF.HECO19)
PSF.HECO$Year<-rep(c(2018,2019),each=120)
PSF.HECO$levels<-paste(PSF.HECO$Treatment,PSF.HECO$DonorSpp,PSF.HECO$Year,sep=".")

#ABOVEGROUND BIOMASS
#Check sample sizes
table(na.omit(PSF.HECO18)$DonorSpp,na.omit(PSF.HECO18)$Treatment)#range 6 to 10
table(na.omit(PSF.HECO19)$DonorSpp,na.omit(PSF.HECO19)$Treatment)#range 6 to 9

#ANOVA both years together
HECO.m1<-lm(lnRatio.mass~Treatment*DonorSpp*factor(Year),data=PSF.HECO)
qqPlot(HECO.m1)#fine
Anova(HECO.m1)#treatment:year p=0.10

#Significant difference from zero
HECO.m2<-lm(lnRatio.mass~levels-1,data=PSF.HECO)
summary(HECO.m2)
#2018: Exclusion BARE significantly positive
#2019: Feedback ARTR significantly positive
#P<0.1: 2018 Exclusion ARTR 

#Visualization
p.mHECO<-ggerrorplot(PSF.HECO, x = "DonorSpp", y = "lnRatio.mass",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="HECO feedback log ratio")+
  geom_hline(yintercept=0,linetype=2)
            

#GERMINATION
#ANOVA both years together
HECO.g1<-lm(lnRatio.germ~Treatment*DonorSpp*factor(Year),data=PSF.HECO)
qqPlot(HECO.g1)#fine
Anova(HECO.g1)#Treatment p=0.0498 , year p=0.00967 
#emmeans(HECO.g1,~Treatment) #exclusion lower
#emmeans(HECO.g1,~Year) #2019 higher

#Significant difference from zero
HECO.g2<-lm(lnRatio.germ~levels-1,data=PSF.HECO)
summary(HECO.g2)
#2019 Control ARTR, Control BARE sig positive
#P<0.1: 2019 Control PSSP, 2019 Feedback BARE positive

#Visualization
p.gHECO<-ggerrorplot(PSF.HECO, x = "DonorSpp", y = "lnRatio.germ",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="HECO feedback log ratio")+
  geom_hline(yintercept=0,linetype=2)

#POSE Home vs. Away [no sig effects, one neg feedback]-----------------------------
PSF.POSE<-rbind.data.frame(PSF.POSE18,PSF.POSE19)
PSF.POSE$Year<-rep(c(2018,2019),each=120)
PSF.POSE$levels<-paste(PSF.POSE$Treatment,PSF.POSE$DonorSpp,PSF.POSE$Year,sep=".")

#ABOVEGROUND BIOMASS
table(na.omit(PSF.POSE18)$DonorSpp,na.omit(PSF.POSE18)$Treatment)#all 10
table(na.omit(PSF.POSE19)$DonorSpp,na.omit(PSF.POSE19)$Treatment)#one 8 others 10

#ANOVA both years together
POSE.m1<-lm(lnRatio.mass~Treatment*DonorSpp*factor(Year),data=PSF.POSE)
qqPlot(POSE.m1)#fine
Anova(POSE.m1)#nothing

#Significant difference from zero
POSE.m2<-lm(lnRatio.mass~levels-1,data=PSF.POSE)
summary(POSE.m2)
#2018: Exclusion BARE significantly negative
#2019:Control BARE, Control HECO, Exclusion HECO, Feedback HECO, Feedback PSSP negative

#Visualization
p.mPOSE_out<-ggerrorplot(PSF.POSE, x = "DonorSpp", y = "lnRatio.mass",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="POSE feedback log ratio")+
  geom_hline(yintercept=0,linetype=2)

#GERMINATION
#ANOVA both years together
POSE.g1<-lm(lnRatio.germ~Treatment*DonorSpp*factor(Year),data=PSF.POSE)
qqPlot(POSE.g1)#fine
Anova(POSE.g1)#nothing 

#Significant difference from zero
POSE.g2<-lm(lnRatio.germ~levels-1,data=PSF.POSE)
summary(POSE.g2)
#2018 Control BARE sig positive

#Visualization
p.gPOSE<-ggerrorplot(PSF.POSE, x = "DonorSpp", y = "lnRatio.germ",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="POSE feedback log ratio")+
  geom_hline(yintercept=0,linetype=2)

#PSSP Home vs. Away [no sig effects or feedback]-----------------------------------
PSF.PSSP<-rbind.data.frame(PSF.PSSP18,PSF.PSSP19)
PSF.PSSP$Year<-rep(c(2018,2019),each=120)
PSF.PSSP$levels<-paste(PSF.PSSP$Treatment,PSF.PSSP$DonorSpp,PSF.PSSP$Year,sep=".")

#ABOVEGROUND BIOMASS
table(na.omit(PSF.PSSP18)$DonorSpp,na.omit(PSF.PSSP18)$Treatment)#9 and 10
table(na.omit(PSF.PSSP19)$DonorSpp,na.omit(PSF.PSSP19)$Treatment)#all 10

#ANOVA both years together
PSSP.m1<-lm(lnRatio.mass~Treatment*DonorSpp*factor(Year),data=PSF.PSSP)
qqPlot(PSSP.m1)#fine
Anova(PSSP.m1)#treatment p=0.0005
emmeans(PSSP.m1,~Treatment)#Exclusion positive whereas Control and Feedback negative

#Significant difference from zero
PSSP.m2<-lm(lnRatio.mass~levels-1,data=PSF.PSSP)
summary(PSSP.m2)
#2019: Control BARE, Control POSE, Feedback HECO significantly negative
#P<0.1: 2018 Feedback ARTR, 2018 Feedback HECO 

#Visualization
p.mPSSP<-ggerrorplot(PSF.PSSP, x = "DonorSpp", y = "lnRatio.mass",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="PSSP feedback log ratio")+
  geom_hline(yintercept=0,linetype=2)

#GERMINATION
#ANOVA both years together
PSSP.g1<-lm(lnRatio.germ~Treatment*DonorSpp*factor(Year),data=PSF.PSSP)
qqPlot(PSSP.g1)#fine
Anova(PSSP.g1)#Treatment P=0.056 Year P=0.0001 
#emmeans(PSSP.g1,~Treatment) #Feedback positive (driven by 2018)
#emmeans(PSSP.g1,~Year) #2018 positive, 2019 negative

#Significant difference from zero
PSSP.g2<-lm(lnRatio.germ~levels-1,data=PSF.PSSP)
summary(PSSP.g2)
#2018: Control BARE, Feedback BARE, Feedback HECO, Feedback POSE sig positive
#P<0.1: 2018 BARE exclusion pos

#Visualization
p.gPSSP<-ggerrorplot(PSF.PSSP, x = "DonorSpp", y = "lnRatio.germ",
            facet.by = c("Year","Treatment"),
            xlab="",
            ylab="PSSP feedback log ratio")+
  geom_hline(yintercept=0,linetype=2)

#Combine pairwise lnRatio plots into giant multipanel plots---------------------
#Combine data
PSF.ARTR18$Transplant<-rep("ARTR",120)
PSF.HECO18$Transplant<-rep("HECO",120)
PSF.POSE18$Transplant<-rep("POSE",120)
PSF.PSSP18$Transplant<-rep("PSSP",120)
PSF.all18<-rbind(PSF.ARTR18,PSF.HECO18,PSF.POSE18,PSF.PSSP18)
PSF.ARTR19$Transplant<-rep("ARTR",120)
PSF.HECO19$Transplant<-rep("HECO",120)
PSF.POSE19$Transplant<-rep("POSE",120)
PSF.PSSP19$Transplant<-rep("PSSP",120)
PSF.all19<-rbind(PSF.ARTR19,PSF.HECO19,PSF.POSE19,PSF.PSSP19)
PSF.all<-rbind.data.frame(PSF.all18,PSF.all19)
PSF.all$Year<-rep(c(2018,2019),each=480)
#Reorder factor levels
PSF.all$Treatment <- factor(PSF.all$Treatment, levels=c("Control", "Feedback","Exclusion"))
PSF.all$DonorSpp <- factor(PSF.all$DonorSpp, levels=c("ARTR", "BARE","HECO", "POSE","PSSP"))

#Biomass PSF
panel.m<-ggerrorplot(data=PSF.all, x="DonorSpp", y="lnRatio.mass", color="Treatment",
                     facet.by=c("Transplant","Year"),
                     palette=c("#E69F00", "#56B4E9", "#009E73"),
                     xlab="Soil environment microsite",
                     ylab="Aboveground biomass plant-soil feedback log ratio",
                     legend="right")+
  geom_hline(yintercept=0,linetype=2)
panel.m

#Germination PSF
panel.g<-ggerrorplot(data=PSF.all, x="DonorSpp", y="lnRatio.germ", color="Treatment",
                     facet.by=c("Transplant","Year"),
                     palette=c("#E69F00", "#56B4E9", "#009E73"),
                     xlab="Soil environment microsite",
                     ylab="Germination plant-soil feedback log ratio",
                     legend="right")+ 
  geom_hline(yintercept=0,linetype=2)
panel.g


#Is (only ARTR-PSSP "feedback" sig neg)-----------------------------------
#Data pre-processing
PSF.all18<-subset(PSF.all18,DonorSpp!="BARE")#we only care about plants
#Make spp pair label
PSF.all18$Pair[(PSF.all18$DonorSpp=="ARTR"&PSF.all18$Transplant=="POSE")|(PSF.all18$Transplant=="ARTR"&PSF.all18$DonorSpp=="POSE")]<-"ARTR-POSE"
PSF.all18$Pair[(PSF.all18$DonorSpp=="ARTR"&PSF.all18$Transplant=="HECO")|(PSF.all18$Transplant=="ARTR"&PSF.all18$DonorSpp=="HECO")]<-"ARTR-HECO"
PSF.all18$Pair[(PSF.all18$DonorSpp=="ARTR"&PSF.all18$Transplant=="PSSP")|(PSF.all18$Transplant=="ARTR"&PSF.all18$DonorSpp=="PSSP")]<-"ARTR-PSSP"
PSF.all18$Pair[(PSF.all18$DonorSpp=="HECO"&PSF.all18$Transplant=="PSSP")|(PSF.all18$Transplant=="HECO"&PSF.all18$DonorSpp=="PSSP")]<-"HECO-PSSP"
PSF.all18$Pair[(PSF.all18$DonorSpp=="HECO"&PSF.all18$Transplant=="POSE")|(PSF.all18$Transplant=="HECO"&PSF.all18$DonorSpp=="POSE")]<-"HECO-POSE"
PSF.all18$Pair[(PSF.all18$DonorSpp=="PSSP"&PSF.all18$Transplant=="POSE")|(PSF.all18$Transplant=="PSSP"&PSF.all18$DonorSpp=="POSE")]<-"PSSP-POSE"
#Combine all 2019 PSF data
PSF.all19<-subset(PSF.all19,DonorSpp!="BARE")#we only care about plants
#Make spp pair label
PSF.all19$Pair[(PSF.all19$DonorSpp=="ARTR"&PSF.all19$Transplant=="POSE")|(PSF.all19$Transplant=="ARTR"&PSF.all19$DonorSpp=="POSE")]<-"ARTR-POSE"
PSF.all19$Pair[(PSF.all19$DonorSpp=="ARTR"&PSF.all19$Transplant=="HECO")|(PSF.all19$Transplant=="ARTR"&PSF.all19$DonorSpp=="HECO")]<-"ARTR-HECO"
PSF.all19$Pair[(PSF.all19$DonorSpp=="ARTR"&PSF.all19$Transplant=="PSSP")|(PSF.all19$Transplant=="ARTR"&PSF.all19$DonorSpp=="PSSP")]<-"ARTR-PSSP"
PSF.all19$Pair[(PSF.all19$DonorSpp=="HECO"&PSF.all19$Transplant=="PSSP")|(PSF.all19$Transplant=="HECO"&PSF.all19$DonorSpp=="PSSP")]<-"HECO-PSSP"
PSF.all19$Pair[(PSF.all19$DonorSpp=="HECO"&PSF.all19$Transplant=="POSE")|(PSF.all19$Transplant=="HECO"&PSF.all19$DonorSpp=="POSE")]<-"HECO-POSE"
PSF.all19$Pair[(PSF.all19$DonorSpp=="PSSP"&PSF.all19$Transplant=="POSE")|(PSF.all19$Transplant=="PSSP"&PSF.all19$DonorSpp=="POSE")]<-"PSSP-POSE"

#Calculate Is for all species pairs
source("./code/func_IScalc.R")
IS2018<-IScalc(PSF.all18)
IS2019<-IScalc(PSF.all19)

#Check sample sizes
table(na.omit(IS2018)$Pair,na.omit(IS2018)$Treatment)#range 3 to 10
table(na.omit(IS2019)$Pair,na.omit(IS2019)$Treatment)#range 4 to 10


#Use means parameterization to figure out difference from zero
IS2018$levels<-paste(IS2018$Treatment,IS2018$Pair,sep=".")
IS2019$levels<-paste(IS2019$Treatment,IS2019$Pair,sep=".")

#Biomass
m.IS2018.m<-lm(Is.m~levels-1,data=IS2018)
summary(m.IS2018.m)#Only Feedback ARTR-PSSP sig neg
m.IS2019.m<-lm(Is.m~levels-1,data=IS2019)
summary(m.IS2018.m)#Only Control ARTR-PSSP sig neg
#Germination
m.IS2018.g<-lm(Is.g~levels-1,data=IS2018)
summary(m.IS2018.g)#Nothing
m.IS2019.g<-lm(Is.g~levels-1,data=IS2019)
summary(m.IS2018.g)#Nothing

#Anova both years together 
ISall<-rbind.data.frame(IS2018,IS2019)
ISall$Year<-rep(c(2018,2019),each=180)
#Biomass
m.IS.m<-lm(Is.m~Treatment*Pair*factor(Year),data=ISall)
qqPlot(m.IS.m)#fine
Anova(m.IS.m)#no sig effects; Pair is 0.06
#Germ
m.IS.g<-lm(Is.g~Treatment*Pair*factor(Year),data=ISall)
qqPlot(m.IS.g)#fine
Anova(m.IS.g)#no sig effects

#Visualize IS----------------------------- 
ISall$Treatment<-factor(ISall$Treatment, levels=c("Control", "Feedback","Exclusion"))

#Biomass
P.IS.m<-ggerrorplot(ISall,x = "Pair", y = "Is.m",
                    orientation = "horizontal",
                    facet.by = "Year",
                    color = "Treatment",
                    palette=c("#E69F00", "#56B4E9", "#009E73"),
                    xlab="Species Pair", ylab="Pairwise Interaction Strength")+
  geom_hline(yintercept=0,linetype=2)
P.IS.m

#Germination
P.IS.g<-ggerrorplot(ISall,x = "Pair", y = "Is.g",
                    orientation = "horizontal",
                    facet.by = "Year",
                    color = "Treatment",
                    palette=c("#E69F00", "#56B4E9", "#009E73"),
                    xlab="Species Pair", ylab="Pairwise Interaction Strength")+
  geom_hline(yintercept=0,linetype=2)
P.IS.g


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
