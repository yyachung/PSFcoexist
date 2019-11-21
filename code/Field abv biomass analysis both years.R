# Combine 2018 and 2019 field data for biomass analysis #
# Update AC 191121

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)
library(emmeans)
library(ggpubr)

#Data import------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Combine data and add year column----------------
sub2018<-all2018[,c("newID","DonorSpp","Transplant","Treatment","Rep","Abv")]
sub2018$Year<-rep(2018,600)
names(sub2018)[1]<-"ID"
sub2019<-all2019[,c("fieldID","DonorSpp","Transplant","Treatment","Rep","Abv")]
sub2019$Year<-rep(2019,600)
names(sub2019)[1]<-"ID"
alldat<-rbind(sub2018,sub2019)
alldat<-subset(alldat,Treatment%in%c("Control","Exclusion","Feedback"))#get rid of the ones with no treatment labels from 2018

#Prelim visualization to check distribution------------------
hist(alldat$Abv)#probably needs a log transform
hist(alldat$Abv[alldat$Year==2018])#same distribution within-year
hist(alldat$Abv[alldat$Year==2019])

#Separate by species-------------------------------------------
ARTR<-subset(alldat,Transplant=="ARTR")
PSSP<-subset(alldat,Transplant=="PSSP")
POSE<-subset(alldat,Transplant=="POSE")
HECO<-subset(alldat,Transplant=="HECO")

#Global lm model for ARTR (Treatment and some interaction effects)------------------------------------------
gARTR<-lm(log(Abv)~Treatment*DonorSpp*factor(Year),data=ARTR)
qqPlot(gARTR$resid)#not bad
hist(gARTR$resid)#also decent
Anova(gARTR,type=3)#Treatment p=0.002636
summary(gARTR)
#TreatmentExclusion p=0.00161 
#TreatmentFeedback:DonorSppPSSP:factor(Year)2019 p=0.0395 

#Emmeans comparisons: contrasting everything else with ARTR soil within treatment and years
emmeans(gARTR, trt.vs.ctrl ~ DonorSpp | Year*Treatment)
#2018 Feedback ARTR - HECO p=0.022 ARTR - PSSP p=0.0311
#2019 Control PSSP - ARTR p=0.0889

#Visualization
boxplot(log(ARTR$Abv)~ARTR$DonorSpp*ARTR$Treatment*ARTR$Year)#still expected pattern
#Estimated marginal means
emmip(gARTR, Treatment ~ DonorSpp | Year)
#Pretty plot
ggerrorplot(ARTR, x = "DonorSpp", y = "Abv", 
            desc_stat = "mean_se", 
            color = "Treatment", 
            facet.by = "Year",
            position = position_dodge(0.3)
            )

#Global lm model for HECO (Treatment effects)------------------------------------------
gHECO<-lm(log(Abv)~Treatment*DonorSpp*factor(Year),data=HECO)
qqPlot(gHECO$resid)#not great
hist(gHECO$resid)#not as bad looking
gHECO.1<-lm(Abv~Treatment*DonorSpp*Year,data=HECO)
qqPlot(gHECO.1$resid)#worse

#Keep with the logged version
Anova(gHECO,type=3)#Treatment p=0.002179
summary(gHECO)#TreatmentExclusion p=0.000582

#Emmeans comparisons: contrasting everything else with HECO soil within treatment and years
emmeans(gHECO, trt.vs.ctrl ~ DonorSpp | Year*Treatment, ref=3)
#2019 Feedback  ARTR - HECO p=0.0494 
#2018 Exclusion BARE - HECO p=0.0803

#Visualization
boxplot(log(HECO$Abv)~HECO$DonorSpp*HECO$Treatment*HECO$Year)#still expected pattern; CIs really wide
emmip(gHECO, Treatment ~ DonorSpp | Year)
#Pretty plot
ggerrorplot(HECO, x = "DonorSpp", y = "Abv", 
            desc_stat = "mean_se", 
            color = "Treatment", 
            facet.by = "Year",
            position = position_dodge(0.3)
)


#Global lm model for POSE (Treatment effects and some interactions)------------------------------------------
gPOSE<-lm(log(Abv)~Treatment*DonorSpp*factor(Year),data=POSE)
qqPlot(gPOSE$resid)#two outliers (179 and 864)
hist(gPOSE$resid)#not great
POSE[c("179","864"),]#plant in 864 is 2 orders of magnitude smaller than others. Try excluding
subPOSE<-POSE[-which(row.names(POSE)=="864"),]
gPOSE.1<-lm(log(Abv)~Treatment*DonorSpp*factor(Year),data=subPOSE)
qqPlot(gPOSE.1$resid)#much better
hist(gPOSE.1$resid)#much better

#Analyze excluded version
Anova(gPOSE.1,type=3)#Treatment p=2.063e-06
summary(gPOSE.1)
#TreatmentExclusion p=3.19e-07 TreatmentFeedback p=0.0110
#DonorSppHECO:factor(Year)2019   p=0.0144

#Emmeans comparisons: among donor spp within treatment and years
emmeans(gPOSE.1, trt.vs.ctrl ~ DonorSpp | Year*Treatment,ref=4)
#2019 Feedback POSE - PSSP p=0.0211 

#Visualization
boxplot(log(subPOSE$Abv)~subPOSE$DonorSpp*subPOSE$Treatment*subPOSE$Year)#still expected pattern; CIs really wide
emmip(gPOSE.1, Treatment ~ DonorSpp | Year)
#Pretty plot
ggerrorplot(subPOSE, x = "DonorSpp", y = "Abv", 
            desc_stat = "mean_se", 
            color = "Treatment", 
            facet.by = "Year",
            position = position_dodge(0.3)
)

#Global lm model for PSSP (Treatment effects and some interactions)------------------------------------------
gPSSP<-lm(log(Abv)~Treatment*DonorSpp*factor(Year),data=PSSP)
qqPlot(gPSSP$resid)#not bad
hist(gPSSP$resid)#also decent
Anova(gPSSP,type=3)#Treatment effect (p=9.438e-10), Treatment:year (p=0.0319)
summary(gPSSP)
#TreatmentExclusion  p= 1.87e-10    TreatmentFeedback p=03.70e-05  
#TreatmentFeedback:Year  p=0.00900 

#Emmeans comparisons: among donor spp within treatment and years
emmeans(gPSSP, trt.vs.ctrl ~ DonorSpp | Year*Treatment, ref=5)
#2019 Control POSE - PSSP p=0.0274;  BARE - PSSP p=0.0587
#2019 Feedback HECO - PSSP p=0.0376

#Visualization
boxplot(log(PSSP$Abv)~PSSP$DonorSpp*PSSP$Treatment*PSSP$Year)#still expected pattern; CIs really wide
emmip(gPSSP, Treatment ~ DonorSpp | Year)#Feedback effect stronger in 2018 than 2019
#Pretty plot
ggerrorplot(PSSP, x = "DonorSpp", y = "Abv", 
            desc_stat = "mean_se", 
            color = "Treatment", 
            facet.by = "Year",
            position = position_dodge(0.3)
)

#General treatment effect visualization---------------------------------------
ggerrorplot(alldat, x = "DonorSpp", y = "Abv", 
            desc_stat = "mean_se", 
            color = "Treatment",
            size=0.7,
            position = position_dodge(0.3)
)
