################################################################
#Code to analyze soil nutrient data
#AC 210207
################################################################

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
setwd("C:/Users/yyach/Box Sync/PSFcoexist/R/PSFcoexist")

library(car)
library(lme4)

#Data import and viz 190725-----------------------------
soil <- read.csv("./data/Midwest labs soil analysis 190725.csv")

#Separate out data categories

soil$Type<-substr(soil$Sample.ID,(nchar(as.character(soil$Sample.ID))-1),nchar(as.character(soil$Sample.ID)))
soil$Source<-substr(soil$Sample.ID,1,4)
soil$Sample.pair<-substr(soil$Sample.ID,1,nchar(as.character(soil$Sample.ID))-3)

#Visualization: LV vs. ST
boxplot(soil$OM~soil$Type)#same
boxplot(soil$P1~soil$Type)#same
boxplot(soil$P2~soil$Type)#same
boxplot(soil$K~soil$Type)#same
boxplot(soil$Mg~soil$Type)#same
boxplot(soil$Ca~soil$Type)#same
boxplot(soil$pH~soil$Type)#ame
boxplot(soil$CEC~soil$Type)#same
boxplot(soil$Nitrate~soil$Type)#live higher

#Visualization: Source
boxplot(soil$OM~soil$Source)#HECO lower
boxplot(soil$P1~soil$Source)#PSSP wider range
boxplot(soil$P2~soil$Source)#same
boxplot(soil$K~soil$Source)#same
boxplot(soil$Mg~soil$Source)#same
boxplot(soil$Ca~soil$Source)#same
boxplot(soil$pH~soil$Source)#same
boxplot(soil$CEC~soil$Source)#HECO lower
boxplot(soil$Nitrate~soil$Source)#same

#Quick stats 190725; mixed model update 210207---------------------------
library(lmerTest)
library(emmeans)

#Organic Matter
m.OM<-lm(OM~Type*Source,data=soil)
qqPlot(m.OM$residuals)
Anova(m.OM,type = 2)#Source sig (HECO low) 

mm.OM<-lmer(OM~Type*Source+(1|Sample.pair),data=soil)
summary(mm.OM)
anova(mm.OM)#Source sig (HECO low)

#P1
m.P1<-lm(P1~Type*Source,data=soil)
qqPlot(m.P1$residuals)#not good
m.P1.log<-lm(log(P1)~Type*Source,data=soil)
qqPlot(m.P1.log$residuals)#better
Anova(m.P1.log,type=3)#nothing

mm.P1<-lmer(log(P1)~Type*Source+(1|Sample.pair),data=soil)
summary(mm.P1)
anova(mm.P1)

#P2
m.P2<-lm(P2~Type*Source,data=soil)
qqPlot(m.P2$residuals)
Anova(m.P2,type = 3)#No diff

mm.P2<-lmer(P2~Type*Source+(1|Sample.pair),data=soil)
summary(mm.P2)
anova(mm.P2)

#K
m.K<-lm(K~Type*Source,data=soil)
qqPlot(m.K$residuals)
Anova(m.K,type = 3)#No diff

mm.K<-lmer(K~Type*Source+(1|Sample.pair),data=soil)
summary(mm.K)
anova(mm.K)

#Mg
m.Mg<-lm(Mg~Type*Source,data=soil)
qqPlot(m.Mg$residuals)#weird
m.Mg.log<-lm(log(Mg)~Type*Source,data=soil)
qqPlot(m.Mg.log$residuals)#better
Anova(m.Mg.log,type=3)#nothing

mm.Mg<-lmer(log(Mg)~Type*Source+(1|Sample.pair),data=soil)
summary(mm.Mg)
anova(mm.Mg)#Type is significant (live higher)

#Ca
m.Ca<-lm(Ca~Type*Source,data=soil)
qqPlot(m.Ca$residuals)
Anova(m.Ca,type = 3)#No diff

mm.Ca<-lmer(Ca~Type*Source+(1|Sample.pair),data=soil)
summary(mm.Ca)
anova(mm.Ca)#Sig interaction between type and source
pairs(emmeans(mm.Ca,~Type|Source))#lower in live HECO vs. ST HECO

#pH
m.pH<-lm(pH~Type*Source,data=soil)
qqPlot(m.pH$residuals)
Anova(m.pH,type = 3)#No diff

mm.pH<-lmer(pH~Type*Source+(1|Sample.pair),data=soil)
summary(mm.pH)
anova(mm.pH)#Type sig
pairs(emmeans(mm.pH,~Type))#Live lower

#CEC
m.CEC<-lm(CEC~Type*Source,data=soil)
qqPlot(m.CEC$residuals)
Anova(m.CEC,type = 3)#Source p=0.024 (HECO low)

mm.CEC<-lmer(CEC~Type*Source+(1|Sample.pair),data=soil)
summary(mm.CEC)
anova(mm.CEC)#Sig interaction between type and source
pairs(emmeans(mm.CEC,~Type|Source))#lower in live HECO vs. ST HECO

#Nitrate
m.N<-lm(Nitrate~Type*Source,data=soil)
qqPlot(m.N$residuals)
Anova(m.N,type = 2)#no diff

mm.N<-lmer(Nitrate~Type*Source+(1|Sample.pair),data=soil)
summary(mm.N)
anova(mm.N)#Sig type
pairs(emmeans(mm.N,~Type))#Live higher


#Paired LV-ST contrasts 190725-----------------------
#Make ID column
soil$Pair<-substr(soil$Sample.ID,1,6)
table(soil$Pair)#worked as expected

#Paired T tests

#Organic Matter
t.test(soil$OM[which(soil$Type=="LV")],soil$OM[which(soil$Type=="ST")],paired=T)#nothing

#P1
t.test(soil$P1[which(soil$Type=="LV")],soil$P1[which(soil$Type=="ST")],paired=T)#nothing

#P2
t.test(soil$P2[which(soil$Type=="LV")],soil$P2[which(soil$Type=="ST")],paired=T)#nothing

#K
t.test(soil$K[which(soil$Type=="LV")],soil$K[which(soil$Type=="ST")],paired=T)#nothing

#Mg
t.test(soil$Mg[which(soil$Type=="LV")],soil$Mg[which(soil$Type=="ST")],paired=T)#nothing

#Ca
t.test(soil$Ca[which(soil$Type=="LV")],soil$Ca[which(soil$Type=="ST")],paired=T)#nothing

#pH
t.test(soil$pH[which(soil$Type=="LV")],soil$pH[which(soil$Type=="ST")],paired=T)#nothing

#CEC
t.test(soil$CEC[which(soil$Type=="LV")],soil$CEC[which(soil$Type=="ST")],paired=T)#nothing

#Nitrate
t.test(soil$Nitrate[which(soil$Type=="LV")],soil$Nitrate[which(soil$Type=="ST")],paired=T)
#nitrate higher in live soils (t = 4.9791, df = 39, p-value = 1.338e-05)
plot(jitter(soil$Nitrate[which(soil$Type=="LV")]),jitter(soil$Nitrate[which(soil$Type=="ST")]),
     xlab="live soil nitrate",ylab="sterile soil nitrate")
abline(a=0,b=1)#add 1:1 line

#Data import and viz 190205-----------------------------
soil <- read.csv("./data/Midwest labs soil analysis 190205.csv")

#Visualization: LV vs. ST
boxplot(soil$OM~soil$Type)#live slightly higher
boxplot(soil$P1~soil$Type)#sterile slightly higher
boxplot(soil$P2~soil$Type)#sterile slightly higher
boxplot(soil$K~soil$Type)#same
boxplot(soil$Mg~soil$Type)#sterile slightly higher
boxplot(soil$Ca~soil$Type)#same
boxplot(soil$pH~soil$Type)#sterile slightly higher
boxplot(soil$CEC~soil$Type)#same
boxplot(soil$Nitrate~soil$Type)#live higher

#Visualization: Source
boxplot(soil$OM~soil$Source)#BARE and BULK lower than under plants
boxplot(soil$P1~soil$Source)#Higher under grasses
boxplot(soil$P2~soil$Source)#lowest under POSE
boxplot(soil$K~soil$Source)#BARE and BULK lower than under plants
boxplot(soil$Mg~soil$Source)#BARE and BULK higher than under plants
boxplot(soil$Ca~soil$Source)#much higher under bulk than others
boxplot(soil$pH~soil$Source)#same
boxplot(soil$CEC~soil$Source)#much higher under bulk than others
boxplot(soil$Nitrate~soil$Source)#same except for PSSP outlier

#Quick stats Live vs. Sterile soil 190205---------------------------
#Organic Matter
m.OM<-lm(OM~Type,data=soil)
qqPlot(m.OM$residuals)
summary(m.OM)#No diff

#P1
m.P1<-lm(P1~Type,data=soil)
qqPlot(m.P1$residuals)
summary(m.P1)#p=0.0506 ST higher

#P2
m.P2<-lm(P2~Type,data=soil)
qqPlot(m.P2$residuals)
summary(m.P2)#No diff

#K
m.K<-lm(K~Type,data=soil)
qqPlot(m.K$residuals)
summary(m.K)#No diff

#Mg
m.Mg<-lm(Mg~Type,data=soil)
qqPlot(m.Mg$residuals)
summary(m.Mg)#No diff

#Ca
m.Ca<-lm(Ca~Type,data=soil)
qqPlot(m.Ca$residuals)#not great
summary(m.Ca)#No diff
#Try transformation
hist(soil$Ca)
hist(log(soil$Ca))#not sure if better
m.logCa<-lm(log(Ca)~Type,data=soil)
qqPlot(m.logCa$residuals)#slightly better
summary(m.logCa)#Still no diff

#pH
m.pH<-lm(pH~Type,data=soil)
qqPlot(m.pH$residuals)
summary(m.pH)#Significantly higher pH in ST (p=0.02)

#CEC
m.CEC<-lm(CEC~Type,data=soil)
qqPlot(m.CEC$residuals)
summary(m.CEC)#No diff

#Nitrate
m.N<-lm(Nitrate~Type,data=soil)
qqPlot(m.N$residuals)#not good (one outlier)
summary(m.N)#ST lower nitrate p=0.0003
#Try transformation
hist(soil$Nitrate)
hist(log(soil$Nitrate))
m.logN<-lm(log(Nitrate)~Type,data=soil)
qqPlot(m.logN$residuals)#fine
summary(m.logN)#ST lower nitrate p<0.0001
#where's the outlier?
soil[which(soil$Nitrate>10),]#The PSSP-LV sample 18ppm
#compare with corresponding ST sample
soil[which(soil$Source=="PSSP"),]#corresponding ST has only 2ppm??

