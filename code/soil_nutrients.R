################################################################
#Code to analyze soil nutrient data
#AC 190207
################################################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
setwd("C:/Users/yyach/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)

#Data import and viz-----------------------------
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

#Quick stats Live vs. Sterile soil---------------------------
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

