#This code analyzes the 2019 greenhouse competition experiment
#AC last updated 191121

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
setwd("C:/Users/yyach/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)

#Import cleaned data, subset, and visualize-----------------------
data<-read.csv("./data/Greenhouse comp data clean 190818.csv", row.names=1)

#Forgot to calculate total number of neighbors at harvest in the data cleaning code
data$N.total<-data$N.dead+data$N.live
data$N.total[is.na(data$N.total)]<-0

#Also, replace NA's in total neighbor mass with 0 for "no-neighbor" situations
data$N.all.mass[is.na(data$N.all.mass)]<-0

#First viz
boxplot(data$F.mass~data$Type*data$FocalNew) #can tentatively tell that solo plants are bigger

#Subset and only keep those with actual focal data
ARTR<-subset(data,FocalNew=="ARTR"&!is.na(F.mass)) #25 obs
PSSP<-subset(data,FocalNew=="PSSP"&!is.na(F.mass)) #53 obs
POSE<-subset(data,FocalNew=="POSE"&!is.na(F.mass)) #67 obs
#rejoin data
data.all<-rbind.data.frame(ARTR,PSSP,POSE)

#Species-specific viz
plot(ARTR$N.total,ARTR$F.mass,col=as.numeric(ARTR$NeighborNew)) #almost "envelope" situation, mostly neighbor spp-agnostic
plot(PSSP$N.total,PSSP$F.mass,col=as.numeric(PSSP$NeighborNew)) #also "envelope" situation, mostly neighbor spp-agnostic
plot(POSE$N.total,POSE$F.mass,col=as.numeric(POSE$NeighborNew)) #also "envelope" situation, mostly neighbor spp-agnostic
#Try using neighbor mass instead
plot(ARTR$N.all.mass,ARTR$F.mass,col=as.numeric(ARTR$NeighborNew)) #still envelope, neighbor species split along mass axis
plot(PSSP$N.all.mass,PSSP$F.mass,col=as.numeric(PSSP$NeighborNew)) #also "envelope" situation, mostly neighbor spp-agnostic
plot(POSE$N.all.mass,POSE$F.mass,col=as.numeric(POSE$NeighborNew)) #this actually doesn't even look like NDD, mostly neighbor spp-agnostic

#Pretty plot
ggplot(data.all, aes(N.total, F.mass))+
  geom_point()+
  xlab("Number of neighbors")+ ylab("Focal mass (mg)")+
  facet_grid(FocalNew~NeighborNew)
ggplot(data.all, aes(N.all.mass, F.mass))+
  geom_point()+
  xlab("Total neighbor mass (mg)")+ ylab("Focal mass (mg)")+
  facet_grid(FocalNew~NeighborNew)

#Log it
plot(ARTR$N.total,log(ARTR$F.mass),col=as.numeric(ARTR$NeighborNew)) #hard to tell what's going on
plot(PSSP$N.total,log(PSSP$F.mass),col=as.numeric(PSSP$NeighborNew)) #pretty flat
plot(POSE$N.total,log(POSE$F.mass),col=as.numeric(POSE$NeighborNew)) #kind of linear?

#Try fitting Ricker-type competition model (need to re-do this)------------------------------
ARTR.lambda<-mean(ARTR$F.mass[ARTR$Type=="Solo"])#23.54286
POSE.lambda<-mean(POSE$F.mass[POSE$Type=="Solo"])#33.4
PSSP.lambda<-mean(PSSP$F.mass[PSSP$Type=="Solo"])#26.45

#ARTR-POSE
ARTRPOSE<-subset(ARTR,NeighborNew%in%c("ARTR","POSE"))
#Fit model with solo data and estimate intercept
m.AR.PO<-lm(log(F.mass)~N.total:NeighborNew,data=ARTRPOSE)
hist(resid(m.AR.PO))#not terrible..
qqPlot(resid(m.AR.PO))#2 outside CI
summary(m.AR.PO)#no significant slopes. Intercept estimate 2.65 (0.38SE) (ln(23.53)=3.16 outside SE)
#Fit model without solo data
m.AR.PO1<-lm(log(F.mass)~N.total:NeighborNew,data=subset(ARTRPOSE,Type!="Solo"))
hist(resid(m.AR.PO1))#not terrible..
qqPlot(resid(m.AR.PO1))#better
summary(m.AR.PO1)
#Intra: -0.28 (p=0.07) Inter:-0.17 (p=0.04) Intercept estimate 3.09 (0.34SE) (ln(23.53)=3.16 inside SE)


#ARTR-PSSP
ARTRPSSP<-subset(ARTR,NeighborNew%in%c("ARTR","PSSP"))
#Fit model with solo data
m.AR.PS<-lm(log(F.mass)~N.total:NeighborNew,data=ARTRPSSP)
hist(resid(m.AR.PS))#not terrible..
qqPlot(resid(m.AR.PS))#3 outside CI
summary(m.AR.PS)#Inter:-0.39 (p=0.04) Intercept estimate 2.70 (0.43SE) (ln(23.53)=3.16 outside SE)
#Fit model without solo data
m.AR.PS1<-lm(log(F.mass)~N.total:NeighborNew,data=subset(ARTRPSSP,Type!="Solo"))
hist(resid(m.AR.PS1))#not terrible..
qqPlot(resid(m.AR.PS1))#better
summary(m.AR.PS1)#Inter:-0.47 (p=0.007) Intercept estimate 3.18 (0.56SE) (ln(23.53)=3.16 inside SE)

#PSSP-POSE
PSSPPOSE<-subset(PSSP,NeighborNew%in%c("PSSP","POSE"))
#Fit model with solo data and estimate intercept
m.PS.PO<-lm(log(F.mass)~N.total:NeighborNew,data=PSSPPOSE)
hist(resid(m.PS.PO))#not great
qqPlot(resid(m.PS.PO))#5 outside CI
summary(m.PS.PO)#no significant slopes. Intercept estimate 2.90 (0.25SE) (ln(26.45)=3.27 outside SE)
#Fit model without solo data
m.PS.PO1<-lm(log(F.mass)~N.total:NeighborNew,data=subset(PSSPPOSE,Type!="Solo"))
hist(resid(m.PS.PO1))#not terrible..
qqPlot(resid(m.PS.PO1))#better still 4 outside CI
summary(m.PS.PO1)
#Intra:-0.13 (p=0.02) Intercept estimate 3.06 (0.23SE) (ln(26.45)=3.27 inside SE)

#PSSP-ARTR
PSSPARTR<-subset(PSSP,NeighborNew%in%c("PSSP","ARTR"))
#Fit model with solo data and estimate intercept
m.PS.AR<-lm(log(F.mass)~N.total:NeighborNew,data=PSSPARTR)
hist(resid(m.PS.AR))#not great
qqPlot(resid(m.PS.AR))#6 outside CI
summary(m.PS.AR)#no significant slopes. Intercept estimate 2.90 (0.23SE) (ln(26.45)=3.27 outside SE)
#Fit model without solo data
m.PS.AR1<-lm(log(F.mass)~N.total:NeighborNew,data=subset(PSSPARTR,Type!="Solo"))
hist(resid(m.PS.AR1))#not terrible..
qqPlot(resid(m.PS.AR1))#better still 5 outside CI
summary(m.PS.AR1)
#Intra:-0.12 (p=0.02) Intercept estimate 3.01 (0.21SE) (ln(26.45)=3.27 outside SE)

#POSE-ARTR
POSEARTR<-subset(POSE,NeighborNew%in%c("POSE","ARTR"))
#Fit model with solo data and estimate intercept
m.PO.AR<-lm(log(F.mass)~N.total:NeighborNew,data=POSEARTR)
hist(resid(m.PO.AR))#decent
qqPlot(resid(m.PO.AR))#decent
summary(m.PO.AR)
#Intra: -0.15 (p=0.003) Inter:-0.28 (p=0.07) Intercept estimate 2.68 (0.22SE) (ln(23.53)=3.51 outside SE)
#Fit model without solo data
m.PO.AR1<-lm(log(F.mass)~N.total:NeighborNew,data=subset(POSEARTR,Type!="Solo"))
hist(resid(m.PO.AR1))#decent
qqPlot(resid(m.PO.AR1))#decent
summary(m.PO.AR1)
#Intra:-0.12 (p=0.03) Intercept estimate 2.40 (0.27SE) (ln(26.45)=3.27 outside SE)

#POSE-PSSP
POSEPSSP<-subset(POSE,NeighborNew%in%c("POSE","PSSP"))
#Fit model with solo data and estimate intercept
m.PO.PS<-lm(log(F.mass)~N.total:NeighborNew,data=POSEPSSP)
hist(resid(m.PO.PS))#not terrible
qqPlot(resid(m.PO.PS))#not great
summary(m.PO.PS)
#Intra: -0.13 (p=0.021) Intercept estimate 2.54 (0.25SE) (ln(23.53)=3.51 outside SE)
#Fit model without solo data
m.PO.PS1<-lm(log(F.mass)~N.total:NeighborNew,data=subset(POSEPSSP,Type!="Solo"))
hist(resid(m.PO.PS1))#decent
qqPlot(resid(m.PO.PS1))#not bad
summary(m.PO.PS1)
#no sig slopes Intercept estimate 2.26 (0.29SE) (ln(26.45)=3.27 outside SE)