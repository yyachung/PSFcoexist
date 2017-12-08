setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R")

#Data
all<-read.csv("Transplant data 171107.csv")
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
