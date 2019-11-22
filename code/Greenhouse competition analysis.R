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

#Try fitting Ricker-type competition model ------------------------------
ARTR.lambda<-mean(ARTR$F.mass[ARTR$Type=="Solo"])#23.54286
POSE.lambda<-mean(POSE$F.mass[POSE$Type=="Solo"])#33.4
PSSP.lambda<-mean(PSSP$F.mass[PSSP$Type=="Solo"])#26.45

#Make empty matrix 9 rows
Ricker.fit<- data.frame(focal=rep(c("ARTR","PSSP","POSE"),3),
                        neighbor=rep(c("ARTR","PSSP","POSE"),each=3),
                        beta=rep(0,9),
                        intcp=rep(0,9),
                        p=rep(0,9),
                        adj.R=rep(0,9))
#Run a big function-fitting loop
for (i in 1:9){
  subdata<-subset(data.all,NeighborNew==Ricker.fit$neighbor[i]&FocalNew==Ricker.fit$focal[i])
  intercept<-log(get(paste(Ricker.fit$focal[i],".lambda",sep="")))
  fit<-lm(I(log(F.mass)-intercept)~0+N.total,data=subdata)
  Ricker.fit$beta[i]<-coef(fit)
  Ricker.fit$intcp[i]<-intercept
  Ricker.fit$p[i]<-summary(fit)[["coefficients"]][4]
  Ricker.fit$adj.R[i]<-summary(fit)[["adj.r.squared"]]
}
#all significant negative slopes except for ARTR-ARTR and PSSP-ARTR
#Intraspecific comp coeffs are not all bigger than interspecific comp coeffs
ggplot(Ricker.fit,aes(focal,neighbor))+
  geom_count()

#Plot all the fits
Ricker.fit.plot<-Ricker.fit
colnames(Ricker.fit.plot)[1:2]<-c("FocalNew","NeighborNew")#make the facet variable the same
ggplot(data.all, aes(N.total, log(F.mass)))+
  geom_point()+
  xlab("Number of neighbors")+ ylab("Focal mass (mg)")+
  geom_abline(data=Ricker.fit.plot, aes(intercept = intcp, slope = beta))+
  facet_grid(FocalNew~NeighborNew)
#Not bad!  
