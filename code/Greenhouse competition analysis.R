#This code analyzes the 2019 greenhouse competition experiment
#AC last updated 190818

setwd("C:/Users/yyach/Box Sync/PSFcoexist/R/PSFcoexist")

#Import cleaned data, subset, and visualize-----------------------
data<-read.csv("./data/Greenhouse comp data clean 190818.csv", row.names=1)

#Forgot to calculate total number of neighbors at harvest in the data cleaning code
data$N.total<-data$N.dead+data$N.live
data$N.total[is.na(data$N.total)]<-0

#Also, replace NA's in total neighbor mass with 0 for "no-neighbor" situations
data$N.all.mass[is.na(data$N.all.mass)]<-0

#First viz
boxplot(data$F.mass~data$Type*data$FocalNew) #can tentatively tell that solo plants are bigger

#Subset
ARTR<-subset(data,FocalNew=="ARTR") #40 obs
PSSP<-subset(data,FocalNew=="PSSP") #68 obs
POSE<-subset(data,FocalNew=="POSE") #72 obs

#Species-specific viz
plot(ARTR$N.total,ARTR$F.mass,col=as.numeric(ARTR$NeighborNew)) #almost "envelope" situation, mostly neighbor spp-agnostic
plot(PSSP$N.total,PSSP$F.mass,col=as.numeric(PSSP$NeighborNew)) #also "envelope" situation, mostly neighbor spp-agnostic
plot(POSE$N.total,POSE$F.mass,col=as.numeric(POSE$NeighborNew)) #also "envelope" situation, mostly neighbor spp-agnostic
#Try using neighbor mass instead
plot(ARTR$N.all.mass,ARTR$F.mass,col=as.numeric(ARTR$NeighborNew)) #still envelope, neighbor species split along mass axis
plot(PSSP$N.all.mass,PSSP$F.mass,col=as.numeric(PSSP$NeighborNew)) #also "envelope" situation, mostly neighbor spp-agnostic
plot(POSE$N.all.mass,POSE$F.mass,col=as.numeric(POSE$NeighborNew)) #this actually doesn't even look like NDD, mostly neighbor spp-agnostic

