########################################
#Code to visualize obs data simulations#
########################################

rm(list=ls())
setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist/code")

# load libraries
library(maptools)
library(sp)
library(rgeos)
library(dplyr)
library(ggplot2)
library(car)
library(ggpubr)

#Generate simulated data using regression code--------
neighborDist <-  data.frame(neighbor = c("ARTR","HECO","POSE","PSSP","BARE"),
                            dist = c(35,13,12,15,NA))  # distance from focal individual's centroid to neighbor's centroid
neighborRadius <- data.frame(neighbor = c("ARTR","HECO","POSE","PSSP","BARE"),
                             radius = c(35,8,6,8,NA))
dist <- c(-2,-1,0,1,2,3) #potential variation in distances (SD for distance from focal to neighbor edge ~2 in field)
dist.x <- expand.grid(dist=dist, rad=dist)

# Fixed vars
focalXY <- c(0,0)
focalSize <- 1   # in cm^2
annuli <- c(seq(2,20,2),seq(25,50,5),seq(60,150,10))
speciesList <- c("ARTR","HECO","POSE","PSSP")
neighborList<-c("ARTR","HECO","POSE","PSSP","BARE")
# function to make annuli around focal plant
makeDiscs <- function(x,y,annuli){
  centroid<-readWKT(paste("POINT(",x," ",y,")",sep=""))
  n <- length(annuli)
  output <- list(n)
  for(ii in 1:n){
    output[[ii]] <-gBuffer(centroid,width=annuli[ii],quadsegs=50)
  }
  return(output)
}
# make rings around focal plant 
rings <- makeDiscs(focalXY[1],focalXY[2],annuli)
# import distance weights
distWts <- read.csv("../data/LongTermData/IdahoModDistanceWeights.csv")

# Loop combinations of test distances
tableList <- list() #create empty list to store tables
for (k in 1:length(dist.x[,1])){
  # Make a new dataframe each k to store output
  newSize.table <- cbind.data.frame(expand.grid(focal=speciesList,neighbor=neighborList),
                                    neighborRadius=rep(0,20),
                                    neighborDist=rep(0,20),
                                    newSize=rep(0,20))
  tableList[[k]] <- newSize.table
  # Loop to regress all
  for (i in 1:20){
    neighborSpp <- newSize.table$neighbor[i]
    focalSpp <- newSize.table$focal[i]
    
    # import focal species' regression parameters
    parms <- read.csv(paste0("../data/LongTermData/",focalSpp,"_growth.csv"))
    
    # calculate crowding for plant neighbors
    W <- rep(0,length(speciesList))
    names(W) <- paste("W.",speciesList,sep="")
    if(neighborSpp!="BARE"){
      # make neighbor
      neighbor <- makeDiscs(focalXY[1],focalXY[2]+neighborDist$dist[neighborSpp]+dist.x$dist[k],
                            neighborRadius$radius[neighborSpp]+dist.x$rad[k])
      
      # calculate overlapping areas of focal annuli and neighbor
      Nrings <- length(annuli)
      nbarea <- numeric(Nrings)
      for(j in 1:Nrings){
        A <-gIntersection(rings[[j]],neighbor[[1]])
        if(!is.null(A)){
          nbarea[j] <- gArea(A) # calculate area only if A is not NULL
        }else{
          nbarea[j] <- 0
        }
      } # next j  
      nbarea <- nbarea-c(0,nbarea[1:(Nrings-1)]) # calculate ring-specific areas
      
      doSpp <- which(speciesList==neighborSpp)
      W[doSpp] <- nbarea%*%distWts[,neighborSpp] 
    }
    
    # get interaction coefficients
    alphas <- parms[1,c("W.ARTR","W.HECO","W.POSE","W.PSSP")]
    
    # do regression (assume average year and average location--no random effects)
    newLogSize <- parms$X.Intercept.[1] + parms$logarea.t0[1]*log(focalSize) + sum(alphas*W) 
    newSize <- exp(newLogSize)
    
    # add to table
    tableList[[k]]$newSize[i] <- newSize
    tableList[[k]]$neighborRadius[i] <- neighborRadius$radius[neighborSpp]+dist.x$rad[k]
    tableList[[k]]$neighborDist[i] <- neighborDist$dist[neighborSpp]+dist.x$dist[k]
  }
  
}

#Make dataframe
neighborFocal.dat <- bind_rows(tableList, .id = "column_label")
#Assume that no focal plants are under neighbor plant canopies
neighborFocal.dat1<-subset(neighborFocal.dat,neighborDist>neighborRadius|neighborDist==neighborRadius)




#Data reaaranging for plotting----------
ARTR<-subset(neighborFocal.dat,focal=="ARTR")
POSE<-subset(neighborFocal.dat,focal=="POSE")
PSSP<-subset(neighborFocal.dat,focal=="PSSP")
HECO<-subset(neighborFocal.dat,focal=="HECO")

ARTR$Type[ARTR$neighbor=="ARTR"]<-"Intra"
ARTR$Type[ARTR$neighbor!="ARTR"]<-"Inter"
ARTR$Type <- factor(ARTR$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

POSE$Type[POSE$neighbor=="POSE"]<-"Intra"
POSE$Type[POSE$neighbor!="POSE"]<-"Inter"
POSE$Type <- factor(POSE$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

PSSP$Type[PSSP$neighbor=="PSSP"]<-"Intra"
PSSP$Type[PSSP$neighbor!="PSSP"]<-"Inter"
PSSP$Type <- factor(PSSP$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

HECO$Type[HECO$neighbor=="HECO"]<-"Intra"
HECO$Type[HECO$neighbor!="HECO"]<-"Inter"
HECO$Type <- factor(HECO$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

all<-rbind.data.frame(ARTR,POSE,PSSP,HECO)

plotdf<-data.frame(aggregate(newSize~focal*neighbor*Type,data=all,mean))
plotdf$sd<-data.frame(aggregate(newSize~focal*neighbor*Type,data=all,sd))[,4]

#Visualize---------------------------
ggplot(data=plotdf,aes(color=neighbor))+
  geom_point(aes(x=Type,y=newSize),position=position_dodge(width = 0.8),size=2)+
  geom_errorbar(mapping=aes(x=Type, ymin=newSize-sd, ymax=newSize+sd), 
                width=0.3,position=position_dodge(width = 0.8))+
  facet_wrap( ~ focal, ncol=2)+
  theme_bw()+
  ylab("size(t+1)")+xlab("")
  
