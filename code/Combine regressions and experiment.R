# AC 190305

rm(list=ls())
#setwd("C:/Repos/PSFcoexist/code")
setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist/code")

# load libraries
library(maptools)
library(sp)
library(rgeos)
library(dplyr)

# Input vars using means from year 2 field data
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

# Loop combinations of test distances----------
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

# Concatenate list of tables ----------------
neighborFocal.dat <- bind_rows(tableList, .id = "column_label")

# Visualize effects of neighbor size and distance----------------------
library(ggplot2)
library(ggpubr)

# Plot function
f.plotNeighborEffect <- function(focalsp,neighborsp,data){
  subdata <- subset(data, focal==focalsp&neighbor==neighborsp)
  p <- ggplot(subdata, aes(x=neighborDist, y=newSize, color=factor(neighborRadius))) +
    geom_point() +
    geom_line() +
    theme_pubr(border=T,legend="right")+
    labs(title=paste0("focal=",focalsp," neighbor=",neighborsp), 
         color="Neighbor\nRadius", x="Neighbor Distance", y="Focal size (t+1)")
  return(p)
}

# Test plot function
f.plotNeighborEffect("PSSP","HECO",neighborFocal.dat)

# Calculate ln ratios from regression data-------------------
# First assume that no focal plants are under neighbor plant canopies
neighborFocal.dat1<-subset(neighborFocal.dat,neighborDist>neighborRadius|neighborDist==neighborRadius)

# Make empty dataframe
obsPSF.table <- cbind.data.frame(focal=rep(0,16),
                                 neighbor=rep(0,16),
                                 lnR.25=rep(0,16),
                                 lnR.75=rep(0,16),
                                 lnR.sd=rep(0,16),
                                 lnR.mean=rep(0,16))

# Calculate conspecific/heterospecific effects
for (i in 1:4){
  focal<-as.character(neighborFocal.dat$focal)[i]
  home <- neighborFocal.dat[which(neighborFocal.dat$focal==focal&neighborFocal.dat$neighbor==focal),]
  away <- neighborFocal.dat[which(neighborFocal.dat$focal==focal&neighborFocal.dat$neighbor!=focal),]
  
  for (j in 1:4){
    neighbor<-as.character(unique(away$neighbor)[j])
    neighbordat<-away[which(away$neighbor==neighbor),]
    combo<-expand.grid(home$newSize,neighbordat$newSize)
    lnRatio<-log(combo[,1]/combo[,2])
    obsPSF.table[(4*(i-1)+j),]<-c(focal,neighbor,quantile(lnRatio,0.25),quantile(lnRatio,0.75),sd(lnRatio),mean(lnRatio)) #focus on the range where 50% of values fall
  }
}
#All ranges negative values

# Process field experiment control data only -------------
all<-read.csv("../data/Transplant data clean 190102.csv")
control<-subset(all,Treatment=="Control")

# Make empty dataframe
expPSF.table <- cbind.data.frame(Transplant=character(),
                                 DonorSpp=character(),
                                 Rep=integer(),
                                 lnRatio=double(),
                                 stringsAsFactors=F)
for (k in 1:4){
  focal <- as.character(neighborFocal.dat$focal)[k]
  
  for (i in 1:10){
    subdata<-control[which(control$Rep==i&control$Transplant==focal),]
    home<-subdata[which(subdata$DonorSpp==focal),]
    away<-subdata[which(subdata$DonorSpp!=focal),]
    
    for (j in 1:4){
      Donor<-unique(away$DonorSpp)[j]
      lnRatio<-log(home$Abv/away$Abv[away$DonorSpp==Donor])
      expPSF.table<-rbind(expPSF.table,list(Transplant=as.character(focal),
                                            DonorSpp=as.character(Donor),
                                            Rep=i,
                                            lnRatio=lnRatio),
                          stringsAsFactors=FALSE)
    }
    
  }
}


# Combine regression and field experiment PSF estimates in visualization -----------------
library(ggpubr)

#Combine experiment and observation data
combPSF.table<-obsPSF.table
for (i in 1:length(obsPSF.table[,1])){
  subdata<-subset(expPSF.table,Transplant==obsPSF.table$focal[i]&DonorSpp==obsPSF.table$neighbor[i])
  combPSF.table$exp.lnR.mean[i]<-mean(subdata$lnRatio,na.rm = T)
  combPSF.table$exp.lnR.sd[i]<-sd(subdata$lnRatio,na.rm = T)
}

#Fix character/numeric issue s
combPSF.table$lnR.25<-as.numeric(combPSF.table$lnR.25)
combPSF.table$lnR.75<-as.numeric(combPSF.table$lnR.75)
combPSF.table$lnR.sd<-as.numeric(combPSF.table$lnR.sd)
combPSF.table$lnR.mean<-as.numeric(combPSF.table$lnR.mean)

#Plot
ggplot(data=combPSF.table)+
  geom_point(aes(x=neighbor,y=exp.lnR.mean))+
  geom_errorbar(mapping=aes(x=neighbor, ymin=exp.lnR.mean-exp.lnR.sd, ymax=exp.lnR.mean+exp.lnR.sd), width=0.2)+
  geom_point(aes(x=neighbor,y=lnR.mean),col="blue")+
  geom_errorbar(mapping=aes(x=neighbor,ymin=lnR.mean-lnR.sd,ymax=lnR.mean+lnR.sd),width=0.2,col="blue")+
  facet_wrap( ~ focal, ncol=2)+
  theme_bw()+
  ylab("ln(home/away)")+xlab("Donor/neighbor species")

#Not sure why "dodge" isn't working but otherwise good.

