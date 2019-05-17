# AC 190222

rm(list=ls())
#setwd("C:/Repos/PSFcoexist/code")
setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist/code")

# load libraries
library(maptools)
library(sp)
library(rgeos)
library(dplyr)

# Input vars
neighborDist <-  data.frame(neighbor = c("ARTR","HECO","POSE","PSSP","BARE"),
                            dist = c(30,7,5,7,NA))  # distance from focal individual's centroid to neighbor's centroid
neighborRadius <- data.frame(neighbor = c("ARTR","HECO","POSE","PSSP","BARE"),
                             radius = c(30,5,3,5,NA))
dist <- c(-2,-1,0,1,2,3) #potential variation in distances
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
f.plotNeighborEffect("PSSP","POSE",neighborFocal.dat)
