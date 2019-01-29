# PBA 1/29/2019

rm(list=ls())
setwd("C:/Repos/PSFcoexist/code")

# load libraries
library(maptools)
library(sp)
library(rgeos)

# Set up
focalSpp <- "PSSP"
focalXY <- c(0,0)
focalSize <- 1   # in cm^2
neighborSpp <- "PSSP"
neighborDist <- 5  # distance from focal individual's centroid to neighbor's centroid
neighborRadius <- 5
annuli <- c(seq(2,20,2),seq(25,50,5),seq(60,150,10))
speciesList <- c("ARTR","HECO","POSE","PSSP")


###
###  1. Calculate spatial neighborhood data for a typical focal seedling
###    experiencing crowding from one neighbor

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

# make neighbor
neighbor <- makeDiscs(focalXY[1],focalXY[2] + neighborDist,neighborRadius)

# see if it worked
plot(neighbor[[1]])
lines(rings[[1]]@polygons[[1]]@Polygons[[1]]@coords,col="red")
lines(rings[[2]]@polygons[[1]]@Polygons[[1]]@coords,col="red")

# calculate overlapping areas of focal annuli and neighbor
Nrings <- length(annuli)
nbarea <- numeric(Nrings)
for(i in 1:Nrings){
  A <-gIntersection(rings[[i]],neighbor[[1]])
  if(!is.null(A)){
    nbarea[i] <- gArea(A) # calculate area only if A is not NULL
  }else{
    nbarea[i] <- 0
  }
} # next i  
nbarea <- nbarea-c(0,nbarea[1:(Nrings-1)]) # calculate ring-specific areas  

# # format output
# output <- as.data.frame(matrix(0,0,Nrings))
# names(output) <- paste(rep(neighborSpp,Nrings),annuli,sep=".")
# output[1,] <- nbarea
 

###
###  2. Do regressions
###

# import distance weights
distWts <- read.csv("../data/LongTermData/IdahoModDistanceWeights.csv")

# import focal species' regression parameters
parms <- read.csv(paste0("../data/LongTermData/",focalSpp,"_growth.csv"))
  
# calculate crowding
W <- rep(0,length(speciesList))
names(W) <- paste("W.",speciesList,sep="")
doSpp <- which(speciesList==neighborSpp)
W[doSpp] <- nbarea%*%distWts[,neighborSpp] 

# get interaction coefficients
alphas <- parms[1,c("W.ARTR","W.HECO","W.POSE","W.PSSP")]

# do regression (assume average year and average location--no random effects)
newLogSize <- parms$X.Intercept.[1] + parms$logarea.t0[1]*log(focalSize) + sum(alphas*W) 
newSize <- exp(newLogSize)


