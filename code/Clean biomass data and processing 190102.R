############################################################
#This code takes the raw data with corrected labels and 
#does cleaning and processing for downstream analyses  
##############################################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)

#Raw data import
raw<-read.csv("./data/Transplant data_labels corrected 190102.csv")

#Derive values of interest---------------
#Area
raw$Area<-raw$Length180709/2*raw$Width180709/2*pi #calculate as ellipse

#Belowground biomass
#Examine dry-wet ratios
dryWet<-raw$Blw/raw$Blw.sub
boxplot(dryWet~raw$Transplant)#2 outliers that are >1
#Take out outliers from raw data
out<-which(dryWet>1)
raw$Blw[out]<-NA
raw$Blw.sub[out]<-NA
#recalculate
dryWet<-raw$Blw/raw$Blw.sub
boxplot(dryWet~raw$Transplant)#species-specific ratios
summary(aov(dryWet~raw$Transplant))#super significant
ARTR.ratio<-mean(na.exclude(dryWet[which(raw$Transplant=="ARTR")]))
HECO.ratio<-mean(na.exclude(dryWet[which(raw$Transplant=="HECO")]))
POSE.ratio<-mean(na.exclude(dryWet[which(raw$Transplant=="POSE")]))
PSSP.ratio<-mean(na.exclude(dryWet[which(raw$Transplant=="PSSP")]))
#Process data
raw$Blw.all<-NA
for (i in 1:length(raw[,1])){
  if (!is.na(raw$Blw[i])){
    raw$Blw.all[i]<-raw$Blw.wet.num[i]*raw$Blw[i]/raw$Blw.sub[i]
    next
  }
  if (!is.na(raw$Blw.wet.num[i])){#use mean spp ratios for ones without subsamples
    if (raw$Transplant[i]=="ARTR"){
      raw$Blw.all[i]<-raw$Blw.wet.num[i]*ARTR.ratio
    }
    if (raw$Transplant[i]=="HECO"){
      raw$Blw.all[i]<-raw$Blw.wet.num[i]*HECO.ratio
    }
    if (raw$Transplant[i]=="POSE"){
      raw$Blw.all[i]<-raw$Blw.wet.num[i]*POSE.ratio
    }
    if (raw$Transplant[i]=="PSSP"){
      raw$Blw.all[i]<-raw$Blw.wet.num[i]*PSSP.ratio
    }
  }

}
#Check for outliers
hist(raw$Blw.all)#looks fine
boxplot(raw$Blw.all~raw$Transplant)#looks fine

#Repro
raw$Repro180514[is.na(raw$Repro180514)]<-0
raw$Repro180709[is.na(raw$Repro180709)]<-0
raw$Repro<-raw$Repro180514+raw$Repro180709
raw$Repro[raw$Repro>1]<-1

#Germination (did not use April census bc low confidence in seedling ID, likely all weeds)
raw$Germ180514[is.na(raw$Germ180514)]<-0
raw$Germ180531[is.na(raw$Germ180531)]<-0
raw$Germ180709[is.na(raw$Germ180709)]<-0
raw$Germ180716[is.na(raw$Germ180716)]<-0
raw$Germ<-raw$Germ180514+raw$Germ180531+raw$Germ180709+raw$Germ180716
hist(raw$Germ)
subset(raw,Germ>10)#only one entry with 11, accompanying notes confirm all definitely ARTR. Leave in as is.

#Survival
raw$Surv<-raw$Surv180709
raw$Surv[is.na(raw$Surv)]<-1

#Height
raw$Height<-raw$Height180709

#OUTPUT clean data with only columns of interest-------------------
clean<-raw[,c(1:6,10,32,36,40:45)]
write.csv(clean,"./data/Transplant data clean 190201.csv")
