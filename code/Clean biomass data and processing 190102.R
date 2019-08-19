############################################################
#This code takes the raw data with corrected labels and 
#does cleaning and processing for downstream analyses 
#AC update 190818
##############################################################

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
setwd("C:/Users/yyach/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)

################Field data 2018/2019 cleaning###############################
#raw data import
raw<-read.csv("./data/Transplant data year 2_labels corrected 190819.csv")

#Area
raw$Area<-raw$Long_190620/2*raw$Short_190620/2*pi #calculate as ellipse

#Survival
raw$Surv<-raw$Surv_190630

#Germination
germ1<-apply(raw[,c("Germ_190424","Germ_190509")],1,FUN="max") #Calculate the max of the 2 first germ surveys since I didn't pull germinants after first survey
raw$Germ<-raw$Germ_190605+germ1 #add the last census
summary(raw$Germ) #looks great

#Reproduction
raw$Repro_190605[is.na(raw$Repro_190605)]<-0
raw$Repro<-raw$Repro_190605+raw$Repro_190620+raw$Repro_190630
hist(raw$Repro)#for the ones counted more than once return to just 1
raw$Repro[raw$Repro>1]<-1
hist(raw$Repro)

#Height
raw$Height<-rep(0,600)
raw$Height[raw$Ruler=="Silver"]<-raw$Height_190620[raw$Ruler=="Silver"]+0.3 #account for bit at bottom of ruler
raw$Height[raw$Ruler=="Pink"]<-raw$Height_190620[raw$Ruler=="Pink"]+0.5 #account for bit at bottom of ruler
raw$Height[raw$Height==0]<-NA #Make sure missing values are not coded as zero

#Fungus
raw$Fungus[is.na(raw$Fungus)]<-0

#Calculate diameter of donor plant
raw$Diam_donor<-raw$Dist_centroid-raw$Dist_edge

#Output
clean<-raw[,c("fieldID","DonorSpp","Transplant","Treatment","Transect","Rep","Position","Offset","Side",
              "Transplant.date","Dist_centroid","Dist_edge","Abv.mass","Area","Surv","Germ","Repro","Height","Diam_donor")]
write.csv(clean,"./data/Transplant data year 2 clean 190819.csv")

#######################Greenhouse Comp 2019 Cleaning##########################
#raw data import
raw<-read.csv("./data/Seedling comp data with harvest 190818.csv")

#Calculate "max neighbors"
N.counts<-raw[,c("N_1","N_3","N_4","N_5","N_6","N_7","N_8","N_9","N_10", #no N_2
                 "N_11","N_12","N_13","N_14","N_15","N_16","N_17")]
raw$N.max<-apply(N.counts,1,function(x){max(x[!is.na(x)])})
raw$N.max[which(raw$N.max=="-Inf")]<-0

#Label intra/inter/solo based on harvest neighbor data
raw$Type<-rep(0,180)
for (i in 1:180) {
  if (raw$FocalNew[i]!=raw$NeighborNew[i]) {
    raw$Type[i]<-"Inter"
    next
  }
  if (raw$FocalNew[i]==raw$NeighborNew[i]) {
    if(is.na(raw$N.all.mass[i])) raw$Type[i]<-"Solo"
    else raw$Type[i]<-"Intra"
  }
}

#Calculate neighbor per cap biomass
raw$N.all.percap<-raw$N.all.mass/(raw$N.dead+raw$N.live)
raw$N.live.percap<-raw$N.all.mass*raw$N.live/(raw$N.dead+raw$N.live)

#OUTPUT clean and useful data for analyses
clean<-raw[,c("ID","Rep","FocalNew","NeighborNew","F.H","N.live","N.dead","F.mass","N.all.mass","N.max","Type","N.all.percap","N.live.percap")]
write.csv(clean,"./data/Greenhouse comp data clean 190818.csv")

#######################Greenhouse Data 2018 Cleaning##########################
#Raw data import
raw<-read.csv("./data/Greenhouse data 190120.csv")

#Derive values of interest---------------

#Germination 
raw$harvested_extra[is.na(raw$harvested_extra)]<-0
#Germination should be calculated as the number of harvested extra +1 except for those that had no germination
#Make a column for germ yes/no
raw$germStatus[raw$g_180605==0]<-0
raw$germStatus[raw$g_180605!=0]<-1
#Calculate final germination
for (i in 1:length(raw[,1])){
  if (raw$germStatus[i]==0) raw$germ[i]<-0
  else raw$germ[i]<-max(raw$harvested_extra[i]+1,raw$g_180605[i])
}

#Survival
raw$newGerm[is.na(raw$newGerm)]<-0
for (i in 1:length(raw[,1])){
  if (raw$germStatus[i]!=0&raw$newGerm[i]!=1) raw$surv[i]<-raw$surv_180713[i]#don't count survival for new germs and no germs
  if (raw$germStatus[i]==0) raw$surv[i]<-NA
}

#Belowground biomass
#Examine dry-wet ratios
dryWet<-raw$blw.dry/raw$blw.sub
boxplot(dryWet~raw$Plant)#not super different
hist(dryWet)
summary.aov(lm(dryWet~raw$Plant))#sig diffs among spp for dry:wet
#Process data
raw$blw<-NA
for (i in 1:length(raw[,1])){
  if(!is.na(raw$blw.dry[i])&is.na(raw$blw.sub[i])) {
    raw$blw[i]<-raw$blw.dry[i]
    next
  }
  if (is.na(raw$blw.wet[i])&!is.na(raw$blw.sub[i])){#use sub dry for those with no full wet samples (n=2)
    raw$blw[i]<-raw$blw.dry[i]
    next
  }
  if (!is.na(raw$blw.sub[i])){
    raw$blw[i]<-raw$blw.wet[i]*raw$blw.dry[i]/raw$blw.sub[i]
  }
  
}

#Calculate RMR (root mass ratio)
raw$RMR<-raw$blw/(raw$abv+raw$blw)
#check distribution
hist(raw$RMR)

#OUTPUT clean data with only columns of interest-------------------
clean<-raw[,c(1:5,16,19,21,23:27)]
write.csv(clean,"./data/Greenhouse data clean 190122.csv")


#######################Field Data 2017/2018 Cleaning##########################
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
