########################################
#Code to visualize all abv biomass data#
########################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)
library(ggpubr)

#Data import------------------------
all<-read.csv("./data/Transplant data clean 190102.csv")

#Subset data by transplant species
ARTR<-subset(all,Transplant=="ARTR")
POSE<-subset(all,Transplant=="POSE")
PSSP<-subset(all,Transplant=="PSSP")
HECO<-subset(all,Transplant=="HECO")

#ARTR--------------------------------
#Add Intra/Inter labels
ARTR$Type[ARTR$DonorSpp=="ARTR"]<-"Intra"
ARTR$Type[ARTR$DonorSpp!="ARTR"]<-"Inter"
ARTR$Type <- factor(ARTR$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#Take out the rows where labels mixed up and cannot ascertain treatment
ARTR<-subset(ARTR,Treatment%in%c("Control","Exclusion","Feedback"))

#Make each donor plot separately
artrBARE<-ggpar(ggline(ARTR[ARTR$DonorSpp=="BARE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.25))
artrHECO<-ggpar(ggline(ARTR[ARTR$DonorSpp=="HECO"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.25))
artrPSSP<-ggpar(ggline(ARTR[ARTR$DonorSpp=="PSSP"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.25))
artrPOSE<-ggpar(ggline(ARTR[ARTR$DonorSpp=="POSE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.25))

#Put plots together
ggarrange(artrBARE , artrHECO, artrPSSP,artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO--------------------------------
#Add Intra/Inter labels
HECO$Type[HECO$DonorSpp=="HECO"]<-"Intra"
HECO$Type[HECO$DonorSpp!="HECO"]<-"Inter"
HECO$Type <- factor(HECO$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#Make each donor plot separately
hecoBARE<-ggpar(ggline(HECO[HECO$DonorSpp=="BARE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.25))
hecoARTR<-ggpar(ggline(HECO[HECO$DonorSpp=="ARTR"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.25))
hecoPSSP<-ggpar(ggline(HECO[HECO$DonorSpp=="PSSP"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.25))
hecoPOSE<-ggpar(ggline(HECO[HECO$DonorSpp=="POSE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.25))

#Put plots together
ggarrange(hecoBARE , hecoARTR, hecoPSSP,hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE--------------------------------
#Add Intra/Inter labels
POSE$Type[POSE$DonorSpp=="POSE"]<-"Intra"
POSE$Type[POSE$DonorSpp!="POSE"]<-"Inter"
POSE$Type <- factor(POSE$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#Make each donor plot separately
poseBARE<-ggpar(ggline(POSE[POSE$DonorSpp=="BARE"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))
poseARTR<-ggpar(ggline(POSE[POSE$DonorSpp=="ARTR"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))
posePSSP<-ggpar(ggline(POSE[POSE$DonorSpp=="PSSP"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))
poseHECO<-ggpar(ggline(POSE[POSE$DonorSpp=="HECO"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))

#Put plots together
ggarrange(poseBARE , poseARTR, posePSSP,poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP--------------------------------
#Add Intra/Inter labels
PSSP$Type[PSSP$DonorSpp=="PSSP"]<-"Intra"
PSSP$Type[PSSP$DonorSpp!="PSSP"]<-"Inter"
PSSP$Type <- factor(PSSP$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#Make each donor plot separately
psspBARE<-ggpar(ggline(PSSP[PSSP$DonorSpp=="BARE"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))
psspARTR<-ggpar(ggline(PSSP[PSSP$DonorSpp=="ARTR"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))
psspPOSE<-ggpar(ggline(PSSP[PSSP$DonorSpp=="POSE"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))
psspHECO<-ggpar(ggline(PSSP[PSSP$DonorSpp=="HECO"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.4))

#Put plots together
ggarrange(psspBARE , psspARTR, psspPOSE,psspHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = POSE","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)
