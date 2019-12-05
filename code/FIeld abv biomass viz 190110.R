########################################
#Code to visualize all abv biomass data#
########################################

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)
library(ggpubr)

#Data import------------------------
#all<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Subset data by transplant species-----------
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

#ARTR aboveground biomass visualization
#Make each donor plot separately
artrBARE<-ggpar(ggline(ARTR[ARTR$DonorSpp=="BARE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.15))#use 0.25 for 2018 data, 0.15 for 2019 data
artrHECO<-ggpar(ggline(ARTR[ARTR$DonorSpp=="HECO"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.15))
artrPSSP<-ggpar(ggline(ARTR[ARTR$DonorSpp=="PSSP"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.15))
artrPOSE<-ggpar(ggline(ARTR[ARTR$DonorSpp=="POSE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Abv",  
                 add = c("mean_se"),
                 color = "Treatment"),
                ylim=c(0,0.15))
#Put plots together
ggarrange(artrBARE , artrHECO, artrPSSP,artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#ARTR germination visualization
#Make each donor plot separately
g.artrBARE<-ggpar(ggline(ARTR[ARTR$DonorSpp=="BARE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.artrHECO<-ggpar(ggline(ARTR[ARTR$DonorSpp=="HECO"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.artrPSSP<-ggpar(ggline(ARTR[ARTR$DonorSpp=="PSSP"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.artrPOSE<-ggpar(ggline(ARTR[ARTR$DonorSpp=="POSE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
#Put plots together
ggarrange(g.artrBARE , g.artrHECO, g.artrPSSP,g.artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)


#ARTR survival visualization
#Make each donor plot separately
s.artrBARE<-ggpar(ggerrorplot(ARTR[ARTR$DonorSpp=="BARE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.artrHECO<-ggpar(ggerrorplot(ARTR[ARTR$DonorSpp=="HECO"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.artrPSSP<-ggpar(ggerrorplot(ARTR[ARTR$DonorSpp=="PSSP"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.artrPOSE<-ggpar(ggerrorplot(ARTR[ARTR$DonorSpp=="POSE"|ARTR$DonorSpp=="ARTR",], x = "Type", y = "Surv",  
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
#Put plots together
ggarrange(s.artrBARE , s.artrHECO, s.artrPSSP,s.artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO--------------------------------
#Add Intra/Inter labels
HECO$Type[HECO$DonorSpp=="HECO"]<-"Intra"
HECO$Type[HECO$DonorSpp!="HECO"]<-"Inter"
HECO$Type <- factor(HECO$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#HECO aboveground biomass visualization
#Make each donor plot separately
hecoBARE<-ggpar(ggline(HECO[HECO$DonorSpp=="BARE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.15))#use 0.25 for 2018 and 0.15 for 2019
hecoARTR<-ggpar(ggline(HECO[HECO$DonorSpp=="ARTR"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.15))
hecoPSSP<-ggpar(ggline(HECO[HECO$DonorSpp=="PSSP"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.15))
hecoPOSE<-ggpar(ggline(HECO[HECO$DonorSpp=="POSE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.15))
#Put plots together
ggarrange(hecoBARE , hecoARTR, hecoPSSP,hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO germination visualization 
#Make each donor plot separately
g.hecoBARE<-ggpar(ggline(HECO[HECO$DonorSpp=="BARE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))#use 1 for 2018, 6 for 2019
g.hecoARTR<-ggpar(ggline(HECO[HECO$DonorSpp=="ARTR"|HECO$DonorSpp=="HECO",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.hecoPSSP<-ggpar(ggline(HECO[HECO$DonorSpp=="PSSP"|HECO$DonorSpp=="HECO",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.hecoPOSE<-ggpar(ggline(HECO[HECO$DonorSpp=="POSE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
#Put plots together
ggarrange(g.hecoBARE ,g.hecoARTR, g.hecoPSSP,g.hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO survival visualization
#Make each donor plot separately
s.hecoBARE<-ggpar(ggerrorplot(HECO[HECO$DonorSpp=="BARE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
            ylim=c(0,1.3))
s.hecoARTR<-ggpar(ggerrorplot(HECO[HECO$DonorSpp=="ARTR"|HECO$DonorSpp=="HECO",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.hecoPSSP<-ggpar(ggerrorplot(HECO[HECO$DonorSpp=="PSSP"|HECO$DonorSpp=="HECO",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.hecoPOSE<-ggpar(ggerrorplot(HECO[HECO$DonorSpp=="POSE"|HECO$DonorSpp=="HECO",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
#Put plots together
ggarrange(s.hecoBARE ,s.hecoARTR, s.hecoPSSP,s.hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE--------------------------------
#Add Intra/Inter labels
POSE$Type[POSE$DonorSpp=="POSE"]<-"Intra"
POSE$Type[POSE$DonorSpp!="POSE"]<-"Inter"
POSE$Type <- factor(POSE$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#POSE aboveground biomass visualization
#Make each donor plot separately
poseBARE<-ggpar(ggline(POSE[POSE$DonorSpp=="BARE"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.3))#use 0.4 in 2018 and 0.3 in 2019
poseARTR<-ggpar(ggline(POSE[POSE$DonorSpp=="ARTR"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.3))
posePSSP<-ggpar(ggline(POSE[POSE$DonorSpp=="PSSP"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.3))
poseHECO<-ggpar(ggline(POSE[POSE$DonorSpp=="HECO"|POSE$DonorSpp=="POSE",], x = "Type", y = "Abv",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,0.3))
#Put plots together
ggarrange(poseBARE , poseARTR, posePSSP,poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE germination visualization
#Make each donor plot separately
g.poseBARE<-ggpar(ggline(POSE[POSE$DonorSpp=="BARE"|POSE$DonorSpp=="POSE",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))#use 4 for 2018, 6 for 2019
g.poseARTR<-ggpar(ggline(POSE[POSE$DonorSpp=="ARTR"|POSE$DonorSpp=="POSE",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.posePSSP<-ggpar(ggline(POSE[POSE$DonorSpp=="PSSP"|POSE$DonorSpp=="POSE",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
g.poseHECO<-ggpar(ggline(POSE[POSE$DonorSpp=="HECO"|POSE$DonorSpp=="POSE",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,6))
#Put plots together
ggarrange(g.poseBARE , g.poseARTR, g.posePSSP,g.poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE survival visualization for 2018
#Make each donor plot separately
s.poseBARE<-ggpar(ggline(POSE[POSE$DonorSpp=="BARE"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv",  
                         add = c("mean_se"),
                         color = "Treatment"),
                  ylim=c(0,1.3))
s.poseARTR<-ggpar(ggline(POSE[POSE$DonorSpp=="ARTR"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv",  
                         add = c("mean_se"),
                         color = "Treatment"),
                  ylim=c(0,1.3))
s.posePSSP<-ggpar(ggline(POSE[POSE$DonorSpp=="PSSP"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv",  
                         add = c("mean_se"),
                         color = "Treatment"),
                  ylim=c(0,1.3))
s.poseHECO<-ggpar(ggline(POSE[POSE$DonorSpp=="HECO"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv",  
                         add = c("mean_se"),
                         color = "Treatment"),
                  ylim=c(0,1.3))
#Put plots together
ggarrange(s.poseBARE , s.poseARTR, s.posePSSP,s.poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE survival visualization for 2019
#Make each donor plot separately
s.poseBARE<-ggpar(ggerrorplot(POSE[POSE$DonorSpp=="BARE"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.poseARTR<-ggpar(ggerrorplot(POSE[POSE$DonorSpp=="ARTR"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.posePSSP<-ggpar(ggerrorplot(POSE[POSE$DonorSpp=="PSSP"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.poseHECO<-ggpar(ggerrorplot(POSE[POSE$DonorSpp=="HECO"|POSE$DonorSpp=="POSE",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
#Put plots together
ggarrange(s.poseBARE , s.poseARTR, s.posePSSP,s.poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP--------------------------------
#Add Intra/Inter labels
PSSP$Type[PSSP$DonorSpp=="PSSP"]<-"Intra"
PSSP$Type[PSSP$DonorSpp!="PSSP"]<-"Inter"
PSSP$Type <- factor(PSSP$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#PSSP aboveground biommas visualization
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

#PSSP germination visualization
#Make each donor plot separately
g.psspBARE<-ggpar(ggline(PSSP[PSSP$DonorSpp=="BARE"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,7))#use 3 in 2018, 7 in 2019
g.psspARTR<-ggpar(ggline(PSSP[PSSP$DonorSpp=="ARTR"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,7))
g.psspPOSE<-ggpar(ggline(PSSP[PSSP$DonorSpp=="POSE"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,7))
g.psspHECO<-ggpar(ggline(PSSP[PSSP$DonorSpp=="HECO"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Germ",  
                       add = c("mean_se"),
                       color = "Treatment"),
                ylim=c(0,7))
#Put plots together
ggarrange(g.psspBARE , g.psspARTR, g.psspPOSE,g.psspHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = POSE","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP survival visualization
#Make each donor plot separately
s.psspBARE<-ggpar(ggerrorplot(PSSP[PSSP$DonorSpp=="BARE"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.psspARTR<-ggpar(ggerrorplot(PSSP[PSSP$DonorSpp=="ARTR"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.psspPOSE<-ggpar(ggerrorplot(PSSP[PSSP$DonorSpp=="POSE"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.psspHECO<-ggpar(ggerrorplot(PSSP[PSSP$DonorSpp=="HECO"|PSSP$DonorSpp=="PSSP",], x = "Type", y = "Surv", 
                              desc_stat = "mean_se", 
                              color = "Treatment", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
#Put plots together
ggarrange(s.psspBARE , s.psspARTR, s.psspPOSE,s.psspHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = POSE","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#Look at species together--------------------------------
allClean<-rbind.data.frame(ARTR,POSE,HECO,PSSP)

abv<-ggline(allClean, x = "DonorSpp", y = "Abv",
       add = c("mean_se"),
       color = "Treatment")
surv<-ggline(allClean, x = "DonorSpp", y = "Surv",
            add = c("mean_se"),
            color = "Treatment")
germ<-ggline(allClean, x = "DonorSpp", y = "Germ",
            add = c("mean_se"),
            color = "Treatment")
#Put plots together
ggarrange(abv , surv, germ,
          common.legend = TRUE, legend = "bottom",
          ncol = 1, nrow = 3)

#2019 only, use distance to donor and donor size as covariate---------------
#Visualize everything
ggscatter(allClean,x="Dist_centroid",y="Abv",color="Treatment",shape="DonorSpp",facet.by="Transplant")
#Big blob no discernable patterns
ggscatter(allClean,x="Diam_donor",y="Abv",color="Treatment",shape="DonorSpp",facet.by="Transplant")
ggscatter(allClean,x="Dist_edge",y="Abv",color="Treatment",shape="DonorSpp",facet.by="Transplant") 
#easier to visually compare since it gets rid of the large size of ARTRs, no strong patterns though

#Separate by transplant species
#ARTR
ggscatter(subset(ARTR,DonorSpp!="BARE"),x="Dist_centroid",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(ARTR,DonorSpp!="BARE"),x="Dist_edge",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(ARTR,DonorSpp!="BARE"),x="Diam_donor",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
#maybe a negative relationship in control POSEs

#HECO
ggscatter(subset(HECO,DonorSpp!="BARE"),x="Dist_centroid",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
#In HECO feedback, potentially a positive trend
ggscatter(subset(HECO,DonorSpp!="BARE"),x="Dist_edge",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(HECO,DonorSpp!="BARE"),x="Diam_donor",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
#In HECO feedback, potentially a positive trend

#POSE
ggscatter(subset(POSE,DonorSpp!="BARE"),x="Dist_centroid",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(POSE,DonorSpp!="BARE"),x="Dist_edge",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(POSE,DonorSpp!="BARE"),x="Diam_donor",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")

#PSSP
ggscatter(subset(PSSP,DonorSpp!="BARE"),x="Dist_centroid",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(PSSP,DonorSpp!="BARE"),x="Dist_edge",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
ggscatter(subset(PSSP,DonorSpp!="BARE"),x="Diam_donor",y="Abv",color="Treatment",facet.by="DonorSpp",scales="free")
