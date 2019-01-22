##############################################
#Code to visualize greenhouse experiment data#
##############################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(car)
library(ggpubr)

#Data import------------------------
all<-read.csv("./data/Greenhouse data clean 190122.csv")

#Subset data by transplant species-----------
ARTR<-subset(all,Plant=="ARTR")
POSE<-subset(all,Plant=="POSE")
PSSP<-subset(all,Plant=="PSSP")
HECO<-subset(all,Plant=="HECO")

#ARTR--------------------------------
#Add Intra/Inter labels
ARTR$Type[ARTR$Source=="ARTR"]<-"Intra"
ARTR$Type[ARTR$Source!="ARTR"]<-"Inter"
ARTR$Type <- factor(ARTR$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#ARTR aboveground biomass visualization
#Make each donor plot separately
artrBARE<-ggpar(ggline(ARTR[ARTR$Source=="BARE"|ARTR$Source=="ARTR",], x = "Type", y = "abv",  
                 add = c("mean_se"),
                 color = "Inoc"),
                ylim=c(0,0.25))
artrHECO<-ggpar(ggline(ARTR[ARTR$Source=="HECO"|ARTR$Source=="ARTR",], x = "Type", y = "abv",  
                 add = c("mean_se"),
                 color = "Inoc"),
                ylim=c(0,0.25))
artrPSSP<-ggpar(ggline(ARTR[ARTR$Source=="PSSP"|ARTR$Source=="ARTR",], x = "Type", y = "abv",  
                 add = c("mean_se"),
                 color = "Inoc"),
                ylim=c(0,0.25))
artrPOSE<-ggpar(ggline(ARTR[ARTR$Source=="POSE"|ARTR$Source=="ARTR",], x = "Type", y = "abv",  
                 add = c("mean_se"),
                 color = "Inoc"),
                ylim=c(0,0.25))
#Put plots together
ggarrange(artrBARE , artrHECO, artrPSSP,artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#ARTR RMR visualization
#Make each donor plot separately
r.artrBARE<-ggpar(ggline(ARTR[ARTR$Source=="BARE"|ARTR$Source=="ARTR",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
r.artrHECO<-ggpar(ggline(ARTR[ARTR$Source=="HECO"|ARTR$Source=="ARTR",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
r.artrPSSP<-ggpar(ggline(ARTR[ARTR$Source=="PSSP"|ARTR$Source=="ARTR",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
r.artrPOSE<-ggpar(ggline(ARTR[ARTR$Source=="POSE"|ARTR$Source=="ARTR",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
#Put plots together
ggarrange(r.artrBARE , r.artrHECO, r.artrPSSP,r.artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#ARTR germination visualization
#Make each donor plot separately
g.artrBARE<-ggpar(ggline(ARTR[ARTR$Source=="BARE"|ARTR$Source=="ARTR",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,4))
g.artrHECO<-ggpar(ggline(ARTR[ARTR$Source=="HECO"|ARTR$Source=="ARTR",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,4))
g.artrPSSP<-ggpar(ggline(ARTR[ARTR$Source=="PSSP"|ARTR$Source=="ARTR",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,4))
g.artrPOSE<-ggpar(ggline(ARTR[ARTR$Source=="POSE"|ARTR$Source=="ARTR",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,4))
#Put plots together
ggarrange(g.artrBARE , g.artrHECO, g.artrPSSP,g.artrPOSE, 
          labels = c("Donor = BARE", "Donor = HECO", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)


#ARTR survival visualization
#Make each donor plot separately
s.artrBARE<-ggpar(ggerrorplot(ARTR[ARTR$Source=="BARE"|ARTR$Source=="ARTR",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.artrHECO<-ggpar(ggerrorplot(ARTR[ARTR$Source=="HECO"|ARTR$Source=="ARTR",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.artrPSSP<-ggpar(ggerrorplot(ARTR[ARTR$Source=="PSSP"|ARTR$Source=="ARTR",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.artrPOSE<-ggpar(ggerrorplot(ARTR[ARTR$Source=="POSE"|ARTR$Source=="ARTR",], x = "Type", y = "surv",  
                              desc_stat = "mean_se", 
                              color = "Inoc", 
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
HECO$Type[HECO$Source=="HECO"]<-"Intra"
HECO$Type[HECO$Source!="HECO"]<-"Inter"
HECO$Type <- factor(HECO$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#HECO aboveground biomass visualization
#Make each donor plot separately
hecoBARE<-ggpar(ggline(HECO[HECO$Source=="BARE"|HECO$Source=="HECO",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.5))
hecoARTR<-ggpar(ggline(HECO[HECO$Source=="ARTR"|HECO$Source=="HECO",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.5))
hecoPSSP<-ggpar(ggline(HECO[HECO$Source=="PSSP"|HECO$Source=="HECO",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.5))
hecoPOSE<-ggpar(ggline(HECO[HECO$Source=="POSE"|HECO$Source=="HECO",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.5))
#Put plots together
ggarrange(hecoBARE , hecoARTR, hecoPSSP,hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO RMR visualization
#Make each donor plot separately
r.hecoBARE<-ggpar(ggline(HECO[HECO$Source=="BARE"|HECO$Source=="HECO",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
r.hecoARTR<-ggpar(ggline(HECO[HECO$Source=="ARTR"|HECO$Source=="HECO",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
r.hecoPSSP<-ggpar(ggline(HECO[HECO$Source=="PSSP"|HECO$Source=="HECO",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
r.hecoPOSE<-ggpar(ggline(HECO[HECO$Source=="POSE"|HECO$Source=="HECO",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.75))
#Put plots together
ggarrange(r.hecoBARE , r.hecoARTR, r.hecoPSSP,r.hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO germination visualization
#Make each donor plot separately
g.hecoBARE<-ggpar(ggline(HECO[HECO$Source=="BARE"|HECO$Source=="HECO",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1.5))
g.hecoARTR<-ggpar(ggline(HECO[HECO$Source=="ARTR"|HECO$Source=="HECO",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1.5))
g.hecoPSSP<-ggpar(ggline(HECO[HECO$Source=="PSSP"|HECO$Source=="HECO",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1.5))
g.hecoPOSE<-ggpar(ggline(HECO[HECO$Source=="POSE"|HECO$Source=="HECO",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1.5))
#Put plots together
ggarrange(g.hecoBARE ,g.hecoARTR, g.hecoPSSP,g.hecoPOSE, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = POSE"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#HECO survival visualization
#Make each donor plot separately
s.hecoBARE<-ggpar(ggerrorplot(HECO[HECO$Source=="BARE"|HECO$Source=="HECO",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
            ylim=c(0,1.3))
s.hecoARTR<-ggpar(ggerrorplot(HECO[HECO$Source=="ARTR"|HECO$Source=="HECO",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.hecoPSSP<-ggpar(ggerrorplot(HECO[HECO$Source=="PSSP"|HECO$Source=="HECO",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.hecoPOSE<-ggpar(ggerrorplot(HECO[HECO$Source=="POSE"|HECO$Source=="HECO",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
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
POSE$Type[POSE$Source=="POSE"]<-"Intra"
POSE$Type[POSE$Source!="POSE"]<-"Inter"
POSE$Type <- factor(POSE$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#POSE aboveground biomass visualization
#Make each donor plot separately
poseBARE<-ggpar(ggline(POSE[POSE$Source=="BARE"|POSE$Source=="POSE",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.2))
poseARTR<-ggpar(ggline(POSE[POSE$Source=="ARTR"|POSE$Source=="POSE",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.2))
posePSSP<-ggpar(ggline(POSE[POSE$Source=="PSSP"|POSE$Source=="POSE",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.2))
poseHECO<-ggpar(ggline(POSE[POSE$Source=="HECO"|POSE$Source=="POSE",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.2))
#Put plots together
ggarrange(poseBARE , poseARTR, posePSSP,poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE RMR visualization
#Make each donor plot separately
r.poseBARE<-ggpar(ggline(POSE[POSE$Source=="BARE"|POSE$Source=="POSE",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
r.poseARTR<-ggpar(ggline(POSE[POSE$Source=="ARTR"|POSE$Source=="POSE",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
r.posePSSP<-ggpar(ggline(POSE[POSE$Source=="PSSP"|POSE$Source=="POSE",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
r.poseHECO<-ggpar(ggline(POSE[POSE$Source=="HECO"|POSE$Source=="POSE",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
#Put plots together
ggarrange(r.poseBARE , r.poseARTR, r.posePSSP,r.poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE germination visualization
#Make each donor plot separately
g.poseBARE<-ggpar(ggline(POSE[POSE$Source=="BARE"|POSE$Source=="POSE",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,7))
g.poseARTR<-ggpar(ggline(POSE[POSE$Source=="ARTR"|POSE$Source=="POSE",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,7))
g.posePSSP<-ggpar(ggline(POSE[POSE$Source=="PSSP"|POSE$Source=="POSE",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,7))
g.poseHECO<-ggpar(ggline(POSE[POSE$Source=="HECO"|POSE$Source=="POSE",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,7))
#Put plots together
ggarrange(g.poseBARE , g.poseARTR, g.posePSSP,g.poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#POSE survival visualization
#Make each donor plot separately
s.poseBARE<-ggpar(ggline(POSE[POSE$Source=="BARE"|POSE$Source=="POSE",], x = "Type", y = "surv",  
                         add = c("mean_se"),
                         color = "Inoc"),
                  ylim=c(0,1.3))
s.poseARTR<-ggpar(ggline(POSE[POSE$Source=="ARTR"|POSE$Source=="POSE",], x = "Type", y = "surv",  
                         add = c("mean_se"),
                         color = "Inoc"),
                  ylim=c(0,1.3))
s.posePSSP<-ggpar(ggline(POSE[POSE$Source=="PSSP"|POSE$Source=="POSE",], x = "Type", y = "surv",  
                         add = c("mean_se"),
                         color = "Inoc"),
                  ylim=c(0,1.3))
s.poseHECO<-ggpar(ggline(POSE[POSE$Source=="HECO"|POSE$Source=="POSE",], x = "Type", y = "surv",  
                         add = c("mean_se"),
                         color = "Inoc"),
                  ylim=c(0,1.3))
#Put plots together
ggarrange(s.poseBARE , s.poseARTR, s.posePSSP,s.poseHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = PSSP","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP--------------------------------
#Add Intra/Inter labels
PSSP$Type[PSSP$Source=="PSSP"]<-"Intra"
PSSP$Type[PSSP$Source!="PSSP"]<-"Inter"
PSSP$Type <- factor(PSSP$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization

#PSSP aboveground biommas visualization
#Make each donor plot separately
psspBARE<-ggpar(ggline(PSSP[PSSP$Source=="BARE"|PSSP$Source=="PSSP",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.4))
psspARTR<-ggpar(ggline(PSSP[PSSP$Source=="ARTR"|PSSP$Source=="PSSP",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.4))
psspPOSE<-ggpar(ggline(PSSP[PSSP$Source=="POSE"|PSSP$Source=="PSSP",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.4))
psspHECO<-ggpar(ggline(PSSP[PSSP$Source=="HECO"|PSSP$Source=="PSSP",], x = "Type", y = "abv",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,0.4))
#Put plots together
ggarrange(psspBARE , psspARTR, psspPOSE,psspHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = POSE","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP RMR visualization
#Make each donor plot separately
r.psspBARE<-ggpar(ggline(PSSP[PSSP$Source=="BARE"|PSSP$Source=="PSSP",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
r.psspARTR<-ggpar(ggline(PSSP[PSSP$Source=="ARTR"|PSSP$Source=="PSSP",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
r.psspPOSE<-ggpar(ggline(PSSP[PSSP$Source=="POSE"|PSSP$Source=="PSSP",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
r.psspHECO<-ggpar(ggline(PSSP[PSSP$Source=="HECO"|PSSP$Source=="PSSP",], x = "Type", y = "RMR",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,1))
#Put plots together
ggarrange(r.psspBARE , r.psspARTR, r.psspPOSE,r.psspHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = POSE","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP germination visualization
#Make each donor plot separately
g.psspBARE<-ggpar(ggline(PSSP[PSSP$Source=="BARE"|PSSP$Source=="PSSP",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,8))
g.psspARTR<-ggpar(ggline(PSSP[PSSP$Source=="ARTR"|PSSP$Source=="PSSP",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,8))
g.psspPOSE<-ggpar(ggline(PSSP[PSSP$Source=="POSE"|PSSP$Source=="PSSP",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,8))
g.psspHECO<-ggpar(ggline(PSSP[PSSP$Source=="HECO"|PSSP$Source=="PSSP",], x = "Type", y = "germ",  
                       add = c("mean_se"),
                       color = "Inoc"),
                ylim=c(0,8))
#Put plots together
ggarrange(g.psspBARE , g.psspARTR, g.psspPOSE,g.psspHECO, 
          labels = c("Donor = BARE", "Donor = ARTR", "Donor = POSE","Donor = HECO"),
          common.legend = TRUE, legend = "bottom",
          hjust=-0.7,
          ncol = 2, nrow = 2)

#PSSP survival visualization
#Make each donor plot separately
s.psspBARE<-ggpar(ggerrorplot(PSSP[PSSP$Source=="BARE"|PSSP$Source=="PSSP",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.psspARTR<-ggpar(ggerrorplot(PSSP[PSSP$Source=="ARTR"|PSSP$Source=="PSSP",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.psspPOSE<-ggpar(ggerrorplot(PSSP[PSSP$Source=="POSE"|PSSP$Source=="PSSP",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
                              position = position_dodge(0.3)),
                  ylim=c(0,1.3))
s.psspHECO<-ggpar(ggerrorplot(PSSP[PSSP$Source=="HECO"|PSSP$Source=="PSSP",], x = "Type", y = "surv", 
                              desc_stat = "mean_se", 
                              color = "Inoc", 
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

abv<-ggline(allClean, x = "Source", y = "abv",
       add = c("mean_se"),
       color = "Inoc")
rmr<-ggline(allClean, x = "Source", y = "RMR",
            add = c("mean_se"),
            color = "Inoc")
surv<-ggline(allClean, x = "Source", y = "surv",
            add = c("mean_se"),
            color = "Inoc")
germ<-ggline(allClean, x = "Source", y = "germ",
            add = c("mean_se"),
            color = "Inoc")
#Put plots together
ggarrange(abv ,rmr, surv, germ,
          common.legend = TRUE, legend = "bottom",
          ncol = 1, nrow = 4)
