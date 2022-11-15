#Code to try and combine all life stage transitions in field PSF experiment
#AC Jun 2022

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
library(ggpubr)
library(ggplot2)
library(tidyr)
library(cowplot)

#Data import------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data
#Take out the rows where labels mixed up and cannot ascertain treatment
all2018<-subset(all2018,Treatment%in%c("Control","Exclusion","Feedback"))

#Estimating params--------
#Add unique treatment labels
all2018$Trans.Don.Treat<-paste(all2018$Transplant,all2018$DonorSpp,all2018$Treatment,sep=".")
all2019$Trans.Don.Treat<-paste(all2019$Transplant,all2019$DonorSpp,all2019$Treatment,sep=".")

#Bootstrap params
source("./code/func_BSDemog.R")
mass.boot.2018.list<-bsDemog(dat=all2018,nb=5000,bs.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
mass.boot.2019.list<-bsDemog(dat=all2019,nb=5000,bs.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")

#Bootstrap the other way (same results)
#source("./code/func_BS1Demog.R")
#mass.boot.2018.list<-bsDemog1(dat=all2018,nb=5000,bs.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
#mass.boot.2019.list<-bsDemog1(dat=all2019,nb=5000,bs.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")

#Assign estimate dataframe
mass.boot.2018<-mass.boot.2018.list[[1]]
mass.boot.2019<-mass.boot.2019.list[[1]]

#Assign bootstrap vector lists
bsVectors2018<-mass.boot.2018.list[[2]]
bsVectors2019<-mass.boot.2019.list[[2]]
names(bsVectors2018)<-mass.boot.2018$treat.combo
names(bsVectors2019)<-mass.boot.2019$treat.combo

# #Non-bootstrap
# source("./code/func_noBSDemog.R")
# mass.2018<-nobsDemog(dat=all2018,levels.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
# mass.2019<-nobsDemog(dat=all2019,levels.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")

#Data clean up and adding labels-----------
#Delete final row from 2018 data since it has no treatment label
mass.boot.2018<-mass.boot.2018[1:60,]
#mass.2018<-mass.2018[1:60,]

#Add labels back to data
mass.boot.2018<-mass.boot.2018 %>% separate(treat.combo, c("Transplant", "Donorspp","Treatment"))
mass.boot.2019<-mass.boot.2019 %>% separate(treat.combo, c("Transplant", "Donorspp","Treatment"))


#Make Intra/Inter label
mass.boot.2018$Type[which(mass.boot.2018$Donorspp!=mass.boot.2018$Transplant)]<-"Inter"
mass.boot.2018$Type[which(mass.boot.2018$Donorspp==mass.boot.2018$Transplant)]<-"Intra"
mass.boot.2018$Type <- factor(mass.boot.2018$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
mass.boot.2018$Treatment <- factor(mass.boot.2018$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization
mass.boot.2019$Type[which(mass.boot.2019$Donorspp!=mass.boot.2019$Transplant)]<-"Inter"
mass.boot.2019$Type[which(mass.boot.2019$Donorspp==mass.boot.2019$Transplant)]<-"Intra"
mass.boot.2019$Type <- factor(mass.boot.2019$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
mass.boot.2019$Treatment <- factor(mass.boot.2019$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization

#Visualize projections-----------------
#BOOTSTRAPPED
#2018
p.demog18<-ggplot(data=mass.boot.2018, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
      facet_wrap(~Transplant,nrow=1)+
      geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                      , position = position_dodge(width = 1))+  
      scale_shape_discrete(name="Microsite")+
      scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
      xlab("")+ylab("")+
      geom_vline(xintercept=1.5,color="black",linetype="dashed")+
      theme_bw()+
      theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#2019
p.demog19<-ggplot(data=mass.boot.2019, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
          facet_wrap(~Transplant,nrow=1)+
          geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                          , position = position_dodge(width = 1))+  
          scale_shape_discrete(name="Microsite")+ 
          scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
          xlab("Soil environment")+ylab("")+
          geom_vline(xintercept=1.5,color="black",linetype="dashed")+
          theme_bw()+
          theme(legend.position = "none",
            panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#Put together in panel
annotate_figure(ggarrange(p.demog18,p.demog19,ncol = 1,labels = c("A","B"),vjust=T,common.legend = TRUE, legend = "right"),
                left = text_grob("Estimated shoot biomass contribution from 100 seeds (g)", rot = 90))     

#Alternative visualization 
ggplot(data=mass.boot.2018, aes(x=Donorspp, y=massContrib.mean, color=Transplant, shape=Transplant))+
  facet_grid(~Treatment)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 0.5))+
  geom_line(aes(group=Transplant), position=position_dodge(width=0.5))+
  scale_color_manual(values=c("#E2D200", "#46ACC8", "#E58601","#00A08A"))+
  xlab("Soil environment microsite")+ylab("Estimated shoot biomass contribution from 100 seeds (g)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggplot(data=mass.boot.2019, aes(x=Donorspp, y=massContrib.mean, color=Transplant, shape=Transplant))+
  facet_grid(~Treatment)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 0.5))+
  geom_line(aes(group=Transplant), position=position_dodge(width=0.5))+
  scale_color_manual(values=c("#E2D200", "#46ACC8", "#E58601","#00A08A"))+
  xlab("Soil environment microsite")+ylab("Estimated shoot biomass contribution from 100 seeds (g)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#New figures for panel plot with better annotations and spacing------------------
#2018
mass.boot.2018$plotX<-paste(mass.boot.2018$Type,mass.boot.2018$Donorspp)
mass.boot.2018$plotX <- factor(mass.boot.2018$plotX, 
                            levels=c("Intra ARTR", "Intra HECO","Intra POSE","Intra PSSP",
                                     "Inter ARTR", "Inter BARE", "Inter HECO", "Inter POSE", "Inter PSSP"))#Flip the order of factors for easier visualization
demog18<-ggplot(data=mass.boot.2018, aes(x=plotX, y=massContrib.mean, color=Treatment, shape=Donorspp))+
  facet_grid(~Transplant,scales="free_x")+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("")+ylab("")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
p.demog18<-demog18+coord_cartesian(clip = "off") + # allows plotting anywhere on the canvas
  draw_label("Intra", x = 1, y = -1.1, size=10)+
  draw_label("Inter", x = 3.5, y = -1.1, size=10) 


#2019
mass.boot.2019$plotX<-paste(mass.boot.2019$Type,mass.boot.2019$Donorspp)
mass.boot.2019$plotX <- factor(mass.boot.2019$plotX, 
                               levels=c("Intra ARTR", "Intra HECO","Intra POSE","Intra PSSP",
                                        "Inter ARTR", "Inter BARE", "Inter HECO", "Inter POSE", "Inter PSSP"))#Flip the order of factors for easier visualization
demog19<-ggplot(data=mass.boot.2019, aes(x=plotX, y=massContrib.mean, color=Treatment, shape=Donorspp))+
  facet_grid(~Transplant,scales="free_x")+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("\nSoil environment")+ylab("")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
p.demog19<-demog19+coord_cartesian(clip = "off") + # allows plotting anywhere on the canvas
  draw_label("Intra", x = 1, y = -1.8, size=10)+
  draw_label("Inter", x = 3.5, y = -1.8, size=10)#not sure why this x value has to be different 

#Put together in panel
annotate_figure(ggarrange(p.demog18,p.demog19,ncol = 1,labels = c("A","B"),vjust=T,common.legend = TRUE, legend = "right"),
                left = text_grob("Estimated shoot biomass contribution from 100 seeds (g)", rot = 90))     


# #NONBOOTSTRAPPED
# #2018 using SEs
# ggplot(data=mass.2018, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
#   facet_wrap(~Transplant,nrow=1)+
#   geom_pointrange(aes(ymin=massContrib.mean-massContrib.SE, ymax=massContrib.mean+massContrib.SE)
#                   , position = position_dodge(width = 1))+  #95%CI is huge
#   scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
#   xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
#   theme_bw()
# #2019
# ggplot(data=mass.2019, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
#   facet_wrap(~Transplant,nrow=1)+
#   geom_pointrange(aes(ymin=massContrib.mean-massContrib.SE, ymax=massContrib.mean+massContrib.SE)
#                   , position = position_dodge(width = 1))+  #Errors are giant now and very asymmetric
#   scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
#   xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
#   theme_bw()
# #Alternative visualization (2018 example)
# ggplot(data=mass.2018, aes(x=Donorspp, y=massContrib.mean, color=Transplant, shape=Transplant))+
#   facet_grid(~Treatment)+
#   geom_pointrange(aes(ymin=massContrib.mean-massContrib.SE, ymax=massContrib.mean+massContrib.SE)
#                   , position = position_dodge(width = 0.5))+
#   geom_line(aes(group=Transplant), position=position_dodge(width=0.5))+
#   scale_color_manual(values=c("#E2D200", "#46ACC8", "#E58601","#00A08A"))+
#   xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds (g)")+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#Visualize bootrstrapped Years combined Controls only--------------------
#2018
control2018<-ggplot(data=mass.boot.2018[mass.boot.2018$Treatment=="Control",], aes(x=Type, y=massContrib.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("")+ylab("2018 \nEst. shoot mass per 100 seeds (g)")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#2019
control2019<-ggplot(data=mass.boot.2019[mass.boot.2019$Treatment=="Control",], aes(x=Type, y=massContrib.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Neighbor type")+ylab("2019 \nEst. shoot mass per 100 seeds (g)")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#put together
ggarrange(control2018, control2019, 
          ncol=1, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

#Germination only
#2018
control2018_germ<-ggplot(data=mass.boot.2018[mass.boot.2018$Treatment=="Control",], aes(x=Type, y=germ.boot.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=germ.boot.95CIlo, ymax=germ.boot.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("")+ylab("2018 \nbootsrap germination")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#2019
control2019_germ<-ggplot(data=mass.boot.2019[mass.boot.2019$Treatment=="Control",], aes(x=Type, y=germ.boot.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=germ.boot.95CIlo, ymax=germ.boot.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Neighbor type")+ylab("2019 \nEst.bootsrap germination")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#put together
ggarrange(control2018_germ, control2019_germ, 
          ncol=1, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

#Survival only
#2018
control2018_surv<-ggplot(data=mass.boot.2018[mass.boot.2018$Treatment=="Control",], aes(x=Type, y=surv.boot.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=surv.boot.95CIlo, ymax=surv.boot.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("")+ylab("2018 \nbootsrap survival")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#2019
control2019_surv<-ggplot(data=mass.boot.2019[mass.boot.2019$Treatment=="Control",], aes(x=Type, y=surv.boot.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=surv.boot.95CIlo, ymax=surv.boot.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Neighbor type")+ylab("2019 \nEst.bootsrap survival")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#put together
ggarrange(control2018_surv, control2019_surv, 
          ncol=1, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

#Biomass only
#2018
control2018_mass<-ggplot(data=mass.boot.2018[mass.boot.2018$Treatment=="Control",], aes(x=Type, y=mass.boot.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=mass.boot.95CIlo, ymax=mass.boot.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("")+ylab("2018 \nbootsrap biomass")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#2019
control2019_mass<-ggplot(data=mass.boot.2019[mass.boot.2019$Treatment=="Control",], aes(x=Type, y=mass.boot.mean, color=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=mass.boot.95CIlo, ymax=mass.boot.95CIhi)
                  , position = position_dodge(width = 1))+  
  #scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Neighbor type")+ylab("2019 \nEst.bootsrap biomass")+
  guides(color=guide_legend("Neighbor spp."))+
  theme_bw()
#put together
ggarrange(control2018_mass, control2019_mass, 
          ncol=1, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

#Permutational differences----------------------------
#We want to see if the 95% CI of all other vs. self comparisons overlaps zero
#2018
bs.diff2018<-data.frame(homeSpp=rep(c("ARTR","HECO","POSE","PSSP"),each=12),
                   awaySpp=c(rep(c("HECO","POSE","PSSP","BARE"),3),
                             rep(c("ARTR","POSE","PSSP","BARE"),3),
                             rep(c("ARTR","HECO","PSSP","BARE"),3),
                             rep(c("ARTR","HECO","POSE","BARE"),3)),
                   Treatment=rep(c("Control","Feedback","Exclusion"),16),
                   diff.mean=rep(0,48),
                   diff.95lo=rep(0,48),
                   diff.95hi=rep(0,48))
for (i in 1:48){
  home<-bsVectors2018[[which(mass.boot.2018$Transplant==bs.diff2018$homeSpp[i]&mass.boot.2018$Donorspp==bs.diff2018$homeSpp[i]&mass.boot.2018$Treatment==bs.diff2018$Treatment[i])]]
  away<-bsVectors2018[[which(mass.boot.2018$Transplant==bs.diff2018$homeSpp[i]&mass.boot.2018$Donorspp==bs.diff2018$awaySpp[i]&mass.boot.2018$Treatment==bs.diff2018$Treatment[i])]]
  diff<-home-away
  bs.diff2018$diff.mean[i]<-mean(diff)
  bs.diff2018$diff.95lo[i]<-quantile(diff,probs=0.025)
  bs.diff2018$diff.95hi[i]<-quantile(diff,probs=0.975)
}

#2019
bs.diff2019<-data.frame(homeSpp=rep(c("ARTR","HECO","POSE","PSSP"),each=12),
                        awaySpp=c(rep(c("HECO","POSE","PSSP","BARE"),3),
                                  rep(c("ARTR","POSE","PSSP","BARE"),3),
                                  rep(c("ARTR","HECO","PSSP","BARE"),3),
                                  rep(c("ARTR","HECO","POSE","BARE"),3)),
                        Treatment=rep(c("Control","Feedback","Exclusion"),16),
                        diff.mean=rep(0,48),
                        diff.95lo=rep(0,48),
                        diff.95hi=rep(0,48))
for (i in 1:48){
  home<-bsVectors2019[[which(mass.boot.2019$Transplant==bs.diff2019$homeSpp[i]&mass.boot.2019$Donorspp==bs.diff2019$homeSpp[i]&mass.boot.2019$Treatment==bs.diff2019$Treatment[i])]]
  away<-bsVectors2019[[which(mass.boot.2019$Transplant==bs.diff2019$homeSpp[i]&mass.boot.2019$Donorspp==bs.diff2019$awaySpp[i]&mass.boot.2019$Treatment==bs.diff2019$Treatment[i])]]
  diff<-home-away
  bs.diff2019$diff.mean[i]<-mean(diff)
  bs.diff2019$diff.95lo[i]<-quantile(diff,probs=0.025)
  bs.diff2019$diff.95hi[i]<-quantile(diff,probs=0.975)
}


#Write out
write.csv(bs.diff2018, "./data/bs.diff2018.220627.csv")
write.csv(bs.diff2019, "./data/bs.diff2019.220627.csv")

#Try to get it all the way to bootstrapped Is-----------------------------
#Calculate using the bsVectors output lists
#2018
b.Is.2018<-data.frame(Pair=rep(c("ARTR-HECO","ARTR-POSE","ARTR-PSSP","HECO-POSE","HECO-PSSP","POSE-PSSP"),3),
                      Treatment=rep(c("Control","Feedback","Exclusion"),each=6))
b.Is.2018$Is<-rep(0,18)
b.Is.2018$Is.95lo<-rep(0,18)
b.Is.2018$Is.95hi<-rep(0,18)

for (i in 1:length(b.Is.2018[,1])){
  sp1<-substr(b.Is.2018$Pair[i],1,4)
  sp2<-substr(b.Is.2018$Pair[i],6,9)
  n11<-paste(sp1,sp1,b.Is.2018$Treatment[i],sep=".")
  n12<-paste(sp1,sp2,b.Is.2018$Treatment[i],sep=".")
  n22<-paste(sp2,sp2,b.Is.2018$Treatment[i],sep=".")
  n21<-paste(sp2,sp1,b.Is.2018$Treatment[i],sep=".")
  Is.vec<-log(bsVectors2018[[n11]]/bsVectors2018[[n12]])+log(bsVectors2018[[n22]]/bsVectors2018[[n21]])
  b.Is.2018$Is[i]<-mean(Is.vec[is.finite(Is.vec)==TRUE])
  b.Is.2018$Is.95lo[i]<-quantile(Is.vec[is.finite(Is.vec)==TRUE],probs = 0.025)
  b.Is.2018$Is.95hi[i]<-quantile(Is.vec[is.finite(Is.vec)==TRUE],probs = 0.975)
}

#2019
b.Is.2019<-b.Is.2018
for (i in 1:length(b.Is.2019[,1])){
  sp1<-substr(b.Is.2019$Pair[i],1,4)
  sp2<-substr(b.Is.2019$Pair[i],6,9)
  n11<-paste(sp1,sp1,b.Is.2019$Treatment[i],sep=".")
  n12<-paste(sp1,sp2,b.Is.2019$Treatment[i],sep=".")
  n22<-paste(sp2,sp2,b.Is.2019$Treatment[i],sep=".")
  n21<-paste(sp2,sp1,b.Is.2019$Treatment[i],sep=".")
  Is.vec<-log(bsVectors2019[[n11]]/bsVectors2019[[n12]])+log(bsVectors2019[[n22]]/bsVectors2019[[n21]])
  b.Is.2019$Is[i]<-mean(Is.vec[is.finite(Is.vec)==TRUE])
  b.Is.2019$Is.95lo[i]<-quantile(Is.vec[is.finite(Is.vec)==TRUE],probs = 0.025)
  b.Is.2019$Is.95hi[i]<-quantile(Is.vec[is.finite(Is.vec)==TRUE],probs = 0.975)
}

#Visualize BS IS----------------------------- 
ISall<-rbind(b.Is.2018,b.Is.2019)
ISall$Year<-rep(c(2018,2019),each=18)
ISall$Treatment<-factor(ISall$Treatment, levels=c("Control", "Feedback","Exclusion"))

P.Is<-ggplot(data=ISall, aes(x=Pair, y=Is, color=Treatment))+
  facet_wrap(~Year,nrow=1)+
  geom_pointrange(aes(ymin=Is.95lo, ymax=Is.95hi),
                  position=position_dodge(width=0.8))+
  coord_flip()+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  theme_bw()+
  labs(x="Species pair", y="Projected pairwise interaction strength")
P.Is
