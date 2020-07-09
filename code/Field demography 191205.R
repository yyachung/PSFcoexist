#Code to try and combine all life stage transitions in field PSF experiment
#AC June 2020

setwd("C:/Users/yc68991/Box Sync/PSFcoexist/R/PSFcoexist")
library(ggpubr)
library(ggplot2)

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
mass.boot.2018<-bsDemog(dat=all2018,nb=5000,bs.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
mass.boot.2019<-bsDemog(dat=all2019,nb=5000,bs.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")

# #Non-bootstrap
# source("./code/func_noBSDemog.R")
# mass.2018<-nobsDemog(dat=all2018,levels.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
# mass.2019<-nobsDemog(dat=all2019,levels.var="Trans.Don.Treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")

#Data clean up and adding labels-----------
#Delete final row from 2018 data since it has no treatment label
mass.boot.2018<-mass.boot.2018[1:60,]
mass.2018<-mass.2018[1:60,]

#Add labels back to data
mass.boot.2018$Transplant<-substr(mass.boot.2018$treat.combo,1,4)
mass.boot.2018$Donorspp<-substr(mass.boot.2018$treat.combo,6,9)
mass.boot.2018$Treatment<-substr(mass.boot.2018$treat.combo,11,100)
mass.boot.2019$Transplant<-substr(mass.boot.2019$treat.combo,1,4)
mass.boot.2019$Donorspp<-substr(mass.boot.2019$treat.combo,6,9)
mass.boot.2019$Treatment<-substr(mass.boot.2019$treat.combo,11,100)

# mass.2018$Transplant<-substr(mass.boot.2018$treat.combo,1,4)
# mass.2018$Donorspp<-substr(mass.boot.2018$treat.combo,6,9)
# mass.2018$Treatment<-substr(mass.boot.2018$treat.combo,11,100)
# mass.2019$Transplant<-substr(mass.boot.2019$treat.combo,1,4)
# mass.2019$Donorspp<-substr(mass.boot.2019$treat.combo,6,9)
# mass.2019$Treatment<-substr(mass.boot.2019$treat.combo,11,100)


#Make Intra/Inter label
mass.boot.2018$Type[which(mass.boot.2018$Donorspp!=mass.boot.2018$Transplant)]<-"Inter"
mass.boot.2018$Type[which(mass.boot.2018$Donorspp==mass.boot.2018$Transplant)]<-"Intra"
mass.boot.2018$Type <- factor(mass.boot.2018$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
mass.boot.2018$Treatment <- factor(mass.boot.2018$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization
mass.boot.2019$Type[which(mass.boot.2019$Donorspp!=mass.boot.2019$Transplant)]<-"Inter"
mass.boot.2019$Type[which(mass.boot.2019$Donorspp==mass.boot.2019$Transplant)]<-"Intra"
mass.boot.2019$Type <- factor(mass.boot.2019$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
mass.boot.2019$Treatment <- factor(mass.boot.2019$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization

# mass.2018$Type[which(mass.boot.2018$Donorspp!=mass.boot.2018$Transplant)]<-"Inter"
# mass.2018$Type[which(mass.boot.2018$Donorspp==mass.boot.2018$Transplant)]<-"Intra"
# mass.2018$Type <- factor(mass.boot.2018$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
# mass.2018$Treatment <- factor(mass.boot.2018$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization
# mass.2019$Type[which(mass.boot.2019$Donorspp!=mass.boot.2019$Transplant)]<-"Inter"
# mass.2019$Type[which(mass.boot.2019$Donorspp==mass.boot.2019$Transplant)]<-"Intra"
# mass.2019$Type <- factor(mass.boot.2019$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
# mass.2019$Treatment <- factor(mass.boot.2019$Treatment, levels=c("Control","Feedback","Exclusion"))#Flip the order of factors for easier visualization

#Visualize projections-----------------
#BOOTSTRAPPED
#2018
ggplot(data=mass.boot.2018, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
  theme_bw()
#2019
ggplot(data=mass.boot.2019, aes(x=Type, y=massContrib.mean, color=Treatment, shape=Donorspp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds")+
  theme_bw()
#Alternative visualization 
ggplot(data=mass.boot.2019, aes(x=Donorspp, y=massContrib.mean, color=Transplant, shape=Transplant))+
  facet_grid(~Treatment)+
  geom_pointrange(aes(ymin=massContrib.95CIlo, ymax=massContrib.95CIhi)
                  , position = position_dodge(width = 0.5))+
  geom_line(aes(group=Transplant), position=position_dodge(width=0.5))+
  scale_color_manual(values=c("#E2D200", "#46ACC8", "#E58601","#00A08A"))+
  xlab("Microsite")+ylab("Estimated shoot biomass contribution from 100 seeds (g)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

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

#Permutational differences 2018 only------------------
np <- 1000

#Make dataframe to store output
bs.diff<-data.frame(matrix(0,nrow=np, ncol=48))
colnames(bs.diff)<-mass.boot.2018$treat.combo[which(mass.boot.2018$Type=="Inter")]

#Permutation and calculate differences
for (i in 1:np){
  #Randomize treatment labels in full dataset
  all2018$bs.treat<-sample(all2018$Trans.Don.Treat, replace=F)
  
  #Estimate bootstrapped contributions
  source("./code/func_BSDemog.R")
  bs.means<-bsDemog(dat=all2018,nb=5000,bs.var="bs.treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
  
  # #Delete final row from 2018 data since it has no treatment label
  # mass.boot.2018<-mass.boot.2018[1:60,]
  # 
  # #Calculate means for each treatment combo ***This doesn't work because the data is already aggregated
  # bs.means<-aggregate.data.frame(mass.boot.2018$massContrib.mean, 
  #                                by=list(mass.boot.2018$bs.treat), FUN="mean")
  # colnames(bs.means)<-c("treat.combo","mean")
  bs.means$transplant<-substring(bs.means$treat.combo,1,4)
  bs.means$treatment<-substring(bs.means$treat.combo,11,100)
  bs.means$donor<-substring(bs.means$treat.combo,6,9)
  bs.means$trans.treat<-paste(bs.means$transplant,bs.means$treatment,sep=".")
  #Calculate difference between means intra/inter comparisons
  for (j in 1:12){
    subdata<-subset(bs.means,trans.treat==unique(bs.means$trans.treat)[j]) #this should have 5 rows
    intra<-subdata$massContrib.mean[which(subdata$transplant==subdata$donor)]
    subdata$diff<-subdata$massContrib.mean-intra
    for (k in 1:length(subdata[,1])){
      if (subdata$diff[k]!=0) {
        bs.diff[i,which(colnames(bs.diff)==subdata$treat.combo[k])] <- subdata$diff[k]
      }
      else next
    }
  }
}
#Took almost 7 houra
#Write out
write.csv(bs.diff, "permute1000.csv")

#Read back in
bs.diff<-read.csv("./data/permute1000.csv",row.names = 1)

##Check distributions
hist(bs.diff[,5]) #looks good

#Get 95% CI
null.diff<-as.data.frame(t(apply( bs.diff , MARGIN=2 , FUN=quantile ,probs = c(0.025,0.975) , na.rm = TRUE )))
#all extremely similar again

#calculate observed diffs
mass.boot.2018$trans.treat<-paste(mass.boot.2018$Transplant,mass.boot.2018$Treatment,sep=".")
for (j in 1:length(unique(mass.boot.2018$trans.treat))){
  subdata<-subset(mass.boot.2018,trans.treat==unique(mass.boot.2018$trans.treat)[j]) #this should have 5 rows
  intra<-subdata$massContrib.mean[which(subdata$Transplant==subdata$Donorspp)]
  subdata$diff<-subdata$massContrib.mean-intra
  for (k in 1:length(subdata[,1])){
    if (subdata$diff[k]!=0) {
      null.diff$obs[which(rownames(null.diff)==subdata$treat.combo[k])] <- subdata$diff[k]
    }
    else next
  }
}

View(null.diff) #Only 1 comparison is outside of CI. All CI's are surprisingly similar...

#Try again with permutational differences, controls only, 2019 data [This code is incomplete and abandoned]----
np <- 10
in.dat<-all2019[all2019$Treatment=="Control",] 

#Make dataframe to store output
perm.diff<-data.frame(matrix(0,nrow=np*4*10, ncol=5))
colnames(perm.diff)<-c("permIndex","Transplant","Contrast1","Contrast2","Diff")
perm.diff$permIndex<-rep(1:np,each=40)
perm.diff$Transplant<-rep(rep(c("ARTR","HECO","POSE","PSSP"),each=10),np)


for (i in 1:np){
  #Randomize treatment labels in full dataset
  in.dat$bs.treat<-sample(in.dat$Trans.Don.Treat, replace=F)
  
  #Estimate bootstrapped contributions
  source("./code/func_BSDemog.R")
  bs.means<-bsDemog(dat=in.dat,nb=5000,bs.var="bs.treat", germ.var="Germ",surv.var = "Surv",growth.var = "Abv")
  bs.means$transplant<-substring(bs.means$treat.combo,1,4)
  bs.means$treatment<-substring(bs.means$treat.combo,11,100)
  bs.means$donor<-substring(bs.means$treat.combo,6,9)
  
  #Calculate contrasts
  subdata<-subset(bs.means,transplant=perm.diff$Transplant[i])
  
  
}