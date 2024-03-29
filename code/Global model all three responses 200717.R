############################################################
#Big global model for germ, surv, biomass                  #
#200717 AC                                                 #
############################################################

library(emmeans)
library(car)
library(lme4)
library(ggpubr)
library(cowplot)
library(lmerTest)
library(knitr)

#Load data and process ---------------------------
all2018<-read.csv("./data/Transplant data clean 190102.csv") #2018 data
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 data

#Clean 2018 data
all2018<-all2018[-which(all2018$Treatment==""),]
all2018$Treatment <- factor(all2018$Treatment)#drop blank factor level

#Add year var and combine
all2018$Year<-rep(2018,length(all2018[,1]))
all2019$Year<-rep(2019,length(all2019[,1]))
#Annoying: the two datasets have different columns
colnames(all2018)
colnames(all2019)
combcol<-colnames(all2018)[colnames(all2018)%in%colnames(all2019)]
sub2018<-all2018[,combcol]
sub2019<-all2019[,combcol]
#Bind datasets together
all<-rbind.data.frame(sub2018,sub2019)
#Remove outlier aboveground biomass
all$Abv[which(all$Abv<0.0002)]<-NA

#Aboveground biomass global model-------------------
m.abv<-lm(log(Abv)~Treatment*DonorSpp*Transplant+factor(Year),data=all)
qqPlot(m.abv$resid)#not perfect but scedasticity is better compared to no transformation
summary(m.abv)
Anova(m.abv)#Treatment p~0, DonorSpp p<0.0001, Transplant p~0, Year p~0, Donor:Transplant p=0.00375
plotdf.abv<-summary(emmeans(m.abv,~DonorSpp|Transplant|Treatment,type = "response"))

#write.csv(plotdf.abv,"./data/Growth coefficients.csv")   
#Try mixed model 
all$block<-paste(all$Year,all$Rep,sep=".")
mixed.abv<-lmer(log(Abv)~Treatment*DonorSpp*Transplant+factor(Year)+(1|block),data=all)
summary(mixed.abv)
Anova(mixed.abv)#results stayed the same

#Germination global model-------------------
#Calculate #failues
summary(all$Germ)#max is 11 
all$Germ.fail<-10-all$Germ
all$Germ.fail[which(all$Germ.fail<0)]<-0
#Fit binomial model
m.germ<-glm(cbind(Germ,Germ.fail)~Treatment*DonorSpp*Transplant+factor(Year),data=all,family="binomial")
plot(m.germ)#decent
Anova(m.germ)#Type II Chisq Analysis of Deviance test
#Treatment, Donorspp, Transplant, Year, sig. Treatment:Transplant p=0.017
summary(m.germ)
plotdf.germ<-summary(emmeans(m.germ,~DonorSpp|Transplant|Treatment,type = "response"))

#Try mixed model 
mixed.germ<-glmer(cbind(Germ,Germ.fail)~Treatment*DonorSpp*Transplant+factor(Year)+(1|block),data=all,family="binomial")
summary(mixed.germ)
Anova(mixed.germ)#results stayed the same

#Survival global model-------------------
m.surv<-glm(Surv~Treatment*DonorSpp*Transplant+factor(Year),data=all,family="binomial")
plot(m.surv)#hard to interpret
summary(m.surv)#hard to interpret
Anova(m.surv)#Chi-sq should still apply. Transplant sig.
#warning message about fitted prob 0/1, understandable since surv probs generally high
plotdf.surv<-summary(emmeans(m.surv,~DonorSpp|Transplant|Treatment,type = "response"))

#Try mixed model 
mixed.surv<-glmer(Surv~Treatment*DonorSpp*Transplant+factor(Year)+(1|block),data=all,family="binomial")
summary(mixed.surv)
Anova(mixed.surv)#results stayed the same

#Make big multi-panel of all results----------------
#Make Intra/Inter label
plotdf.abv$Type<-NA
plotdf.abv$Type[which(as.character(plotdf.abv$DonorSpp)!=as.character(plotdf.abv$Transplant))]<-"Inter"
plotdf.abv$Type[which(as.character(plotdf.abv$DonorSpp)==as.character(plotdf.abv$Transplant))]<-"Intra"
plotdf.abv$Type <- factor(plotdf.abv$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
plotdf.abv$Treatment <- factor(plotdf.abv$Treatment, levels=c("Control", "Feedback","Exclusion"))#Flip the order of factors for easier visualization
plotdf.surv$Type<-NA
plotdf.surv$Type[which(as.character(plotdf.surv$DonorSpp)!=as.character(plotdf.surv$Transplant))]<-"Inter"
plotdf.surv$Type[which(as.character(plotdf.surv$DonorSpp)==as.character(plotdf.surv$Transplant))]<-"Intra"
plotdf.surv$Type <- factor(plotdf.surv$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
plotdf.surv$Treatment <- factor(plotdf.surv$Treatment, levels=c("Control", "Feedback","Exclusion"))#Flip the order of factors for easier visualization
plotdf.germ$Type<-NA
plotdf.germ$Type[which(as.character(plotdf.germ$DonorSpp)!=as.character(plotdf.germ$Transplant))]<-"Inter"
plotdf.germ$Type[which(as.character(plotdf.germ$DonorSpp)==as.character(plotdf.germ$Transplant))]<-"Intra"
plotdf.germ$Type <- factor(plotdf.germ$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
plotdf.germ$Treatment <- factor(plotdf.germ$Treatment, levels=c("Control", "Feedback","Exclusion"))#Flip the order of factors for easier visualization

#Make subplots
abv<-ggplot(data=plotdf.abv, aes(x=Type, y=response, color=Treatment, shape=DonorSpp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("Soil environment")+ylab("Aboveground biomass (g)")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
abv

surv<-ggplot(data=plotdf.surv, aes(x=Type, y=prob, color=Treatment, shape=DonorSpp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=prob-SE, ymax=prob+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  ylim(0,1)+
  xlab("Soil environment")+ylab("Survival probability")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
surv

germ<-ggplot(data=plotdf.germ, aes(x=Type, y=prob, color=Treatment, shape=DonorSpp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=prob-SE, ymax=prob+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("Soil environment")+ylab("Germination probability")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
germ

ggarrange(germ, surv, abv, ncol = 1, nrow = 3, align = "v",labels = c("A", "B","C"),
          common.legend = TRUE, legend = "right")

#Alternative multi=panel plot with better spacing and annotations------------------
#Make Intra/Inter label
plotdf.abv$Type<-NA
plotdf.abv$Type[which(as.character(plotdf.abv$DonorSpp)!=as.character(plotdf.abv$Transplant))]<-"Inter"
plotdf.abv$Type[which(as.character(plotdf.abv$DonorSpp)==as.character(plotdf.abv$Transplant))]<-"Intra"
plotdf.abv$Type <- factor(plotdf.abv$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
plotdf.abv$Treatment <- factor(plotdf.abv$Treatment, levels=c("Control", "Feedback","Exclusion"))#Flip the order of factors for easier visualization
plotdf.surv$Type<-NA
plotdf.surv$Type[which(as.character(plotdf.surv$DonorSpp)!=as.character(plotdf.surv$Transplant))]<-"Inter"
plotdf.surv$Type[which(as.character(plotdf.surv$DonorSpp)==as.character(plotdf.surv$Transplant))]<-"Intra"
plotdf.surv$Type <- factor(plotdf.surv$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
plotdf.surv$Treatment <- factor(plotdf.surv$Treatment, levels=c("Control", "Feedback","Exclusion"))#Flip the order of factors for easier visualization
plotdf.germ$Type<-NA
plotdf.germ$Type[which(as.character(plotdf.germ$DonorSpp)!=as.character(plotdf.germ$Transplant))]<-"Inter"
plotdf.germ$Type[which(as.character(plotdf.germ$DonorSpp)==as.character(plotdf.germ$Transplant))]<-"Intra"
plotdf.germ$Type <- factor(plotdf.germ$Type, levels=c("Intra", "Inter"))#Flip the order of factors for easier visualization
plotdf.germ$Treatment <- factor(plotdf.germ$Treatment, levels=c("Control", "Feedback","Exclusion"))#Flip the order of factors for easier visualization

#Abv
plotdf.abv$plotX<-paste(plotdf.abv$Type,plotdf.abv$DonorSpp)
plotdf.abv$plotX <- factor(plotdf.abv$plotX, 
                           levels=c("Intra ARTR", "Intra HECO","Intra POSE","Intra PSSP",
                                    "Inter ARTR", "Inter BARE", "Inter HECO", "Inter POSE", "Inter PSSP"))#Flip the order of factors for easier visualization
abv<-ggplot(data=plotdf.abv, aes(x=plotX, y=response, color=Treatment, shape=DonorSpp))+
  facet_grid(~Transplant,scales="free_x")+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("\n Soil environment")+ylab("Aboveground biomass (g)")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
abvp<-abv+coord_cartesian(clip = "off") + # allows plotting anywhere on the canvas
  draw_label("Intra", x = 1, y = -0.012, size=10)+
  draw_label("Inter", x = 3.5, y = -0.012, size=10)

#Surv
plotdf.surv$plotX<-paste(plotdf.surv$Type,plotdf.surv$DonorSpp)
plotdf.surv$plotX <- factor(plotdf.surv$plotX, 
                           levels=c("Intra ARTR", "Intra HECO","Intra POSE","Intra PSSP",
                                    "Inter ARTR", "Inter BARE", "Inter HECO", "Inter POSE", "Inter PSSP"))#Flip the order of factors for easier visualization
surv<-ggplot(data=plotdf.surv, aes(x=plotX, y=prob, color=Treatment, shape=DonorSpp))+
  facet_grid(~Transplant,scales="free_x")+
  geom_pointrange(aes(ymin=prob-SE, ymax=prob+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("\n Soil environment")+ylab("Survival probability")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  ylim(0,1.05)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
survp<-surv+coord_cartesian(clip = "off") + # allows plotting anywhere on the canvas
  draw_label("Intra", x = 1, y = -0.08, size=10)+
  draw_label("Inter", x = 3.5, y = -0.08, size=10)

#Germ
plotdf.germ$plotX<-paste(plotdf.germ$Type,plotdf.germ$DonorSpp)
plotdf.germ$plotX <- factor(plotdf.germ$plotX, 
                            levels=c("Intra ARTR", "Intra HECO","Intra POSE","Intra PSSP",
                                     "Inter ARTR", "Inter BARE", "Inter HECO", "Inter POSE", "Inter PSSP"))#Flip the order of factors for easier visualization
germ<-ggplot(data=plotdf.germ, aes(x=plotX, y=prob, color=Treatment, shape=DonorSpp))+
  facet_grid(~Transplant,scales="free_x")+
  geom_pointrange(aes(ymin=prob-SE, ymax=prob+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("\n Soil environment")+ylab("Germination probability")+
  geom_vline(xintercept=1.5,color="black",linetype="dashed")+
  theme_bw()+
  ylim(0,0.4)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
germp<-germ+coord_cartesian(clip = "off") + # allows plotting anywhere on the canvas
  draw_label("Intra", x = 1, y = -0.03, size=10)+
  draw_label("Inter", x = 3.5, y = -0.03, size=10)

#Plot together
ggarrange(germp, survp, abvp, ncol = 1, nrow = 3, align = "v",labels = c("A", "B","C"),
          common.legend = TRUE, legend = "right")

#Correlating responses across vital rates----------------
ggplot(data=all, aes(x=Germ/10, y=Abv, color=Treatment, shape=DonorSpp))+
  facet_grid(~Transplant)+
  geom_point()+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_shape_discrete(name="Soil environment")+
  xlab("Germination probability")+ylab("Aboveground biomass")

