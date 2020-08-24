############################################################
#Big global model for germ, surv, biomass                  #
#200717 AC                                                 #
############################################################

setwd("C:/Users/yyach/Desktop/Box Sync/PSFcoexist/R/PSFcoexist")
library(emmeans)
library(car)
library(lme4)
library(ggpubr)

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

#Aboveground biomass global model-------------------
m.abv<-lm(log(Abv)~Treatment*DonorSpp*Transplant+factor(Year),data=all)
qqPlot(m.abv$resid)#not perfect but scedasticity is better compared to no transformation
Anova(m.abv)#Treatment p~0, DonorSpp p<0.0001, Transplant p~0, Year p~0, Donor:Transplant p=0.00375
plotdf.abv<-summary(emmeans(m.abv,~DonorSpp|Transplant|Treatment,type = "response"))
                    
#Germination global model-------------------
m.germ<-glm(Germ~Treatment*DonorSpp*Transplant+factor(Year),data=all,family="poisson")
plot(m.germ)#decent
Anova(m.germ)#Type II Chisq test
#Treatment, Donorspp, Transplant, Year, sig. Treatment:Transplant p=0.062
plotdf.germ<-summary(emmeans(m.germ,~DonorSpp|Transplant|Treatment,type = "response"))

#Survival global model-------------------
m.surv<-glm(Surv~Treatment*DonorSpp*Transplant+factor(Year),data=all,family="binomial")
plot(m.surv)#hard to interpret
summary(m.surv)#hard to interpret
Anova(m.surv)#Chi-sq should still apply. Transplant sig.
#warning message about fitted prob 0/1, understandable since surv probs generally high
plotdf.surv<-summary(emmeans(m.surv,~DonorSpp|Transplant|Treatment,type = "response"))

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
  theme_bw()+
  theme(legend.position = "none")
abv

surv<-ggplot(data=plotdf.surv, aes(x=Type, y=prob, color=Treatment, shape=DonorSpp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=prob-SE, ymax=prob+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("Soil environment")+ylab("Survival probability")+
  theme_bw()+
  theme(legend.position = "none")
surv

germ<-ggplot(data=plotdf.germ, aes(x=Type, y=rate, color=Treatment, shape=DonorSpp))+
  facet_wrap(~Transplant,nrow=1)+
  geom_pointrange(aes(ymin=rate-SE, ymax=rate+SE)
                  , position = position_dodge(width = 1))+  
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  scale_shape_discrete(name="Microsite")+
  xlab("Soil environment")+ylab("Germination rate (out of 10)")+
  theme_bw()+
  theme(legend.position = "none")
germ

ggarrange(germ, surv, abv, ncol = 1, nrow = 3, align = "v",labels = c("A", "B","C"),
          common.legend = TRUE, legend = "right")

