################################################################
#Code to analyze germination data during the experiment
################################################################

setwd("C:/Users/A02268715/Box Sync/PSFcoexist/R/PSFcoexist")
library(visreg)
library(car)
library(lme4)
library(emmeans)
library(ggpubr)

#Load and subset data---------------------------------------------
all<-read.csv("./data/Transplant data 180712.csv")

#Subset germination data only
germ<-all[,c("newID","DonorSpp","Transplant","Treatment","Germ180418","Germ180514","Germ180531","Germ180709")]

#Data processing------------------------------------------------
#Fill in with zeros
germ[is.na(germ)] <- 0

#Sum all germination records
germ$Total<-rowSums(germ[,6:8])#did not pull germs after first census

#Visualize distribution
hist(germ$Total)#lots of zeros

#Check sample sizes
aggregate(Total~DonorSpp*Transplant*Treatment, germ, function(x) sum(x != 0))
#HECO has the most missing data, followed by HECO
#No completely missing cells though

#Check means
aggregate(Total~DonorSpp*Transplant*Treatment, germ, FUN="mean")
#Pretty uniformly low...

#Try germination analysis---------------------------------------
#Global model
m1<-glm(Total~DonorSpp*Transplant*Treatment,data=germ,family=poisson())
hist(resid(m1))#pretty decent
summary(m1)#difficult to interpret, some indications of borderline significant effects
Anova(m1)#DonorSpp, Transplant, Treatment are sig
#Look at pairwise comparisons of main effects
m1.emm.donor <- emmeans(m1, "DonorSpp")
pairs(m1.emm.donor)#generally less germination on BARE, and better germ on ARTR than POSE
plot(m1.emm.donor, comparisons = TRUE)
boxplot(germ$Total~germ$DonorSpp)
m1.emm.plant <- emmeans(m1, "Transplant")
pairs(m1.emm.plant)#all pairwise contrasts are significant
plot(m1.emm.plant, comparisons = TRUE)
boxplot(germ$Total~germ$Transplant)
m1.emm.treat <- emmeans(m1, "Treatment")
pairs(m1.emm.treat)#control less than feedback
plot(m1.emm.treat, comparisons = TRUE)
boxplot(germ$Total~germ$Treatment)

#Try species separate since they have such different germ rates
ARTR<-subset(germ,Transplant=="ARTR")
POSE<-subset(germ,Transplant=="POSE")
PSSP<-subset(germ,Transplant=="PSSP")
HECO<-subset(germ,Transplant=="HECO")
#ARTR
m.ARTR<-glm(Total~DonorSpp*Treatment,data=ARTR,family=poisson())
hist(resid(m.ARTR))#pretty decent
summary(m.ARTR)#difficult to interpret, some indications of borderline significant effects
Anova(m.ARTR)#Treatment p=0.03, donor p=0.06, interaction p=0.19
m.ARTR.treat <- emmeans(m.ARTR, "Treatment")
pairs(m.ARTR.treat)#control less than feedback
plot(m.ARTR.treat, comparisons = TRUE)
#visualize PSF...no obvious fungal-mediated feedbacks
p.ARTR <- ggboxplot(ARTR, x = "Treatment", y = "Total",
               color = "DonorSpp", 
               add = "jitter", shape = "DonorSpp")
#POSE
m.POSE<-glm(Total~DonorSpp*Treatment,data=POSE,family=poisson())
hist(resid(m.POSE))#pretty decent
summary(m.POSE)#difficult to interpret, some indications of borderline significant effects
Anova(m.POSE)#Treatment p=0.13, donor p<0.0001, interaction p=0.72
m.POSE.donor <- emmeans(m.POSE, "DonorSpp")
pairs(m.POSE.donor)#Low germination in bare
plot(m.POSE.donor, comparisons = TRUE)
#visualize PSF...no obvious fungal-mediated feedbacks
p.POSE <- ggboxplot(POSE, x = "Treatment", y = "Total",
                    color = "DonorSpp", 
                    add = "jitter", shape = "DonorSpp")
#PSSP
m.PSSP<-glm(Total~DonorSpp*Treatment,data=PSSP,family=poisson())
hist(resid(m.PSSP))#pretty decent
summary(m.PSSP)#difficult to interpret, some indications of borderline significant effects
Anova(m.PSSP)#Treatment p=0.18, donor p=0.003, interaction p=0.09
m.PSSP.donor <- emmeans(m.PSSP, "DonorSpp")
pairs(m.PSSP.donor)#Low germination in bare
plot(m.PSSP.donor, comparisons = TRUE)
#visualize PSF...no obvious fungal-mediated feedbacks
p.PSSP <- ggboxplot(PSSP, x = "Treatment", y = "Total",
                    color = "DonorSpp", 
                    add = "jitter", shape = "DonorSpp")
#HECO
m.HECO<-glm(Total~DonorSpp*Treatment,data=HECO,family=poisson())
hist(resid(m.HECO))#not the greatest
summary(m.HECO)#donor effects
Anova(m.HECO)#Treatment p=0.66, donor p=0.003, interaction p=0.91
m.HECO.donor <- emmeans(m.HECO, "DonorSpp")
pairs(m.HECO.donor)#Low germination in bare
plot(m.HECO.donor, comparisons = TRUE)
#visualize PSF...no obvious fungal-mediated feedbacks, very low germination overall
p.HECO <- ggboxplot(HECO, x = "Treatment", y = "Total",
                    color = "DonorSpp", 
                    add = "jitter", shape = "DonorSpp")
