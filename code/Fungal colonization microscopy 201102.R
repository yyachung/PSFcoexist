####################################################
# Fungal colonization analyses from microscopy data#
# AC 201204                                        #
####################################################

setwd("C:/Users/yyach/Desktop/Box Sync/PSFcoexist/R/PSFcoexist")
library(emmeans)
library(car)
library(lme4)
library(ggpubr)

#Load data ------------------
#microdat<-read.csv("./data/Fungal Colonization Data 201102.csv")
microdat<-read.csv("./data/Fungal Colonization Data 201204.csv")
all2019<-read.csv("./data/Transplant data year 2 clean 190819.csv") #2019 transplant data to match IDs

#Match metadata ------------------
#First figure out if microscopy data has duplicate IDs
table(microdat$ID)#There are 2x 95's, 2x 406's, 2x394, 2x567 , 2x449
#Look at the duplicates
microdat[microdat$ID%in%c(95,406,394,567,449),]#Look pretty different, no notes
#Exclude these data for now
microdat2<-microdat[-which(microdat$ID%in%c(95,406,394,567,449)),]

#Merge data
total <- merge(microdat2,all2019[,2:5],by.x="ID",by.y="fieldID")

#Preliminary data visualization-------------
#Check replication and sample sizes
table(total[,11:13])#not super even coverage. More POSE and PSSP samples, ARTR and HECO only replicated 1-3 times
table(total[,"Transplant"])#HECO and ARTR about 1/2 of POSE and PSSP
table(total[,"DonorSpp"])#Decently even
table(total[,"Treatment"])#Fewer Feedback but pretty even

#Boxplots
boxplot((total$Blue+total$DSE)~total$Treatment)#okay the pattern makes sense, phew
boxplot((total$Blue+total$DSE)~total$Transplant)#pretty even. POSE is low but has high var
boxplot((total$Blue+total$DSE)~total$DonorSpp)#pretty even, HECO is slightly higher

#Calculate colonization rates----------------
total$Blue.rate<-total$Blue/total$Total
total$DSE.rate<-total$DSE/total$Total
total$BSF.rate<-total$BSF/total$Total
total$Extra.rate<-total$Extra.radical/total$Total
total$Ves.rate<-total$Vesicles/total$Total

#Do ANOVAs for each category of structures-------------
#Blue Aseptate
m.blue<-lm(Blue.rate~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.blue$residuals)#looks great
Anova(m.blue)#Treatment, transplant, donor main effects; sig Treat:Trans
pairs(emmeans(m.blue, ~Treatment))
#Control>Exclusion p<0.0001, Control>Feedback p<0.0001, Feedback about the same as exclusion
#So the fungi don't like living in bags?
pairs(emmeans(m.blue, ~Transplant))#ARTR=HECO=PSSP>POSE
ggerrorplot(data = total, x="Transplant",y="Blue.rate",
            facet="Treatment")
#HECO is STILL the outier here, with high colonization of blue in the Exclusion but low in Feedback
table(total$Treatment[total$Transplant=="HECO"])#Hmm sample sizes look fine
pairs(emmeans(m.blue, ~DonorSpp))#No sig direct pairwise differences
ggerrorplot(data = total, x="DonorSpp",y="Blue.rate")#Looks like HECO has high and BARE has low donor effects (pairwise P=0.07)
ggerrorplot(data = total, x="Transplant",y="Blue.rate")#POSE low everyone else high
#Plot 3-way interaction just to see
ggerrorplot(data = total, x="DonorSpp",y="Blue.rate",
            facet.by=c("Transplant","Treatment"))


#DSE
m.dse<-lm(DSE.rate~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.dse$residuals)#looks bad
hist(m.dse$residuals)#not too bad here, but long tail?
#log-transform
m.dse2<-lm(log(DSE.rate+1)~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.dse2$residuals)#not much better
m.dse3<-lm(log(DSE.rate+0.01)~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.dse3$residuals)#better
#ANOVA
Anova(m.dse3)#Treatment, Transplant sig , Donor p=0.052
pairs(emmeans(m.dse3, ~Treatment))#All significantly different from each other, pattern as expected
ggerrorplot(data = total, x="Transplant",y="DSE.rate",
            facet="Treatment")#Treatment effect most obvious
ggerrorplot(data = total, x="Transplant",y="DSE.rate")#HECO high
pairs(emmeans(m.dse3, ~Transplant))#No sig direct pairwise differences 
ggerrorplot(data = total, x="DonorSpp",y="DSE.rate")#HECO high, ARTR/BARE low
pairs(emmeans(m.dse3, ~DonorSpp))#No sig direct pairwise differences (ARTR-HECO 0.051)
#Plot 3-way interaction just to see
ggerrorplot(data = total, x="DonorSpp",y="DSE.rate",
            facet.by=c("Transplant","Treatment"))

#BSF
m.bsf<-lm(BSF.rate~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.bsf$residuals)#looks bad
#log transform
m.bsf1<-lm(log(BSF.rate+0.01)~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.bsf1$residuals)#good
Anova(m.bsf1)#treatment sig
ggerrorplot(data = total, x="Transplant",y="BSF.rate",
            facet="Treatment")#overall really low colonization, Exclusion and feedback equally low

#Extra-radical fungi
m.extra<-lm(Extra.rate~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.extra$residuals)#not great
#log transform
m.extra1<-lm(log(Extra.rate+0.01)~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.extra1$residuals)#good
Anova(m.extra1)#Treatment, Donor, and Treatment:Transplant sig
pairs(emmeans(m.extra1, ~Treatment))#Weirdly control and exclusion both sig higher than feedback???
pairs(emmeans(m.extra1, ~Treatment|Transplant))#Still driven by weird patterns in HECO
pairs(emmeans(m.extra1, ~DonorSpp))#HECO=PSSP>POSE=BARE=ARTR
ggerrorplot(data = total, x="Transplant",y="Extra.rate",
            facet="Treatment")
#Still weird outlier for HECO exclusion where it's the highest??

#Vesicles
m.ves<-lm(Ves.rate~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.ves$residuals)#not great
#log transform
m.ves1<-lm(log(Ves.rate+0.01)~Treatment*Transplant*DonorSpp-Treatment:Transplant:DonorSpp,data = total)
qqPlot(m.ves1$residuals)#better
Anova(m.ves1)#Sig treatment, transplant, donor
ggerrorplot(data = total, x="Transplant",y="Ves.rate",
            facet="Treatment")
pairs(emmeans(m.ves1, ~Transplant))#ARTR>POSE sig, rank: ARTR>HECO>PSSP>POSE
pairs(emmeans(m.ves1, ~DonorSpp))#rank: HECO>PSSP=POSE=ARTR=BARE; HECO>ARTR and HECO>BARE sig
ggerrorplot(data = total, x="Transplant",y="Ves.rate") #ARTR high, POSE low
ggerrorplot(data = total, x="DonorSpp",y="Ves.rate")#HECO high, everything else low

#MANOVA---------------------------------
mandata<-na.omit(total)
manova<-manova(cbind(mandata$Blue.rate,mandata$DSE.rate,mandata$BSF.rate,mandata$Extra.rate,mandata$Ves.rate)~mandata$Treatment*mandata$Transplant*mandata$DonorSpp-mandata$Treatment:mandata$Transplant:mandata$DonorSpp)
summary(manova)#treatment, transplant, donor sig
summary.aov(manova) #broken down for each response
emmeans(manova, pairwise ~ Treatment)#Control significantly higher than exclusion or feedback, feedback is insignificantly higher than exclusion

#Investigate correlation among fungal features--------------
cor(na.omit(total[,14:18]))
#Strongest correlation is between Blue and Vesicles
library(ggcorrplot)
cor_pmat(na.omit(total[,14:18]))#all very significant
#Visualize
library(GGally)
ggpairs(na.omit(total[,14:18]))
#Sure there's a correlation but the data points are actually not super convincing

#Prelim viz of fungi vs. biomass----------------
#Add Abv mass variable
mass <- merge(total,all2019[,c(2,14)],by.x="ID",by.y="fieldID")
#Get rid of the one POSE outlier
mass <- mass[-which(mass$Abv<0.001),]

#Prelim viz
ggplot(mass, aes(x=Blue.rate, y=log(Abv),col=Treatment, shape=Transplant)) + geom_point()
ggplot(mass, aes(x=log(DSE.rate+0.1), y=log(Abv),col=Treatment, shape=Transplant)) + geom_point()
ggplot(mass, aes(x=log(Ves.rate+0.1), y=log(Abv),col=Treatment, shape=Transplant)) + geom_point()
#Hard to tell if there's anything going on

#Fungi vs. biomass regression------------------
#Use fungal colonization rates instead of treatments
#Blue aseptate
massBlue<-lm(Abv~Blue.rate*Transplant*DonorSpp, data=mass)
qqPlot(massBlue$residuals)#Terrible
#log transform
massBlue2<-lm(log(Abv)~Blue.rate*Transplant*DonorSpp, data=mass)
qqPlot(massBlue2$residuals)#Better but still not great
summary(massBlue2)
#Blue has negative slope (p=0.017) and a bunch of sig interactive effects
#Mostly to do with 1)low biomass in HECO, and 2) positive association in ARTR with BARE donor
Anova(massBlue2)#Blue and Transplant sig
#Visualize
ggplot(mass, aes(x=Blue.rate, y=log(Abv),col=Treatment)) + 
  geom_point()+
  facet_grid(cols=vars(Transplant),rows=vars(DonorSpp))

#DSE
massDSE<-lm(log(Abv)~DSE.rate*Transplant*DonorSpp, data=mass)
qqPlot(massDSE$residuals)#Similarly not great
summary(massDSE)#Only sig Transplant effects on weight
Anova(massDSE)#DSE and Transplant sig
#Visualize
ggplot(mass, aes(x=DSE.rate, y=log(Abv),col=Treatment, shape=DonorSpp)) + 
  geom_point()+
  facet_wrap(~Transplant, nrow=2)

#BSF
massBSF<-lm(log(Abv)~BSF.rate*Transplant*DonorSpp, data=mass)
qqPlot(massBSF$residuals)#Similarly not great
summary(massBSF)#BSF rate is p=0.052 (negative), sig interaction with PSSP donor (p=0.02) and HECO/HECO (p=0.47)
Anova(massBSF)#Transplant sig
#Visualize
ggplot(mass, aes(x=BSF.rate, y=log(Abv))) + 
  geom_point()+
  facet_grid(cols=vars(Transplant),rows=vars(DonorSpp))

#Extraradical
massExtra<-lm(log(Abv)~Extra.rate*Transplant*DonorSpp, data=mass)
qqPlot(massExtra$residuals)#not bad
summary(massExtra)#no sig effects
Anova(massExtra)#Transplant sig
#Visualize
ggplot(mass, aes(x=Extra.rate, y=log(Abv))) + 
  geom_point()+
  facet_grid(cols=vars(Transplant),rows=vars(DonorSpp))

#Vesicles
massVes<-lm(log(Abv)~Ves.rate*Transplant*DonorSpp, data=mass)
qqPlot(massVes$residuals)#not fantastic
summary(massVes)#higher in POSE and PSSP, significant neg slope in HECO/HECO
Anova(massVes)#Ves.rate, transplant, ves.rate:Donor sig
#Visualize
ggplot(mass, aes(x=Ves.rate, y=log(Abv))) + 
  geom_point()+
  facet_grid(cols=vars(Transplant),rows=vars(DonorSpp))

#Try PCA---------------------------
noNA<- na.omit(mass)
PCA<-prcomp(noNA[,14:18], center=T, scale=T)
summary(PCA)#the first 2 PC's only explain 62% of variance
PCA #First PC is positively associated with all, second PC is mostly negatively assoviated with extra-radical
biplot(PCA)#Giant clump that all goes in the same direction
#This may suffer from horseshoe effect
#But can't do NMDS with empty rows

#Pretty plot
library(ggfortify)
autoplot(PCA)
autoplot(PCA, data = noNA, colour = 'Transplant',frame=T)#Everything is mixed together
autoplot(PCA, data = noNA, colour = 'DonorSpp',frame=T)#Everything is mixed together
autoplot(PCA, data = noNA, colour = 'Treatment',frame=T)#Everything is mixed together

