PSFcalc<- function(dat,plantID, germ.var, mass.var){
  #Write empty data frame for output
  PSF<-data.frame(DonorSpp=character(),
                  Rep=integer(), 
                  Treatment=character(), 
                  lnRatio.mass=double(),
                  lnRatio.germ=double(),
                  stringsAsFactors=FALSE) 
  #Run loop to calculate PSFs and populate new dataframe
  for (i in 1:10){
    subdata<-subset(dat,Rep==i)
    home<-subset(subdata,DonorSpp==plantID)
    away<-subset(subdata,DonorSpp!=plantID)
    for (j in 1:4){
      Donor<-unique(away$DonorSpp)[j]
      #mass
      feedback.ratio.m<-log(home[home$Treatment=="Feedback",mass.var]/away[away$DonorSpp==Donor&home$Treatment=="Feedback",mass.var])
      exclusion.ratio.m<-log(home[home$Treatment=="Exclusion",mass.var]/away[away$DonorSpp==Donor&home$Treatment=="Exclusion",mass.var])
      control.ratio.m<-log(home[home$Treatment=="Control",mass.var]/away[away$DonorSpp==Donor&home$Treatment=="Control",mass.var])
      #germ (add 1 to all numbers to correct for zeros)
      feedback.ratio.g<-log((home[home$Treatment=="Feedback",germ.var]+1)/(away[away$DonorSpp==Donor&home$Treatment=="Feedback",germ.var]+1))
      exclusion.ratio.g<-log((home[home$Treatment=="Exclusion",germ.var]+1)/(away[away$DonorSpp==Donor&home$Treatment=="Exclusion",germ.var]+1))
      control.ratio.g<-log((home[home$Treatment=="Control",germ.var]+1)/(away[away$DonorSpp==Donor&home$Treatment=="Control",germ.var]+1))
      #No surv because ln ratios would not make sense
      PSF<-rbind(PSF,list(DonorSpp=rep(as.character(Donor),3),
                          Rep=rep(i,3),
                          Treatment=c("Feedback","Exclusion","Control"),
                          lnRatio.mass=c(feedback.ratio.m,exclusion.ratio.m,control.ratio.m),
                          lnRatio.germ=c(feedback.ratio.g,exclusion.ratio.g,control.ratio.g)),
                      stringsAsFactors=FALSE)
    }
  }
  return(PSF)
}
  


