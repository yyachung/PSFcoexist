PSFcalc<- function(dat,plantID){
  #Write empty data frame for output
  PSF<-data.frame(DonorSpp=character(),
                  Rep=integer(), 
                  Treatment=character(), 
                  lnRatio=double(),
                  stringsAsFactors=FALSE) 
  #Run loop to calculate PSFs and populate new dataframe
  for (i in 1:10){
    subdata<-subset(dat,Rep==i)
    home<-subset(subdata,DonorSpp==plantID)
    away<-subset(subdata,DonorSpp!=plantID)
    for (j in 1:4){
      Donor<-unique(away$DonorSpp)[j]
      feedback.ratio<-log(home$Abv[home$Treatment=="Feedback"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Feedback"])
      exclusion.ratio<-log(home$Abv[home$Treatment=="Exclusion"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Exclusion"])
      control.ratio<-log(home$Abv[home$Treatment=="Control"]/away$Abv[away$DonorSpp==Donor&home$Treatment=="Control"])
      PSF<-rbind(PSF,list(DonorSpp=rep(as.character(Donor),3),
                          Rep=rep(i,3),
                          Treatment=c("Feedback","Exclusion","Control"),
                          lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                      stringsAsFactors=FALSE)
    }
  }
  return(PSF)
}
  


