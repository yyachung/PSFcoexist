IScalc<-function(dat){
  #make empty dataframe
  Is<-data.frame(Pair=character(),
                 Rep=integer(), 
                 Treatment=character(), 
                 Is.m=double(),
                 Is.g=double(),
                 stringsAsFactors=FALSE) 
  for (i in 1:10){
    subdata<-subset(dat,Rep==i)
    for (j in 1:6){
      pair<-unique(subdata$Pair)[j]
      pairdata<-subset(subdata,Pair==pair)
      #Add lnRatios for Is; if any element is NA, will return NA
      feedback.ratio.m<-sum(pairdata$lnRatio.mass[pairdata$Treatment=="Feedback"])
      exclusion.ratio.m<-sum(pairdata$lnRatio.mass[pairdata$Treatment=="Exclusion"])
      control.ratio.m<-sum(pairdata$lnRatio.mass[pairdata$Treatment=="Control"])
      feedback.ratio.g<-sum(pairdata$lnRatio.germ[pairdata$Treatment=="Feedback"])
      exclusion.ratio.g<-sum(pairdata$lnRatio.germ[pairdata$Treatment=="Exclusion"])
      control.ratio.g<-sum(pairdata$lnRatio.germ[pairdata$Treatment=="Control"])
      Is<-rbind(Is,list(Pair=rep(as.character(pair),3),
                        Rep=rep(i,3),
                        Treatment=c("Feedback","Exclusion","Control"),
                        Is.m=c(feedback.ratio.m,exclusion.ratio.m,control.ratio.m),
                        Is.g=c(feedback.ratio.g,exclusion.ratio.g,control.ratio.g)),
                stringsAsFactors=FALSE)
    }
    
  }
  return(Is)
}