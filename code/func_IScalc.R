IScalc<-function(dat){
  #make empty dataframe
  Is<-data.frame(Pair=character(),
                 Rep=integer(), 
                 Treatment=character(), 
                 lnRatio=double(),
                 stringsAsFactors=FALSE) 
  for (i in 1:10){
    subdata<-subset(dat,Rep==i)
    for (j in 1:6){
      pair<-unique(subdata$Pair)[j]
      pairdata<-subset(subdata,Pair==pair)
      #Add lnRatios for Is; if any element is NA, will return NA
      feedback.ratio<-sum(pairdata$lnRatio[pairdata$Treatment=="Feedback"])
      exclusion.ratio<-sum(pairdata$lnRatio[pairdata$Treatment=="Exclusion"])
      control.ratio<-sum(pairdata$lnRatio[pairdata$Treatment=="Control"])
      Is<-rbind(Is,list(Pair=rep(as.character(pair),3),
                        Rep=rep(i,3),
                        Treatment=c("Feedback","Exclusion","Control"),
                        lnRatio=c(feedback.ratio,exclusion.ratio,control.ratio)),
                stringsAsFactors=FALSE)
    }
    
  }
  return(Is)
}