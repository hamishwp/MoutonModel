# directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
directory<-paste0(getwd(),"/")
setwd(directory)

library(tidyverse)
library(magrittr)
library(dissPackage3)
library(pracma)
source(paste0(directory,'Rcode/SMC.R'))

GetSoaySheep <-function(directory="/home/patten/Documents/Coding/Oxford/MoutonModel/",oneSex=T){

  SHEEP <- read.csv(paste0(directory,"/SHEEP.csv"), header = 2)
  
  # How many years are included in the data:
  SHEEP$sheep.yr %>% unique %>% length
  
  minyr<-min(SHEEP$sheep.yr)
  # Add some handy columns:
  SHEEP$census.number <- SHEEP$sheep.yr-minyr+1
  SHEEP$size <- log(SHEEP$wtt1)
  SHEEP$prev.size <- log(SHEEP$wt)
  SHEEP$rec1.wt<-log(SHEEP$rec1.wt)
  SHEEP$rec2.wt<-log(SHEEP$rec2.wt)
  SHEEP$off.born<-SHEEP$rec
  SHEEP$reproduced<-as.integer(SHEEP$off.born>0)
  SHEEP$survived<-SHEEP$surv
  
  SHEEP$wt<-SHEEP$wtt1<-SHEEP$sheep.yr<-SHEEP$pop.size<-
    SHEEP$surv<-SHEEP$mean.off.mass<-SHEEP$rec<-NULL
    
  # Try and find out what proportion of the population of females gets detected:
  # We make the assumption that detectability of females is 1 (and have confirmed
  # with experts that this is reasonable). So in fact, we know the true population
  # size at each time step, and the state space model is less appropriate.
  # However, we would like to see what the methodology produces in this context.
  # THE CAPTURED NUMBER IS LESS THAN THE DETECTED NUMBER
  detectedLiveFemales <- SHEEP%>%filter(survived==1)%>%group_by(census.number)%>%
    summarise(detectedNum=length(size),.groups = 'drop_last')%>%pull(detectedNum)%>%unname()
  
  return(list(solveDF=SHEEP,
              detectedNum=detectedLiveFemales))
  
}

GetSoaySheep_binned <-function(SHEEP,shift=0.49,oneSex=T,nbks=10){
  
  # Get breakpoints (histogram bin intervals), count data, and growth rate per census:
  # breaks<-quantile(c(SHEEP$prev.size,SHEEP$size),probs = 0:(nbks-1)/(nbks-1),na.rm = T)%>%unname
  # SHEEP$breaks<-seq(from=min(c(SHEEP$solveDF$prev.size,SHEEP$solveDF$size),na.rm = T),
  #             to=max(c(SHEEP$solveDF$prev.size,SHEEP$solveDF$size),na.rm = T),length.out=nbks)
  SHEEP$breaks<-ggplot_build(ggplot()+geom_histogram(data=SHEEP$solveDF,mapping=aes(size),bins = nbks))$data[[1]]$x
  SHEEP$sizes <- SHEEP$breaks[-(nbks)] + shift*diff(SHEEP$breaks)
  # Check how this pans out
  previousState<-vectorToCounts(SHEEP$solveDF$size, SHEEP$breaks)
  # Merge all bins that are in the lowest 5th percentile
  minnie<-SHEEP$breaks[min(which(cumsum(previousState)/sum(previousState)>0.05))]
  others<-seq(from=minnie,to=max(SHEEP$breaks),length.out=(nbks-1))
  SHEEP$breaks<-c(min(SHEEP$breaks),others)
  SHEEP$sizes <- SHEEP$breaks[-(nbks)] + shift*diff(SHEEP$breaks)
  # breaks<-histss(SHEEP$size[!is.na(SHEEP$size)],nbks)$breaks
  SHEEP$breaks[c(1,nbks)]<-c(-Inf,Inf)
  # breaks <- c(0, 16, 20, 24, 28, 40)
  SHEEP$COUNTS <- getSizeDistns(SHEEP$solveDF, SHEEP$breaks)
  SHEEP$priorProbs<-rowSums(SHEEP$COUNTS)/sum(SHEEP$COUNTS)
  
  SHEEP$obsProbTime <- apply(SHEEP$COUNTS, 2, sum)/SHEEP$detectedNum
  
  # Get the pars of the initial size dist:
  # SHEEP$sheepStartSize <- SHEEP$detectedNum
  # SHEEP$sheepMuPar <- SHEEP$COUNTS[,1] %>% countsToProbs
  
  return(SHEEP)
  
}
