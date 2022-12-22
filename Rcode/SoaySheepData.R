
GetSoaySheep <-function(directory="/home/patten/Documents/Coding/Oxford/MoutonModel/",oneSex=T){

  SHEEP <- read.csv(paste0(directory,"SHEEP.csv"), header = 2)
  
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
  
  L<-min(c(SHEEP$prev.size, SHEEP$size),na.rm = T)
  U<-max(c(SHEEP$prev.size, SHEEP$size),na.rm = T)
  
  return(list(solveDF=SHEEP,
              detectedNum=detectedLiveFemales,
              L=L,U=U))
  
}

calcBreaks<-function(SHEEP,nbks,regbinspace=T){
  # Get breakpoints (histogram bin intervals), count data, and growth rate per census:
  if(regbinspace){
    breaks<-ggplot_build(ggplot()+geom_histogram(data=SHEEP$solveDF,mapping=aes(size),bins = nbks))$data[[1]]$x
    # Merge all bins that are in the lowest 5th percentile
    minnie<-quantile(SHEEP$solveDF$size,1/nbks,na.rm = T)%>%unname()
    breaks<-seq(from=minnie,to=max(breaks),length.out=nbks)
    # others<-seq(from=minnie,to=max(breaks),length.out=(nbks-1))
    # breaks<-c(min(breaks),others)
    if(breaks[nbks]>SHEEP$U) breaks[nbks]<-max(c(SHEEP$solveDF$size,SHEEP$solveDF$prev.size),na.rm = T)-1e-5
    if(breaks[1]>SHEEP$L) breaks[1]<-min(c(SHEEP$solveDF$size,SHEEP$solveDF$prev.size),na.rm = T)+1e-5
  } else {
    breaks<-quantile(c(SHEEP$solveDF$prev.size,SHEEP$solveDF$size),probs = 0:(nbks-1)/(nbks-1),na.rm = T)%>%unname
  }
  return(breaks)
}

GetSoaySheep_binned <-function(SHEEP,shift=0.5,oneSex=T,nbks=10,regbinspace=F){
  
  SHEEP$breaks<-calcBreaks(SHEEP,nbks,regbinspace)
  
  SHEEP$sizes <- SHEEP$breaks[-(nbks)] + shift*diff(SHEEP$breaks)
  breaks<-SHEEP$breaks; breaks[c(1,nbks)]<-c(-Inf,Inf)
  SHEEP$COUNTS <- getSizeDistns(SHEEP$solveDF, SHEEP$breaks)
  SHEEP$cNames<-c(vapply(paste0("yr_",as.character(SHEEP$solveDF$sheep.yr),"__"), 
                        function(yname) paste0(yname,paste0("sz_",as.character(signif(SHEEP$sizes,4)))),
                        character(length(SHEEP$sizes))))
  
  SHEEP$priorProbs<-rowSums(SHEEP$COUNTS)/sum(SHEEP$COUNTS)
  SHEEP$obsProbTime <- apply(SHEEP$COUNTS, 2, sum)/SHEEP$detectedNum
  
  # What form? array? We have size bins, variable and year, so 3D matrix?
  
  formData<-function(SHEEP, censy){
    # Filter by census first
    solveDF<-SHEEP$solveDF[SHEEP$solveDF$census.number==censy,]
    D<-length(SHEEP$breaks)
    # Initialise output array - dimension are c(size bins, validation variables)
    output<-array(NA,dim = c(D-1,5))
    # Sheep that survived from last census, based on previous size bin
    output[,1]<-vectorToCounts(c(solveDF$prev.size[solveDF$survived==1]), SHEEP$breaks)
    # All sheep (including newborns)
    output[,2]<-vectorToCounts(c(solveDF$size), SHEEP$breaks)
    # Sheep that were also parents in this censuss
    output[,3]<-vectorToCounts(c(solveDF$size[solveDF$reproduced==1]), SHEEP$breaks)
    # Growth of survived individuals - to infer change from one size bin to another
    output[,4]<-vectorToCounts(c(solveDF$size[solveDF$survived==1]), 
                               SHEEP$breaks) #- output[,1]
    # Number of offspring per size bin
    output[,5]<-vectorToCounts(c(solveDF$rec1.wt,solveDF$rec2.wt), SHEEP$breaks)

    return(output)
  }
  
  cen<-unique(SHEEP$solveDF$census.number)
  
  SHEEP$SumStats<-sapply(2:length(cen),
                         FUN = function(i) formData(SHEEP,cen[i]),
                         simplify = F) %>% unlist() %>%
                         array(dim=c(nbks-1,5,length(cen)))

  return(SHEEP)
  
}

###################################################################################
####################### Load data and the initial values ##########################
###################################################################################
# Read in the Soay sheep data
lSHEEP<-GetSoaySheep(directory,oneSex=oneSex)
# Make the inital values using the piecemeal GLM approach
x0<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb)))%>%relist(skeleton = skeleton)
# Number of parameters
Np<-length(unlist(x0))
# Provide initial estimate of covariance matrix using the confidence intervals from the GLM:
propCOV<-diag(unlist((do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb,CI=T))))$sd))/Np
# propCOV<-diag(Np)*(2.38)^2/Np
lSHEEP$breaks<-calcBreaks(lSHEEP,nbks,regbinspace = regbinspace)
# for truncated normal distributions
if(normsampler=="sampleDTN") {
  IPMLTP$DTN<-data.frame(L=lSHEEP$L,U=lSHEEP$U)
}
# Calculate the shift factor that offsets the size class mid-point
if(!manshift) {shift<-CalcShift_Kernel(x0 = Sample2Physical(x0,IPMLTP),IPMLTP = IPMLTP,nbks = nbks,
                                       breaks = lSHEEP$breaks,
                                       halfpop = oneSex,L = lSHEEP$L,U = lSHEEP$U)
} else shift<-0.5
print(paste0("Grid shift = ",shift, " for ",nbks," number of breaks." ))
# Get the sheep counts and sizes from the actual data (required even if simulated data is used)
lSHEEP<-GetSoaySheep_binned(lSHEEP,shift=shift,oneSex=T,nbks=nbks,regbinspace=regbinspace)  

