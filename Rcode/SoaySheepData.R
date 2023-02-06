
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
propCOV<-diag(unlist((do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb,CI=T))))$sd))
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

