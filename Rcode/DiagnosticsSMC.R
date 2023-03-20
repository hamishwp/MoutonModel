
# Time the different components of the model-SMC to see where the code can be improved
TimeSMC<-function(x0,IPMLTP,timer=100000){
  
  proposed<-Sample2Physical(x0,IPMLTP)
  
  survFunc <- IPMLTP$survFunc
  survPars <- proposed$survPars
  growthSamp <- IPMLTP$growthSamp
  growthPars <- proposed$growthPars
  reprFunc <- IPMLTP$reprFunc
  reprPars <- proposed$reprPars
  offNumSamp <- IPMLTP$offNumSamp
  offNumPars <- proposed$offNumPars
  offSizeSamp <- IPMLTP$offSizeSamp
  breaks <- IPMLTP$breaks
  offSizePars <- proposed$offSizePars
  Schild<-proposed$Schild
  sizes<-IPMLTP$sizes
  oneSex <- IPMLTP$oneSex
  
  previousState<-IPMLTP$Y[,3]
  
  # Initialise output data.frame
  outer<-data.frame()
  # Get the size distribution of individuals that survive:
  D<-length(previousState); reppie<-rep(0,D)
  # Generate the survival probability [0,1] for each 'size' bin
  # Using previousState (X - total #sheep per bin), predict number of OBSERVED sheep per bin
  # Generate a sample population at the bin midpoint (at 'sizes' values) but repeating the values
  ptm <- proc.time()[3]
  for(i in 1:timer) newSizesI <- rep(rbinom(n=D, size=previousState, prob = survFunc(sizes, survPars)),x=sizes)  
  newS_fin<-(proc.time()[3] - ptm)
  # What would the growth of such a population be in one year?
  ptm <- proc.time()[3]
  for(i in 1:timer) newSizes <- growthSamp(newSizesI, growthPars)
  grow_fin<-(proc.time()[3] - ptm)
  # Bin these by breaks
  ptm <- proc.time()[3]
  for(i in 1:timer) newCounts<-vectorToCounts(newSizes, breaks)
  newC_fin<-(proc.time()[3] - ptm)
  # Get the probability that each class reproduces based on number of sheep (not observed but actual)
  # stop("Something is wrong here, what is reprFunc meant to do and why did it have z=previousState")
  ptm <- proc.time()[3]
  for(i in 1:timer) reprProbs <- reprFunc(sizes, reprPars)
  repr_fin<-(proc.time()[3] - ptm)
  # Parent counts
  # stop("what is reprProbs meant to do? weightedSelection?")
  ptm <- proc.time()[3]
  for(i in 1:timer) reprCounts <- rbinom(D, newCounts, reprProbs)
  repC_fin<-(proc.time()[3] - ptm)
  # Check if any births occurred
  if (sum(reprCounts)==0) {
    
    vallies<-array(c(newCounts,
                     reprCounts,
                     reppie),
                   dim = c(length(sizes),3),
                   dimnames = list(round(sizes,digits = 2),
                                   c("NoSurv","NoParents","NoOff")))
    # Calculate the distance metric (timed)
    ptm <- proc.time()[3]
    MultiMod<-function(i) multinomialObs(IPMLTP$SumStats[,i,rep(3,3)], 
                                         abind(vallies[,i],vallies[,i],vallies[,i],along = 2),
                                         0.9, logy=T)
    sw<-rowSums(sapply(1:3,MultiMod))
    dist_fin<-(proc.time()[3] - ptm)
    
    return(list(
      timings=data.frame(
        newS_fin=newS_fin,
        grow_fin=grow_fin,
        newC_fin=newC_fin,
        repr_fin=repr_fin,
        repC_fin=repC_fin,
        offC_fin=NA,
        offN_fin=NA,
        born_fin=NA,
        surC_fin=NA,
        repS_fin=NA,
        ofSS_fin=NA,
        dist_fin=dist_fin
      ),
      values = vallies,sw = sw))
  } else {
    # Modify the number of reproducing sheep to retain females only
    if(oneSex) {
      # Are all sheep included in the study female?
      ptm <- proc.time()[3]
      for(i in 1:timer) offCounts<-rbinom(D, reprCounts, 0.5)
      offC_fin<-(proc.time()[3] - ptm)
    } else offCounts<-reprCounts
    # Number of offspring born per PARENT size class
    ptm <- proc.time()[3]
    for(i in 1:timer) {
      offCountsS<-offCounts
      offCountsS[offCountsS!=0]<-offNumSamp(offCounts[offCounts!=0],offNumPars)
    }
    offCounts<-offCountsS
    offN_fin<-(proc.time()[3] - ptm)
    # Save the total number born
    ptm <- proc.time()[3]
    for(i in 1:timer) bornCount<-sum(offCounts)
    born_fin<-(proc.time()[3] - ptm)
    # How many survive to the next year?
    ptm <- proc.time()[3]
    for(i in 1:timer) {
      offCountsS<-offCounts
      offCountsS <- rbinom(D, offCounts, Schild)
    }
    offCounts<-offCountsS
    surC_fin<-(proc.time()[3] - ptm)
    # Convert to parent sizes
    ptm <- proc.time()[3]
    for(i in 1:timer) reprSizes <- rep(sizes,offCounts)
    repS_fin<-(proc.time()[3] - ptm)
    # Sample the size of the offspring
    ptm <- proc.time()[3]
    for(i in 1:timer) offSizes <- offSizeSamp(reprSizes,offSizePars) 
    ofSS_fin<-(proc.time()[3] - ptm)
  }
  # Return the ABCSIR object with all the count data, for this year only. This will be combined with rbind later
  vallies<-array(c(newCounts,
                   reprCounts,
                   vectorToCounts(offSizes, breaks)),
                 dim = c(length(sizes),3),
                 dimnames = list(round(sizes,digits = 2),
                                 c("NoSurv","NoParents","NoOff")))
  # Calculate the distance metric (timed)
  ptm <- proc.time()[3]
  MultiMod<-function(i) multinomialObs(IPMLTP$SumStats[,i,rep(3,3)], 
                                       abind(vallies[,i],vallies[,i],vallies[,i],along = 2),
                                       0.9, logy=T)
  sw<-rowSums(sapply(1:3,MultiMod))
  dist_fin<-(proc.time()[3] - ptm)
  
  return(list(
    timings=data.frame(
      newS_fin=newS_fin,
      grow_fin=grow_fin,
      newC_fin=newC_fin,
      repr_fin=repr_fin,
      repC_fin=repC_fin,
      offC_fin=offC_fin,
      offN_fin=offN_fin,
      born_fin=born_fin,
      surC_fin=surC_fin,
      repS_fin=repS_fin,
      ofSS_fin=ofSS_fin,
      dist_fin=dist_fin
    ),
    values = vallies, sw = sw))
  
}

tmp<-TimeSMC(x0,IPMLTP,timer=1)

# Test out the model-SMC to check for numerical underflow (ESS too low)
TestPF<-function(x0,IPMLTP, samplez=10){
  
  Sd<-IPMLTP$SumStats
  proposed<-Sample2Physical(x0,IPMLTP)
  sampleStatePar <- list(survFunc = IPMLTP$survFunc,
                         survPars = proposed$survPars,
                         growthSamp = IPMLTP$growthSamp, 
                         growthPars = proposed$growthPars,
                         reprFunc = IPMLTP$reprFunc,
                         reprPars = proposed$reprPars,
                         offNumSamp = IPMLTP$offNumSamp,
                         offNumPars = proposed$offNumPars,
                         offSizeSamp = IPMLTP$offSizeSamp,
                         breaks = IPMLTP$breaks,
                         offSizePars = proposed$offSizePars,
                         Schild=proposed$Schild,
                         oneSex = IPMLTP$oneSex,
                         sizes=IPMLTP$sizes)
  
  t <- dim(Sd)[3]
  
  outperf<-data.frame()
  for(i in 1:samplez){
    # Setup initial states
    prevStates <- do.call(IPMLTP$mu, IPMLTP$muPar)
    # Setup weight matrix and standardised weight matrix:
    output <- list(distance=0, sw=rep(1,IPMLTP$b),
                   shat=array(NA, dim=c(nrow(Sd),3,t)))
    # Update weights for first time step:
    wArgs <- list(Sd=Sd, pobs=IPMLTP$obsProbPar, NoParts=IPMLTP$b)
    
    performy<-data.frame()
    
    for (time in 1:t){
      wArgs$time<-time
      # Importance resample from the previous states 
      particleIndices <- sample(1L:wArgs$NoParts, wArgs$NoParts, replace = T, 
                                prob = output$sw)
      sampledStates <-  prevStates[, particleIndices]
      # IPM push forward
      wArgs$Sstar <- vectorisedSamplerIPM_ABCSIR(sampledStates, sampleStatePar)
      # prevStates is the total population, by size bin
      prevStates<-wArgs$Sstar[,1,]+wArgs$Sstar[,3,]
      # Convert to observed from latent space & calculate the objective function vector
      output2 <- IPMLTP$obsProb(output,wArgs)
      if(any(is.na(output2$sw))) {saveRDS(list(wArgs=wArgs,outputN=output2,outputO=output),"./tmpPF.Rdata"); stop("")} else output<-output2
      
      performy%<>%rbind(data.frame(time=time,distance=output$distance,ESS=1/(sum(output$sw^2))))
      
    }
    outperf%<>%rbind(cbind(performy,data.frame(attempt=i)))
  }
  # return(abs(output2$shat-IPMLTP$SumStats))
  return(list(outperf=outperf,shat=output$shat, Sstar=wArgs$Sstar))
}

Perfy<-TestPF(x0,IPMLTP,samplez = 1)
mean(Perfy$outperf$distance)

xxx<-x0
xxx[2]<--xxx[2]
Perfy2<-TestPF(xxx,IPMLTP,samplez = 3)
mean(Perfy2$outperf$distance)

# Test the multinomial observation model
obsModel<-"MultinomObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)

PerfMNO%>%group_by(time)%>%summarise(avESS=mean(ESS),avDist=mean(distance))

# Test the poisson-multinomial observation model
obsModel<-"multinomPoisObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfPMNO<-TestPF(x0,IPMLTP,samplez = 3)

PerfPMNO%>%group_by(time)%>%summarise(avESS=mean(ESS),avDist=mean(distance))

# Test the poisson observation model
obsModel<-"PoisObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfPO<-TestPF(x0,IPMLTP,samplez = 3)

PerfPO%>%group_by(time)%>%summarise(avESS=mean(ESS),avDist=mean(distance))

# Test the binomial observation model
obsModel<-"BinomObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfBO<-TestPF(x0,IPMLTP,samplez = 3)

PerfBO%>%group_by(time)%>%summarise(avESS=mean(ESS),avDist=mean(distance))

# Test the adaptive MAD-based distance model
obsModel<-"MADadaptdist"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMADad<-TestPF(x0,IPMLTP,samplez = 3)

PerfMADad%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))

# Test the mean absolute error distance model
obsModel<-"MAEdist"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMAE<-TestPF(x0,IPMLTP,samplez = 3)

PerfMAE%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
















xxx<-summaries$mean
xxl<-summaries$lower
xxu<-summaries$upper

obsModel<-"MultinomObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(x0-0.5,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(x0+0.5,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))


# Test the poisson-multinomial observation model
obsModel<-"multinomPoisObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(x0-0.5,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(x0+0.5,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))


# Test the poisson observation model
obsModel<-"PoisObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(x0-0.5,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(x0+0.5,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))


# Test the binomial observation model
obsModel<-"BinomObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(x0-0.5,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(x0+0.5,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))

# Test the adaptive MAD-based distance model
obsModel<-"MADadaptdist"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(xxl,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(xxu,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))
PerfMNO$shat[,,10]
PerfMNOm05$shat[,,10]
PerfMNOp05$shat[,,10]
PerfMNOxxx$shat[,,10]

# Test the adaptive MAD-based distance model
obsModel<-"MAEdistVar"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(xxl,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(xxu,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))
PerfMNO$shat[,,10]
PerfMNOm05$shat[,,10]
PerfMNOp05$shat[,,10]
PerfMNOxxx$shat[,,10]

# Test the mean absolute error distance model
obsModel<-"MAEdist"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(xxl,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
# PerfMNOp05<-TestPF(xxu,IPMLTP,samplez = 3)

PerfMNO$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOm05$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOp05$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOxxx$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
# PerfMNO$outperf%>%group_by(time)%>%summarise(avESS=mean(ESS),avDist=mean(distance))


obsModel<-"KSlike"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(xxl,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNO$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOm05$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOxxx$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))

c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))
PerfMNO$shat[,,10]
PerfMNOm05$shat[,,10]
PerfMNOp05$shat[,,10]
PerfMNOxxx$shat[,,10]


obsModel<-"Fudger"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(xxl,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(xxu,IPMLTP,samplez = 3)
PerfMNO$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOm05$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOp05$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))
PerfMNOxxx$outperf%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))


obsModel<-"multinomMAE"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PerfMNO<-TestPF(x0,IPMLTP,samplez = 3)
PerfMNOm05<-TestPF(x0-0.5,IPMLTP,samplez = 3)
PerfMNOxxx<-TestPF(xxx,IPMLTP,samplez = 3)
PerfMNOp05<-TestPF(x0+0.5,IPMLTP,samplez = 3)
c(sum(PerfMNO$outperf$distance),sum(PerfMNOm05$outperf$distance),sum(PerfMNOp05$outperf$distance),sum(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$distance),median(PerfMNOm05$outperf$distance),median(PerfMNOp05$outperf$distance),median(PerfMNOxxx$outperf$distance))
c(median(PerfMNO$outperf$ESS),median(PerfMNOm05$outperf$ESS),median(PerfMNOp05$outperf$ESS),median(PerfMNOxxx$outperf$ESS))
c(min(PerfMNO$outperf$ESS),min(PerfMNOm05$outperf$ESS),min(PerfMNOp05$outperf$ESS),min(PerfMNOxxx$outperf$ESS))

PerfMAE%>%filter(time==max(time))%>%summarise(avESS=mean(ESS),avDist=mean(distance))


















############### TEST OUT THE DIFFERENCES IN THE METRICS ###############

# Test the poisson-multinomial observation model
obsModel<-"multinomPoisObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
PMN<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])
# Test the multinomial observation model
obsModel<-"MultinomObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
MN<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])
# Test the poisson observation model
obsModel<-"PoisObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
Pp<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])
obsModel<-"BinomObs"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
Bn<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])
obsModel<-"MADadaptdist"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
MADa<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])
obsModel<-"MAEdist"
source(paste0(directory,'Rcode/ObsDistance.R'))
IPMLTP$obsProb<-obsfun
MAE<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])

metry<-data.frame(PoisMultnom=PMN$sw,Multnom=MN$sw,Poisson=Pp$sw,
                  Binom=Bn$sw,MADadapt=MADa$sw,MAE=MAE$sw)

# saveRDS(metry,"./Results/metric_function_comparisons.Rdata")
GGally::ggpairs(metry[!is.infinite(metry$Multnom) & !is.infinite(metry$Binom),])




# Run particleFilter with varying NoParts
# First prepare everything for sampling, including initial values
proposed<-Sample2Physical(x0,IPMLTP)
stateSpaceSampArgs <- list(survFunc = IPMLTP$survFunc, survPars = proposed$survPars,
                           growthSamp = IPMLTP$growthSamp, growthPars = proposed$growthPars,
                           reprFunc = IPMLTP$reprFunc, reprPars = proposed$reprPars, 
                           offNumSamp = IPMLTP$offNumSamp, offNumPars = proposed$offNumPars,
                           offSizeSamp = IPMLTP$offSizeSamp, breaks = IPMLTP$breaks,
                           offSizePars = proposed$offSizePars, Schild=proposed$Schild,
                           sizes=IPMLTP$sizes, oneSex = IPMLTP$oneSex)

numsam<-2000
# Now we calculate the standard deviation of the samples, when using specific numbers of particles - NoParts
tmp<-c(unlist(mclapply(1:numsam, function(i) {particleFilter(Sd=lSHEEP$SumStats, mu=mufun, muPar=muPar, obsProb = obsfun,
                                                           sampleState = vectorisedSamplerIPM_ABCSIR,
                                                           sampleStatePar = stateSpaceSampArgs,
                                                           obsProbPar = obsProbPar, 
                                                           fixedObsProb=fixedObsProb,
                                                           NoParts = 500)$distance},mc.cores = ncores)))
# We can see how many samples we need to achieve a stable standard deviation by plotting the cumulative SD (around 1800)
plot(vapply(seq_along(tmp), function(i) sd(tmp[1:i]), 1),xlab = "Number of Samples",ylab = "Cumulative S.D.",main = "How many samples required for S.D. ~ 600")

parties<-(1:30)*500
coster<-sapply(parties,function(pp) {
  muPar$n<-pp
  c(unlist(mclapply(1:numsam, function(i) {particleFilter(Sd=lSHEEP$SumStats, mu=mufun, muPar=muPar, obsProb = obsfun,
                                                          sampleState = vectorisedSamplerIPM_ABCSIR,
                                                          sampleStatePar = stateSpaceSampArgs,
                                                          obsProbPar = obsProbPar, 
                                                          fixedObsProb=fixedObsProb,
                                                          NoParts = pp)$distance},mc.cores=ncores)))
})
# Calculate the running (k=3) standard deviation of the distance by number of particles
runSD<-data.frame(parties=parties[3:length(parties)],runSD=vapply(3:length(parties), function(i) sd(log(-apply(coster,2,mean))[(i-3):i]), 1))
# To find the ideal number of particles, fit various linear models trained with 
# the data using sequentially less of the lower no. particle values
# and predict when the gradient of the linear model equals zero.
# Note that this method most likely over-estimates the required number of particles.
tmp<-data.frame(parties=parties[1:(nrow(runSD)-1)],grad=vapply(1:(nrow(runSD)-1),function(i) (lm(runSD~parties,runSD[i:nrow(runSD),]))$coefficients[2],1))
# Predict the required number of mSMC particles!
SMC_parts<-round(unname(predict(lm(parties~poly(grad,4),tmp),newdata = data.frame(grad=c(0)))))

p<-runSD%>%ggplot(aes(parties,runSD))+geom_point()+geom_vline(xintercept = SMC_parts)+
  xlab("Number of Particles in mSMC")+ylab("Running S.D (k=3)")+
  ggtitle("Convergence in No. Particles by S.D")+theme(plot.title = element_text(hjust = 0.5))
q<-data.frame(meany=apply(coster,2,mean),parties=parties)%>%ggplot(aes(parties,meany))+geom_point()+geom_vline(xintercept = SMC_parts)+
  xlab("Number of Particles in mSMC")+ylab("Average Distance")+
  ggtitle("Convergence in No. Particles w.r.t Distance")+theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(directory,"Plots/Hamish/Convergence/NoPart_mSMC.png"),
       gridExtra::grid.arrange(q,p,nrow=1),width = 10,height = 4); gridExtra::grid.arrange(q,p,nrow=1)

saveRDS(list(SMC_parts=SMC_parts,parties=parties,coster=coster,runSD=runSD),paste0(directory,"Results/SMCPARTS_calc.RData"))