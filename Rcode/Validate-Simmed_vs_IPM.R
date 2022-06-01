directory<-paste0(getwd(),"/")

list.of.packages <- c("xtable","magrittr","doParallel","Rfast","mc2d", "abind")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source(paste0(directory,'/Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'/Rcode/SimulateData.R'))
source(paste0(directory,'/Rcode/ModelSpecIPM.R'))
source(paste0(directory,'/Rcode/piecemealFunctions.R'))
source(paste0(directory,'/Rcode/SMC.R'))
library(dissPackage3)
library(xtable)
library(mcmcse)
library(tidyverse)
library(magrittr)

# ############################ GENERATE THE DATA #################################
# 
# # Get the observations:
# set.seed(102010)
# simmedData <- do.call(simulateIBM, simPars)
# Define the breaks for the 4 group specification:
breaks <- seq(1.5, 3.55, l=6)[-2]
shift=0.49
D<-length(breaks)
sizes <- breaks[-D] + shift*diff(breaks)
breaks[c(1, D)] <- c(-Inf, Inf)
oneSex<-T

########## GENERATE THE SIMULATED DATA ##########
startValues <- list(
  survPars = c(-9.65, 3.77),
  growthPars = c(1.41, 0.56, 0.08),
  # reprPars = c(-7.23, 2.6),
  reprPars = c(0, 0.2),
  offNumPars = 1,
  offSizePars = c(0.36, 0.71, 0.16),
  Schild = 0.873,
  obsProbPar = 0.999 # not too close to 50 since this will hurt the chain
)
vals<-unlist(startValues)
simPars <- list(n=100, t=5000,
                # set survival details:
                survFunc = linLogit, survPars = vals[1:2],
                # set growth details:
                growthSamp = sampleDTN,
                growthPars = c(vals[3:5], 1.5, 3.55),
                # set reproduction details:
                reprFunc = linLogit, reprPars = vals[6:7],
                # set offspring number and size distribution details:
                offNum=vals[8], offNumSamp=returnConstant, offSizeSamp = sampleDTN,
                offSizePars = c(vals[9:11], 1.5, 3.55),
                # Child survival probability:
                Schild=vals[12], obsProb=vals[13],
                # set other miscelaneous parameters:
                Start = 2.7, thresh=1000, OneGend = oneSex,
                popPrint = F, verbose=F)

# simmedData<-readRDS(paste0(directory,"RDSobjects/simmedData"))
set.seed(102010)
simmedData <- do.call(simulateIBM, simPars)

# Make a list of the parameters for the particle filter:
ungulatePars <- list(
  survFunc = linLogit, survPars = c(-9.65, 3.77),
  growthSamp = sampleDTN, growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
  reprFunc = linLogit, reprPars = c(-7.23, 2.6), 
  offSizeSamp = sampleDTN, offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
  offNumSamp = returnConstant, offNumPars = 1,
  Schild=qlogis(0.873), oneSex = TRUE, breaks = breaks, shift=qlogis(0.44))


popSizeAtStart <- subset(simmedData, simmedData$census.number==1)$size %>%
  na.omit %>% length

maxCN <- max(simmedData$census.number)

simSize <- 200
popSizeMatrix <- matrix(NA, nrow=simSize, ncol=length(obsPopSizes)) 

projectStateSpace(sampleStateIPM, ungulatePars,
                  multnomMu(popSizeAtStart, muPs), maxCN-1) %>% apply(2, sum)








repeats <- 20
simTime <- 13
sizeSet <- c(5, 10, 25, 50, 200)
results <- matrix(NA, nrow = length(sizeSet) + 1, ncol = simTime)

set.seed(101010)

simParsOrig <- list(n=100, t=5000,
                    # set survival details:
                    survFunc = linLogit, survPars = c(-9.65, 3.77),
                    # set growth details:
                    growthSamp = sampleDTN,
                    growthPars = c(1.41, 0.56, 0.08, 1.5, 3.55),
                    # set reproduction details:
                    reprFunc = linLogit, reprPars = c(-7.23, 2.6),
                    # set offspring number and size distribution details:
                    offNumPars=1, offNumSamp=returnConstant,
                    offSizeSamp = sampleDTN,
                    offSizePars = c(0.36, 0.71, 0.16, 1.5, 3.55),
                    # Child survival probability:
                    Schild=0.873, obsProb=0.999,
                    # set other miscelaneous parameters:
                    Start = 2.7, thresh=1000, OneGend = T,
                    popPrint = F, verbose=F)

simIBM <- do.call(simulateIBM, simParsOrig)

shift<-0.6
nbks<-6
breaks <- seq(1.5, 3.55, l=nbks)

sizeDists <- getSizeDistns(simIBM, breaks)
sizeDists <- sizeDists[,(ncol(sizeDists)-simTime+1):ncol(sizeDists)]
results[1,] <- apply(sizeDists, 2, sum)

simPars <- c(simParsOrig[3:14],list(oneSex=T,checks = FALSE,verbose = FALSE))
# %<>% c(list(oneSex = simPars$OneGend, breaks = breaks, shift=qlogis(0.49)))
obsProb<-simParsOrig$obsProb;simPars$obsProb<-NULL

for (j in 1:length(sizeSet)){
  resultsj <- matrix(NA, nrow=repeats, ncol=simTime)
  for (i in 1:repeats){
    breaks <- seq(1.5, 3.55, l=sizeSet[j]+1)
    
    simPars$breaks<-breaks
    simPars$sizes<-breaks[-(length(breaks))] + shift*diff(breaks)
    
    startSize <- results[1, 1]
    sizeDists <- getSizeDistns(simIBM, breaks)
    sizeDists <- sizeDists[,(ncol(sizeDists)-simTime+1):ncol(sizeDists)]
    startDistn <- sizeDists[,1]/sum(sizeDists[,1])
    resultsj[i,] <- projectStateSpace(sampleStateIPM_red, simPars,
                                      multnomMu(popSize = startSize, probs = startDistn,pobs = simParsOrig$obsProb), simTime - 1) %>% apply(2, sum)
  }
  
  results[j+1, ] <- apply(resultsj, 2, mean)
}

plot(results[1,], type='l', col=adjustcolor(2, alpha=0.5), lwd=2,
     ylab='Population size', lty=2,
     xlab="Census number")

for (i in 1:length(sizeSet)){
  lines(results[i+1,], col=adjustcolor(col=i+1, alpha=0.5), lwd=2, lty=2)
}

legend('topleft', fill = adjustcolor(col = 2:(length(sizeSet)+1), alpha = 0.5),
       legend = as.character(sizeSet), bty = 'n', horiz = T,
       title = "Number of Size Classes")










# IPM required parameters and class info
# breaks <- seq(1.5, 3.55, l=6)[-2]
breaks <- seq(1.5, 3.55, l=8)
shift=0.49
D<-length(breaks)
sizes <- breaks[-D] + shift*diff(breaks)
breaks[c(1, D)] <- c(-Inf, Inf)
oneSex<-T
# Simulate IBM population parameterisation
set.seed(101010)
simParsOrig <- list(n=200, t=50000,
                    # set survival details:
                    survFunc = linLogit, survPars = c(-9.65, 3.77),
                    # set growth details:
                    growthSamp = sampleDTN,
                    growthPars = c(1.41, 0.56, 0.08, 1.5, 3.55),
                    # set reproduction details:
                    reprFunc = linLogit, reprPars = c(-7.23, 2.6),
                    # set offspring number and size distribution details:
                    offNumPars=1, offNumSamp=returnConstant,
                    offSizeSamp = sampleDTN,
                    offSizePars = c(0.36, 0.71, 0.16, 1.5, 3.55),
                    # Child survival probability:
                    Schild=0.873, obsProb=0.999,
                    # set other miscelaneous parameters:
                    Start = 2.9, thresh=15000, OneGend = T,
                    popPrint = F, verbose=F)

# Simulate the population, starting from every individual weighing exp(2.7) kgs
# simIBM <- do.call(simulateIBM, simParsOrig)
# startSizes<-sample(simIBM$size[simIBM$census.number==max(simIBM$census.number) & !is.na(simIBM$size)],size = 200)
# simParsOrig$Start<-startSizes
# simIBM <- do.call(simulateIBM, simParsOrig)
# Get the size distributions for the IPM
# sizeDists <- getSizeDistns(simIBM, breaks)

survFunc <- linLogit
survPars <- c(-9.65, 3.77)
growthSamp <- sampleDTN
growthPars <- c(1.41, 0.56, 0.08, 1.5, 3.55)
reprFunc <- linLogit
reprPars <- c(-7.23, 2.6)
offNumPars<-1
offNumSamp<-returnConstant
offSizeSamp <- sampleDTN
offSizePars <- c(0.36, 0.71, 0.16, 1.5, 3.55)
Schild<-0.873
obsProb<-0.999
OneGend = T
L<-1.5
U<-3.55


popinc<-kernelOneVar(m = 500, growthFunc = doublyTruncatedNormal,
                     growthPars = growthPars, survFunc = survFunc,
                     survPars = survPars, repFunc = reprFunc,
                     repPars = reprPars, offNum = offNumPars,
                     offSizeFunc =  doublyTruncatedNormal,
                     offSizePars = offSizePars, L = L, U = U,
                     childSurv = Schild,shift=0.5,halfPop = 0.5) %>%
  eigen %>% `$`(values) %>% `[`(1) %>% Re


vals<-c(10:500)
nsims<-length(vals)

output<-mclapply(vals,FUN = function(i) {
  optimise(f=function(shift) {
    
    tmp<-kernelOneVar(m = i, growthFunc = doublyTruncatedNormal,
                 growthPars = growthPars, survFunc = survFunc,
                 survPars = survPars, repFunc = reprFunc,
                 repPars = reprPars, offNum = offNumPars,
                 offSizeFunc =  doublyTruncatedNormal,
                 offSizePars = offSizePars, L = L, U = U,
                 childSurv = Schild,shift=shift,halfPop = 0.5) %>%
      eigen %>% `$`(values) %>% `[`(1) %>% Re
    
    return(abs(tmp-popinc))
      
  },
           lower = 0.,upper = 1.0)
},mc.cores = 32)

output<-data.frame(nbreaks=vals,
                    shift=unlist(output)[seq.int(from=1,to=2L*nsims,by=2)],
                    AbsError=unlist(output)[seq.int(from=2,to=2L*nsims,by=2)])

p<-ggplot(output)+geom_point(aes(nbreaks,shift,colour=AbsError))
q<-ggplot(output)+geom_point(aes(nbreaks,AbsError,colour=shift))
gridExtra::grid.arrange(p,q,nrow=1)




simIBM <- do.call(simulateIBM, simParsOrig)


IPMShift_optim<-function(shift=0.49, ll=6, kn=50,cost=T){
  
  # breaks<-quantile(c(simIBM$size),probs = 0:(ll-1)/(ll-1),na.rm = T)%>%unname
  breaks<-seq(from=min(simIBM$size,na.rm = T),to=max(simIBM$size,na.rm = T),length.out=ll)
  D<-length(breaks)
  sizes <- breaks[-D] + shift*diff(breaks)
  breaks[c(1, D)] <- c(-Inf, Inf)
  
  sizeDists <- getSizeDistns(simIBM, breaks)
  
  DF<-simIBM[simIBM$census.number==1,]
  previousState<-sizeDists[,1]
  
  OutprevIPM<-prevIPM<-sizeDists
  OutprevIPM[]<-NA
  
  # isamps<-is.na(previous.census$size)
  # isizes<-!is.na(simIBM$size)
  
  for (k in 1:kn){
    
    prevIPM[]<-NA
    previousState<-sizeDists[,1]
    
    # # Randomly sample sizes for missing size sheep from the population
    # previous.census$size[isamps]<-sample(simIBM$size[isizes],sum(isamps),replace = T)
    # # Assume all the NA sheep were still alive but just not weighed
    # previous.census$survived[isamps]<-1
    
    ######## BUILD 3 TYPES - 1) IBM, 2) IPM 3) IPM based on IBM at previous iteration ########
    for (i in 1:max(simIBM$census.number,na.rm = T)){
      
      D<-length(previousState)
      
      ### IPM ### - updating the existing sheep by survival and growth
      newIPM <-  growthSamp(rep(rbinom(survFunc(sizes, survPars),
                                       n=D, size=previousState),x=sizes),growthPars)
      
      # Determine who reproduces (out of the survivors):
      #¼¼¼¼¼¼¼¼¼¼ The difference here is that the IBM retains information about the parent-child heritage ¼¼¼¼¼¼¼¼¼¼#
      ### IPM ###
      # reprobIPM <- reprFunc(newIPM, reprPars)
      reprobIPM <- reprFunc(rep(sizes,previousState), reprPars)
      reprIPM <- weightedSelection(newIPM, reprobIPM) # equivalent of parentsDF$size
      
      # Determine the number of offspring for each parent & Determine the sizes of the offspring:
      #¼¼¼¼¼¼¼¼¼¼ Difference is from using rep(reprIBM, NoffIBM) instead of NoffIBM directly in offIBM
      ### IPM ###
      offIPM<-offSizeSamp(reprIPM[!is.na(reprIPM)],offSizePars)
      offIPM<-offIPM[!is.na(offIPM)]
      
      # Determine the survival of the offspring:
      ### IPM ###
      offIPM <- weightedSelection(offIPM,rep(Schild,length(offIPM)))
      
      # Determine the sex of the offspring:
      #¼¼¼¼¼¼¼¼¼¼ Difference found between onesex and OneGend
      if(OneGend) {
        ### IPM ###
        offIPM <- weightedSelection(offIPM,rep(0.5,length(offIPM)))
      } 
      ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2222
      
      # Find out which children make it to census: # NOTE: IPMs are already done (above)!
      ### IPM ###
      # Update the binned sheep info of current and newly born sheep
      prevIPM[,i]<-previousState<-vectorToCounts(c(offIPM, newIPM), breaks)
      
    }
    
    if(k==1){
      OutprevIPM<-prevIPM
    } else {
      OutprevIPM<-OutprevIPM+prevIPM/k
    }
    
  }
  
  if(cost)  {
    colly<-colSums(OutprevIPM)
    return(mean((colly[2:length(colly)]/colly[1:(length(colly)-1)] - popinc)^2))
    # return(mean(abs(OutprevIPM[,5:ncol(OutprevIPM)]/OutprevIPM[,4:(ncol(OutprevIPM)-1)] - popinc)))
  }
  
  return(list(OutprevIPM=OutprevIPM))
  
}

vals<-c(3,4,5,6,7,8,9,10,15,20,25,30,40,50,80,100,150,200,250,300,400,500)
nsims<-length(vals)

outs<-list()
for (rr in 1:10){
  
  print(rr)
  output<-mclapply(vals,FUN = function(i) {
    ran<-runif(1,0.1,0.2)
    optimise(f=function(shift) IPMShift_optim(shift, ll=i, kn=6700,cost=T),
                                                   lower = ran,upper = (1-ran))
    },mc.cores = nsims)
  names(output)<-vals
  
  output<-data.frame(nbreaks=vals,
                     shift=unlist(output)[seq.int(from=1,to=2L*nsims,by=2)],
                     AbsError=unlist(output)[seq.int(from=2,to=2L*nsims,by=2)])
  
  outs<-c(outs,list(output))
  
}

outing<-data.frame()
for (i in 1:length(outs)){
  
  outing%<>%rbind(outs[[i]])
  
}

p<-ggplot(output)+geom_point(aes(nbreaks,shift,colour=AbsError))+scale_x_log10()
q<-ggplot(output)+geom_point(aes(nbreaks,AbsError,colour=shift))+scale_x_log10()
gridExtra::grid.arrange(p,q,nrow=1)




output2<-mclapply(vals,FUN = function(i) {
  optimise(f=function(shift) IPMShift_optim(shift, ll=i, kn=6700,cost=T),
           lower = 0.2,upper = 0.8)
},mc.cores = nsims)
names(output2)<-vals

output2<-data.frame(nbreaks=vals,
                    shift=unlist(output2)[seq.int(from=1,to=2L*nsims,by=2)],
                    AbsError=unlist(output2)[seq.int(from=2,to=2L*nsims,by=2)])

ggplot(output2)+geom_point(aes(nbreaks,shift,colour=AbsError))+scale_x_log10()
q<-ggplot(output2)+geom_point(aes(nbreaks,AbsError,colour=shift))+scale_x_log10()
gridExtra::grid.arrange(p,q,nrow=1)





output3<-mclapply(vals,FUN = function(i) {
  optimise(f=function(shift) IPMShift_optim(shift, ll=i, kn=7000,cost=T),
           lower = 0.2,upper = 0.8)
},mc.cores = nsims)
output3<-data.frame(nbreaks=vals,
                    shift=unlist(output3)[seq.int(from=1,to=2L*nsims,by=2)],
                    AbsError=unlist(output3)[seq.int(from=2,to=2L*nsims,by=2)])
p<-ggplot(output3)+geom_point(aes(nbreaks,shift,colour=AbsError))+scale_x_log10()
q<-ggplot(output3)+geom_point(aes(nbreaks,AbsError,colour=shift))+scale_x_log10()
gridExtra::grid.arrange(p,q,nrow=1)


# plot(output$nbreaks,output$shift)
# plot(output$nbreaks,output$AbsError)












































funcy<-function(shift,nbreaks){

  # breaks <- seq(1.5, 3.55, l=nbreaks)
  breaks<-quantile(c(simIBM$size[simIBM$census.number==9]),probs = 0:(nbreaks-1)/(nbreaks-1),na.rm = T)%>%unname
  D<-length(breaks)
  sizes <- breaks[-D] + shift*diff(breaks)
  breaks[c(1, D)] <- c(-Inf, Inf)

  sizeDists <- getSizeDistns(simIBM, breaks)
  previousState<-sizeDists[,9]%/%2L
  saveState<-previousState
  
  incr<-0
  # for (j in 1:5){
    
    previousState<-saveState
    
    for (t in 1:3){
      
      nextState<-sampleStateIPM_red(previousState, survFunc, survPars,
                                    growthSamp, growthPars, reprFunc, reprPars, 
                                    offNumSamp, offNumPars, offSizeSamp, offSizePars,
                                    Schild, breaks, oneSex = TRUE, checks = FALSE,
                                    verbose = FALSE, sizes)
      
      incr<-incr+(sum(nextState)/sum(previousState)-popinc)^2
      
      # previousState<-nextState
      
    }
    
  # }

  return(incr)

}

funcy(0.5, nbreaks=10)

optimise(f=function(shift) funcy(shift, nbreaks=10),lower = 0.2,upper = 0.8)






vals<-c(3,4,5,6,7,8,9,10,15,20,25,30,40,50,80,100,150,200,250,300,400,500)
nsims<-length(vals)
output<-mclapply(vals,FUN = function(i) optimise(f=function(shift) funcy(shift, nbreaks=i),
                                                 lower = 0.2,upper = 0.8),mc.cores = nsims)
names(output)<-vals

output<-data.frame(nbreaks=vals,
                   shift=unlist(output)[seq.int(from=1,to=2L*nsims,by=2)],
                   AbsError=unlist(output)[seq.int(from=2,to=2L*nsims,by=2)])
plot(output$nbreaks,output$shift)
plot(output$nbreaks,output$AbsError)


# DF <- data.frame(individual=1:n, size=simIBM$size[simIBM$census.number==1], survived=1, 
#                  census.number=1, parent.size=Start, reproduced=NA, 
#                  prev.size=NA, off.survived=NA, off.born=NA)

checkShift<-function(shift=0.49, ll=6, kn=20){
  
  breaks <- seq(1.5, 3.55, l=ll)
  D<-length(breaks)
  sizes <- breaks[-D] + shift*diff(breaks)
  breaks[c(1, D)] <- c(-Inf, Inf)
  
  sizeDists <- getSizeDistns(simIBM, breaks)
  
  DF<-simIBM%>%filter(census.number==1)
  previous.census <- DF
  previousState<-sizeDists[,1]
  
  OutprevIBM<-OutprevIPM<-OutprevIPMtie<-sizeDists
  prevIBM<-prevIPM<-prevIPMtie<-sizeDists
  OutprevIBM[]<-OutprevIPM[]<-OutprevIPMtie[]<-NA
  
  for (k in 1:kn){
  
    prevIBM[]<-prevIPM[]<-prevIPMtie[]<-NA
    ######## BUILD 3 TYPES - 1) IBM, 2) IPM 3) IPM based on IBM at previous iteration ########
    for (i in 1:max(simIBM$census.number,na.rm = T)){
      
      survivorsDF <- subset(previous.census, previous.census$survived==1)
      current.pop.size <- nrow(survivorsDF)
      D<-length(previousState)
      tiedState <- c(getSizeDistns(survivorsDF, breaks))
      Dt<-length(tiedState)
      
      # Determine who survives this time
      ### IBM ###
      survIBM <- rbinom(current.pop.size, 1, survFunc(survivorsDF$size, survPars) )
      ### IPM ###
      survIPM <- survFunc(sizes, survPars) # This never changes
      
      # Determine their new sizes:
      ### IBM ###
      newIBM<-rep(NA,length(survIBM))
      newIBM[survIBM==1] <- growthSamp(survivorsDF$size[survIBM==1],growthPars)
      ### IPM ###
      newIPM <- survIPM %>% rbinom(n=D, size=previousState) %>%
        rep(x=sizes) %>% growthSamp(growthPars)
      ### IPM - CROSS ###
      newIPMtie <- survIPM %>% rbinom(n=Dt, size=tiedState) %>%
        rep(x=sizes) %>% growthSamp(growthPars)
      
      # Determine who reproduces (out of the survivors):
      #¼¼¼¼¼¼¼¼¼¼ The difference here is that the IBM retains information about the parent-child heritage ¼¼¼¼¼¼¼¼¼¼#
      ### IBM ###
      reprIBMid <- rep(NA, current.pop.size)
      # reprProbs <- reprFunc(newIBM[!is.na(newIBM)], reprPars)
      reprProbs <- reprFunc(na.omit(newIBM), reprPars)
      reprIBMid[survIBM==1] <- rbinom(sum(survIBM==1), 1, reprProbs)
      nRepr <- sum(reprIBMid, na.rm = T)
      # reprIBM <- newIBM[reprIBMid==1 & !is.na(reprIBMid)]
      reprIBM <- survivorsDF$size[reprIBMid==1 & !is.na(reprIBMid)]
      # parentsDF <- subset(survivorsDF, reprIBM==1 & !is.na(reprIBM))
      ### IPM ###
      reprobIPM <- reprFunc(newIPM, reprPars)
      reprIPM <- weightedSelection(newIPM, reprobIPM) # equivalent of parentsDF$size
      ### IPM - CROSS ###
      reprobIPMtie <- reprFunc(newIPMtie, reprPars)
      reprIPMtie <- weightedSelection(newIPMtie, reprobIPMtie) # equivalent of parentsDF$size
      
      # Determine the number of offspring for each parent:
      ### IBM ###
      NoffIBM <- offNumSamp(reprIBM, offNumPars)
      parentIDs <- rep(survivorsDF$individual[reprIBMid==1 & !is.na(reprIBMid)], NoffIBM)
      parentSizes<-rep(reprIBM, NoffIBM)
      ### IPM ###
      NoffIPM <- offNumSamp(reprIPM, offNumPars)
      ### IPM - CROSS ###
      NoffIPMtie <- offNumSamp(reprIPMtie, offNumPars)
      
      # Determine the sizes of the offspring:
      #¼¼¼¼¼¼¼¼¼¼ Difference is from using rep(reprIBM, NoffIBM) instead of NoffIBM directly in offIBM
      ### IBM ###
      offIBM<-rep(NA,length(parentSizes))
      offIBM[!is.na(parentSizes)]<-offSizeSamp(parentSizes[!is.na(parentSizes)], offSizePars)
      offIBM<-offIBM[!is.na(offIBM)]
      ### IPM ###
      offIPM<-offSizeSamp(rep(x = reprIPM, NoffIPM),offSizePars)
      ### IPM - CROSS ###
      offIPMtie<-offSizeSamp(rep(reprIPMtie, NoffIPMtie),offSizePars)
      
      # Determine the survival of the offspring:
      ### IBM ###
      survivingChildren <- rep(0,length(offIBM))
      survivingChildren[!is.na(offIBM)]<-rbinom(length(offIBM[!is.na(offIBM)]), 1, Schild)
      ### IPM ###
      offIPM <- weightedSelection(offIPM,rep(Schild,length(offIPM)))
      ### IPM - CROSS ###
      offIPMtie <- weightedSelection(offIPMtie,rep(Schild,length(offIPMtie)))
      
      # Determine the sex of the offspring:
      #¼¼¼¼¼¼¼¼¼¼ Difference found between onesex and OneGend
      if(OneGend) {
        ### IBM ###
        gender <- rbinom(length(offIBM), 1, 0.5)
        ### IPM ###
        offIPM <- weightedSelection(offIPM,rep(0.5,length(offIPM)))
        ### IPM - CROSS ###
        offIPMtie <- weightedSelection(offIPMtie,rep(0.5,length(offIPMtie)))
      } else gender<-rep(1,length(offIBM))
      ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2222
      
      # Find out which children make it to census: # NOTE: IPMs are already done (above)!
      ### IBM ###
      censusedChildren <- which(gender==1 & survivingChildren==1)
      # offIBM <- offIBM[censusedChildren]
      
      # Update the binned sheep info of current and newly born sheep
      prevIBM[,i]<-vectorToCounts(c(offIBM[censusedChildren], newIBM[survIBM]), breaks)
      prevIPM[,i]<-previousState<-vectorToCounts(c(offIPM, newIPM), breaks)
      prevIPMtie[,i]<-vectorToCounts(c(offIPMtie, newIPMtie), breaks)
      
      # Calculate how many surviving children per parent with a helper:
      helper <- function(x){
        ifelse(!x %in% parentIDs,
               NA,
               sum(x==parentIDs & survivingChildren==1))
      }#helper
      
      helper2 <- function(x){
        ifelse(!x %in% parentIDs,
               NA,
               sum(x==parentIDs))
      }#helper2
      off.survived <- sapply(survivorsDF$individual, helper)
      off.born <- sapply(survivorsDF$individual, helper2)
      
      # Create the DF for the parents:
      currentDF <- data.frame(individual=survivorsDF$individual,
                              size=newIBM,
                              survived=survIBM,
                              census.number=i,
                              parent.size=NA,
                              reproduced=reprIBMid,
                              prev.size=survivorsDF$size,
                              off.survived = off.survived, off.born = off.born)
      # Offspring
      if (length(censusedChildren)==0) {previous.census <- currentDF
      } else {
        newIDStart <- (max(DF$individual)+1)
        newIDs <- newIDStart:(length(censusedChildren)+newIDStart-1)
        offspringDF <- data.frame(individual = newIDs,
                                  size = offIBM[censusedChildren],
                                  survived = 1, census.number = i,
                                  parent.size = parentSizes[censusedChildren],
                                  reproduced = NA, prev.size = NA,
                                  off.survived = NA, off.born = NA)
        
        # Update the previous.census data.frame:
        previous.census <- rbind(currentDF, offspringDF)
      }
      
      DF <- rbind(DF, previous.census)
      
    }
    
    if(k==1){
      OutprevIBM<-prevIBM
      OutprevIPM<-prevIPM
      OutprevIPMtie<-prevIPMtie
    } else {
      OutprevIBM<-OutprevIBM+prevIBM/k
      OutprevIPM<-OutprevIPM+prevIPM/k
      OutprevIPMtie<-OutprevIPMtie+prevIPMtie/k
    }
    
  }
  
  return(list(OutprevIBM=OutprevIBM,OutprevIPM=OutprevIPM,
              OutprevIPMtie=OutprevIPMtie,sizeDists=sizeDists))
  
}




# previous.census<-simIBM[simIBM$census.number==1,]
# iiis<-which(is.na(previous.census$size) & !is.na(previous.census$age))
# iis<-which(is.na(previous.census$size) & is.na(previous.census$age))
# iouts<-which(is.na(simIBM$size) | is.na(simIBM$age))






# lSHEEP<-GetSoaySheep(directory,oneSex=oneSex)
# simIBM<-lSHEEP$solveDF; rm(lSHEEP)
# indys<-simIBM$census.number==1 & is.na(simIBM$size)
# simIBM$size[indys]<-sample(simIBM$size[!(simIBM$census.number==1 | is.na(simIBM$size))],sum(indys),replace = F)
# simIBM$survived[is.na(simIBM$survived)]<-1

simIBM <- do.call(simulateIBM, simParsOrig)
rangeCens<-(max(simIBM$census.number)-50):max(simIBM$census.number)
simIBM<-simIBM[simIBM$census.number%in%rangeCens & 
                  !is.na(simIBM$size),]
simIBM$census.number<-simIBM$census.number-min(simIBM$census.number)+1
simIBM$census.number<-abs(max(simIBM$census.number)-simIBM$census.number)+1
names(simIBM)[1]<-"id"

# indy<-1:nrow(simIBM)
# simIBM<-simIBM[sample(indy,sum(lSHEEP$solveDF$census.number==1),replace = F),]
# simIBM$census.number<-1
# names(simIBM)[1]<-"id"
# gpmodel<-GauPro::GauPro(simIBM$age[-iouts],simIBM$size[-iouts], parallel=FALSE)
# plot(previous.census$age[-iis],previous.census$size[-iis])
# lines(
#   rnorm(n = 12*15,gpmodel$predict(1:11),
#         gpmodel$predict(previous.census$age[iiis], se=T)$se/sqrt(length(previous.census$age[-iouts])))
# 
#             , add=T, col=2)

# FOR THE REAL SHEEP DATASET - CONTAINS NA VALUES
checkShift_optim<-function(shift=0.49, ll=6, kn=50,cost=T){
  
  breaks<-quantile(c(simIBM$size),probs = 0:(ll-1)/(ll-1),na.rm = T)%>%unname
  # breaks<-seq(from=min(simIBM$size,na.rm = T),to=max(simIBM$size,na.rm = T),length.out=ll)
  D<-length(breaks)
  sizes <- breaks[-D] + shift*diff(breaks)
  breaks[c(1, D)] <- c(-Inf, Inf)
  
  sizeDists <- getSizeDistns(simIBM, breaks)
  
  DF<-simIBM[simIBM$census.number==1,]
  previous.census <- DF
  previousState<-sizeDists[,1]

  OutprevIBM<-OutprevIPM<-prevIBM<-prevIPM<-sizeDists
  OutprevIBM[]<-OutprevIPM[]<-NA

  # isamps<-is.na(previous.census$size)
  # isizes<-!is.na(simIBM$size)
    
  for (k in 1:kn){
    
    prevIBM[]<-prevIPM[]<-NA
    previous.census <- DF
    previousState<-sizeDists[,1]
    
    # # Randomly sample sizes for missing size sheep from the population
    # previous.census$size[isamps]<-sample(simIBM$size[isizes],sum(isamps),replace = T)
    # # Assume all the NA sheep were still alive but just not weighed
    # previous.census$survived[isamps]<-1
  
    ######## BUILD 3 TYPES - 1) IBM, 2) IPM 3) IPM based on IBM at previous iteration ########
    for (i in 1:max(simIBM$census.number,na.rm = T)){
      
      # survivorsDF <- previous.census[previous.census$survived==1 & !is.na(previous.census$survived),]
      survivorsDF <- subset(previous.census, previous.census$survived==1)
      
      
      
      
      current.pop.size <- nrow(survivorsDF)
      D<-length(previousState)
      tiedState <- c(getSizeDistns(survivorsDF, breaks))
      Dt<-length(tiedState)

      ### IBM ###
      newIBM<-survIBM<-rep(NA,current.pop.size)
      survIBM[!is.na(survivorsDF$size)] <- rbinom(current.pop.size, 1, 
                                                  survFunc(survivorsDF$size[!is.na(survivorsDF$size)], survPars) )
      newIBM[survIBM==1 & !is.na(survIBM)] <- growthSamp(survivorsDF$size[survIBM==1 & !is.na(survIBM)],growthPars)
      ### IPM ###
      newIPM <-  growthSamp(rep(rbinom(survFunc(sizes, survPars),
                            n=D, size=previousState),x=sizes),growthPars)
      
      # Determine who reproduces (out of the survivors):
      #¼¼¼¼¼¼¼¼¼¼ The difference here is that the IBM retains information about the parent-child heritage ¼¼¼¼¼¼¼¼¼¼#
      ### IBM ###
      reprIBMid <- rep(NA, current.pop.size)
      reprIBMid[survIBM==1 & !is.na(survIBM)] <- rbinom(sum(survIBM,na.rm = T), 1, reprFunc(na.omit(newIBM), reprPars))
      # nRepr <- sum(reprIBMid, na.rm = T)
      # reprIBM <- newIBM[reprIBMid==1 & !is.na( reprIBMid)]
      iiis<-reprIBMid==1 & !is.na(reprIBMid)
      reprIBM <- survivorsDF$size[iiis]
      # parentsDF <- subset(survivorsDF, reprIBM==1 & !is.na(reprIBM))
      ### IPM ###
      reprobIPM <- reprFunc(newIPM, reprPars)
      reprIPM <- weightedSelection(newIPM, reprobIPM) # equivalent of parentsDF$size
      
      # Determine the number of offspring for each parent:
      ### IBM ###
      NoffIBM <- offNumSamp(reprIBM, offNumPars)
      # parentIDs <- rep(survivorsDF$id[iiis], NoffIBM)
      
      # Determine the sizes of the offspring:
      #¼¼¼¼¼¼¼¼¼¼ Difference is from using rep(reprIBM, NoffIBM) instead of NoffIBM directly in offIBM
      ### IBM ###
      offIBM<-offIPM<-rep(NA,length(reprIBM))
      
      
      
      
      
      
      
      offIBM[!is.na(reprIBM)]<-offSizeSamp(reprIBM[!is.na(reprIBM)], offSizePars)
      # This is same length as the parent sizes, not the number of offspring actually born (important when offNum!=1)
      
      
      
      
      
      
      
      
      
      
      
      offIBM<-offIBM[!is.na(offIBM)]
      ### IPM ###
      offIPM<-offSizeSamp(reprIPM[!is.na(reprIPM)],offSizePars)
      offIPM<-offIPM[!is.na(offIPM)]
      
      # Determine the survival of the offspring:
      ### IBM ###
      survivingChildren <- rep(0,length(offIBM))
      survivingChildren[!is.na(offIBM)]<-rbinom(length(offIBM[!is.na(offIBM)]), 1, Schild)
      ### IPM ###
      offIPM <- weightedSelection(offIPM,rep(Schild,length(offIPM)))
      
      # Determine the sex of the offspring:
      #¼¼¼¼¼¼¼¼¼¼ Difference found between onesex and OneGend
      if(OneGend) {
        ### IBM ###
        gender <- rbinom(length(offIBM), 1, 0.5)
        ### IPM ###
        offIPM <- weightedSelection(offIPM,rep(0.5,length(offIPM)))
      } else gender<-rep(1,length(offIBM))
      ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2222
      
      # Find out which children make it to census: # NOTE: IPMs are already done (above)!
      ### IBM ###
      censusedChildren <- which(gender==1 & survivingChildren==1)
      # offIBM <- offIBM[censusedChildren]
      
      # Update the binned sheep info of current and newly born sheep
      prevIBM[,i]<-vectorToCounts(c(offIBM[censusedChildren], newIBM[survIBM]), breaks)
      prevIPM[,i]<-previousState<-vectorToCounts(c(offIPM, newIPM), breaks)
      
      # Create the DF for the parents:
      currentDF <- data.frame(id=survivorsDF$id,
                              size=newIBM)
      # Offspring
      if (length(censusedChildren)==0) {previous.census <- currentDF
      } else {
        newIDStart <- (max(previous.census$id)+1)
        newIDs <- newIDStart:(length(censusedChildren)+newIDStart-1)
        offspringDF <- data.frame(id = newIDs,
                                  size = offIBM[censusedChildren])
        
        # Update the previous.census data.frame:
        previous.census <- rbind(currentDF, offspringDF)
      }
      
    }
    
    if(k==1){
      OutprevIBM<-prevIBM
      OutprevIPM<-prevIPM
    } else {
      OutprevIBM<-OutprevIBM+prevIBM/k
      OutprevIPM<-OutprevIPM+prevIPM/k
    }
    
  }
  
  if(cost) return(sum(abs(colSums(OutprevIBM)-colSums(OutprevIPM))))
  
  return(list(OutprevIBM=OutprevIBM,OutprevIPM=OutprevIPM,
              sizeDists=sizeDists))
  
}



# funcyCheck<-function(shift,ll) {
#   CostF<-0
#   
#   for(j in 1:30){
#     
#     output<-checkShift_optim(shift,ll,kn = 1)
#     print(paste0("SizeDists: ",sum(colSums(output$sizeDists))))
#     print(paste0("prevIBM: ",sum(colSums(output$OutprevIBM))," compare ",100/sum(output$sizeDists)*sum(abs(colSums(output$OutprevIBM)-colSums(output$sizeDists)))))
#     print(paste0("prevIPM: ",sum(colSums(output$OutprevIPM))," compare ",100/sum(output$sizeDists)*sum(abs(colSums(output$OutprevIPM)-colSums(output$sizeDists)))))
#     print(paste0("prevIPMtie: ",sum(colSums(output$OutprevIPMtie))," compare ",100/sum(output$sizeDists)*sum(abs(colSums(output$OutprevIPMtie)-colSums(output$sizeDists)))))
#     
#     if(j==1){
#       OutprevIBM<-output$OutprevIBM
#       OutprevIPM<-output$OutprevIPM
#       OutprevIPMtie<-output$OutprevIPMtie
#     } else {
#       OutprevIBM<-OutprevIBM+output$OutprevIBM/j
#       OutprevIPM<-OutprevIPM+output$OutprevIPM/j
#       OutprevIPMtie<-OutprevIPMtie+output$OutprevIPMtie/j
#     }
#     
#   }
#   
#   return(list(OutprevIBM=OutprevIBM,OutprevIPM=OutprevIPM,
#               OutprevIPMtie=OutprevIPMtie,sizeDists=sizeDists))
#   
# }

# output<-checkShift_optim(0.49,8)
# 
# funcy<-function(shift,ll,nnn=30) {
#   CostF<-0
#   for(j in 1:nnn){
#     output<-checkShift(shift,ll,kn = 1)
#     CostF<-CostF+sum(abs(colSums(output$OutprevIBM)-colSums(output$OutprevIPM)))
#     print(paste0("Shift = ",shift,", Cost = ",CostF))
#   }
#   return(CostF/nnn)
# }


# out5<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=5,nnn=75)
# out6<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=6,nnn=75)
# out7<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=7,nnn=75)
# out8<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=8,nnn=75)
# out9<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=9,nnn=75)
# out10<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=10,nnn=75)
# out15<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=15,nnn=75)
# out20<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=20,nnn=75)
# out30<-optimise(f=funcy,lower = 0.35,upper = 0.65,ll=30,nnn=75)
# out50<-optimise(f=funcy,lower = 0.25,upper = 0.5,ll=50,nnn=75)
# out100<-optimise(f=funcy,lower = 0.25,upper = 0.5,ll=100,nnn=75)
# out200<-optimise(f=funcy,lower = 0.25,upper = 0.5,ll=200,nnn=75)

IBM<-IPM<-c()
# for (j in 1:100){
#   output<-checkShift(out10$minimum,10,1)
#   IBM<-IBM+sum(output$OutprevIBM)
#   IPM<-IPM+sum(output$OutprevIPM)
# }

funcyError<-function(i){
    output<-checkShift(out10$minimum,10,1)
    return(c(sum(output$OutprevIBM),sum(output$OutprevIPM)))
}

nsims<-150

output<-mclapply(1:nsims,FUN = funcyError,mc.cores = 32)
output<-data.frame(IBM=unlist(output)[seq.int(from=1,to=2L*nsims,by=2)],
           IPM=unlist(output)[seq.int(from=2,to=2L*nsims,by=2)])

mean(output$IPM-output$IBM)
plot(output$IPM,output$IBM)
# mean(IBM-sum(output$sizeDists))
# mean(IPM-sum(output$sizeDists))


plot(colSums(prevIBM),colSums(sizeDists)); abline(b=1,a=0)
plot(colSums(prevIPM),colSums(sizeDists)); abline(b=1,a=0)
plot(colSums(prevIBM),colSums(prevIPMtie)); abline(b=1,a=0)
plot(colSums(prevIBM),colSums(prevIPM)); abline(b=1,a=0)
plot(colSums(prevIPM),colSums(prevIPMtie)); abline(b=1,a=0)




