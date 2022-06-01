# vals<-unlist(x0)

x0<-vals<-c(-8.25, 3.77,
            1.41, 0.56, log(0.08),
            -7.23, 2.6,
            log(0.06),
            0.36, 0.71, log(0.16),
            qlogis(0.9),
            log(50),
            log(10)
)
if(fixedObsProb) x0<-x0[1:12]
# Calculate the expected observation probability
obsMean<-vals[13]/(vals[13]+vals[14])

if(!manshift) {shift<-CalcShift_Kernel(vals,IPMLTP,nbks,oneSex,lSHEEP$L,lSHEEP$U)
} else shift<-0.5

for (i in 1:length(IPMLTP$links))  vals[i] <- match.fun(IPMLTP$links[i])(vals[i])
vals%<>%relist(skeleton=IPMLTP$skeleton)

if(normsampler=="sampleDTN") {
  vals$growthPars%<>%c(lSHEEP$L,lSHEEP$U)
  vals$offSizePars%<>%c(lSHEEP$L,lSHEEP$U)
}

if(fixedObsProb) { lSHEEP$obsProbTime <- rep(obsMean,yearing+1)
} else lSHEEP$obsProbTime <- rbeta(yearing+1,vals$obsProbPar[1],vals$obsProbPar[2])

simPars <- list(n=100, t=300,
                # set survival details:
                survFunc = IPMLTP$survFunc, survPars = vals$survPars,
                # set growth details:
                growthSamp = IPMLTP$growthSamp,
                growthPars = vals$growthPars,
                # set reproduction details:
                reprFunc = IPMLTP$reprFunc, reprPars = vals$reprPars,
                # set offspring number and size distribution details:
                offNumSamp=IPMLTP$offNumSamp, offNumPars=vals$offNumPars,
                offSizeSamp = IPMLTP$offSizeSamp,
                offSizePars = vals$offSizePars,
                # Child survival probability:
                Schild=vals$Schild,
                # set other miscelaneous parameters:
                Start = 2.7, thresh=10000, OneGend = TRUE,
                popPrint = F, verbose=F)

simmedData <- do.call(simulateIBM, simPars)

# simPars$t<-yearing
# simPars$Start<-sample(simmedData$size,poptot)
# simPars$n<-poptot
# simmedData <- do.call(simulateIBM, simPars)  

Starts<-sample(simmedData$size[!is.na(simmedData$size) & simmedData$census.number>round(0.5*max(simmedData$census.number))],poptot)

lSHEEP$breaks<-ggplot_build(ggplot()+geom_histogram(data=data.frame(size=Starts),mapping=aes(size),bins = nbks))$data[[1]]$x
lSHEEP$sizes <- lSHEEP$breaks[-(nbks)] + shift*diff(lSHEEP$breaks)

previousState<-vectorToCounts(Starts, lSHEEP$breaks)

minnie<-lSHEEP$breaks[min(which(cumsum(previousState)/sum(previousState)>0.05))]
others<-seq(from=minnie,to=max(lSHEEP$breaks),length.out=(nbks-1))
lSHEEP$breaks<-c(min(lSHEEP$breaks),others)
lSHEEP$sizes <- lSHEEP$breaks[-(nbks)] + shift*diff(lSHEEP$breaks)

previousState<-vectorToCounts(Starts, lSHEEP$breaks)

stateSpaceSampArgs <- list(survFunc = IPMLTP$survFunc, survPars = vals$survPars,
                           growthSamp = IPMLTP$growthSamp, growthPars = vals$growthPars,
                           reprFunc = IPMLTP$reprFunc, reprPars = vals$reprPars,
                           offNumSamp = IPMLTP$offNumSamp, offNumPars = vals$offNumPars,
                           offSizeSamp = IPMLTP$offSizeSamp, breaks = lSHEEP$breaks,
                           offSizePars = vals$offSizePars, Schild=vals$Schild,
                           sizes=lSHEEP$sizes, oneSex = oneSex)

lSHEEP$COUNTS<-projectStateSpace(sampleStateIPM_red, stateSpaceSampArgs, previousState, yearing)
lSHEEP$breaks[c(1,nbks)]<-c(-Inf,Inf)

# lSHEEP$COUNTS<-matrix(NA,nrow = (nbks-1), ncol = yearing)
# for(t in 1:yearing){
#   
#   Weight<-simmedData$size[sample(which(simmedData$census.number==t),
#                                  round(sum(simmedData$census.number==t)*lSHEEP$obsProbTime[t]),replace = F)]
#   lSHEEP$COUNTS[,t]<-vectorToCounts(Weight, lSHEEP$breaks)
#   
# }
lSHEEP$priorProbs<-rowSums(lSHEEP$COUNTS)/sum(lSHEEP$COUNTS)
saveRDS(list(vals=vals,lSHEEP=lSHEEP),paste0("Results/SimulatedData_",namer,".Rdata"))

# popdyn<-eigen(kernelOneVar(m = 1000, growthFunc = IPMLTP$growthFunc,
#                      growthPars = vals$growthPars, survFunc = IPMLTP$survFunc,
#                      survPars = vals$survPars, repFunc = IPMLTP$reprFunc,
#                      repPars = vals$reprPars, offNum = vals$offNumPars,
#                      offSizeFunc =  IPMLTP$offSizeFunc,
#                      offSizePars = vals$offSizePars, L = L, U = U,
#                      childSurv = vals$Schild,shift=0.5,halfPop = T) )
# 
# popinc<-Re(popdyn$values[1])
# 
# lSHEEP$COUNTS<-matrix(NA,nrow = (nbks-1), ncol = yearing)
# 
# for(t in 1:yearing){
#   
#   Weight<-sample(seq(lSHEEP$L,lSHEEP$U,length.out=1000),size=(poptot*popinc^t),replace = T,prob = rowsums(abs(popdyn$vectors)))  
#   # obsprob remove some!
#   Weight<-Weight[sample(1:length(Weight),round(length(Weight)*lSHEEP$obsProbTime[t]),replace = F)]
#   
#   # Need counts (eigenvector sampled and then binned)
#   lSHEEP$COUNTS[,t]<-vectorToCounts(Weight, lSHEEP$breaks)
# 
# }
# 
# lSHEEP$priorProbs<-rowSums(lSHEEP$COUNTS)/sum(lSHEEP$COUNTS)
# saveRDS(lSHEEP,paste0("Results/SimulatedData_",namer,".Rdata"))