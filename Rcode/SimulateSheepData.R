lSHEEP<-list(L=1.131, U=3.532)
# for truncated normal distributions
if(normsampler=="sampleDTN") {
  IPMLTP$DTN<-data.frame(L=lSHEEP$L,U=lSHEEP$U)
}

x0<-vals<-c(-7.25, 3.77,
            1.49111883, 0.53069364, log(0.08806918),
            -4.619443,  1.369697,
            log(0.067204),
            0.6887019, 0.5931042, 0.2073659,
            qlogis(0.4989097),
            log(50),
            log(10))
Np<-length(unlist(x0))

if(fixedObsProb) {
  # Remove the last two parameters
  x0<-x0[1:12]; Np<-length(unlist(x0))
  # Calculate the expected observation probability
  obsMean<-vals[13]/(vals[13]+vals[14])
  # Temporal observation probability must be kept fixed
  lSHEEP$obsProbTime <- rep(obsMean,yearing+1)
} else lSHEEP$obsProbTime <- rbeta(yearing+1,vals$obsProbPar[1],vals$obsProbPar[2])

# Convert to physical coordinates
vals%<>%Sample2Physical(IPMLTP)

# Stable population values and distribution from the eigenvalues & vectors
popmod<-kernelOneVar(m = 500, growthFunc = IPMLTP$growthFunc,
                     growthPars = vals$growthPars, survFunc = IPMLTP$survFunc,
                     survPars = vals$survPars, repFunc = IPMLTP$reprFunc,
                     repPars = vals$reprPars, offNum = vals$offNumPars,
                     offSizeFunc =  IPMLTP$offSizeFunc,
                     offSizePars = vals$offSizePars, L = lSHEEP$L, U = lSHEEP$U,
                     childSurv = vals$Schild,shift=0.5,halfPop = 0.5)%>%eigen

popinc<-Re(popmod$values[1])
eigvec<-data.frame(size=seq(lSHEEP$L,lSHEEP$U,length.out=500),prob=-Re(popmod$vectors[,1]))

###################### SIMULATE THE INDIVIDUAL BASED MODEL ######################

Starties<-sample(eigvec$size,poptot,prob = abs(eigvec$prob),replace = T)

simPars <- list(n=poptot, t=yearing,
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
                Schild=vals$Schild, obsProb=lSHEEP$obsProbTime,
                # set other miscelaneous parameters:
                Start = Starties, thresh=10000, OneGend = TRUE,
                popPrint = F, verbose=F)

# Let's simulate some sheeepieeeesss!
lSHEEP <- do.call(simulateIBM, simPars)

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

x0<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb)))%>%relist(skeleton = skeleton)
propCOV<-2*diag(unlist((do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb,CI=T))))$sd))

# x0<-x0+runif(length(x0),-3,3)

#################################################################################

###################### SIMULATE THE STATE SPACE MODEL ######################
# breaks<-lSHEEP$breaks; breaks[c(1,nbks)]<-c(-Inf,Inf)
# previousState<-vectorToCounts(Starties, breaks)
# 
# stateSpaceSampArgs <- list(survFunc = IPMLTP$survFunc, survPars = vals$survPars,
#                            growthSamp = IPMLTP$growthSamp, growthPars = vals$growthPars,
#                            reprFunc = IPMLTP$reprFunc, reprPars = vals$reprPars,
#                            offNumSamp = IPMLTP$offNumSamp, offNumPars = vals$offNumPars,
#                            offSizeSamp = IPMLTP$offSizeSamp, breaks = lSHEEP$breaks,
#                            offSizePars = vals$offSizePars, Schild=vals$Schild,
#                            sizes=lSHEEP$sizes, oneSex = oneSex)
# 
# simmies<-50
# tmp<-lapply(1:simmies,function(i) projectStateSpace(sampleStateIPM_red, stateSpaceSampArgs, previousState, yearing))
# counter<-tmp[[1]]
# for(i in 2:simmies) counter<-counter+tmp[[i]]
# 
# check<-round(counter/simmies)
# diff(colSums(check)/100)+1
# popinc
# 
# poptrack<-data.frame()
# for(shsh in seq(0.05,0.95,length.out=200)){
#   breaks<-lSHEEP$breaks
#   sizes <- breaks[-(nbks)] + shsh*diff(breaks)
#   breaks[c(1,nbks)]<-c(-Inf,Inf)
#   
#   stateSpaceSampArgs <- list(survFunc = IPMLTP$survFunc, survPars = vals$survPars,
#                              growthSamp = IPMLTP$growthSamp, growthPars = vals$growthPars,
#                              reprFunc = IPMLTP$reprFunc, reprPars = vals$reprPars,
#                              offNumSamp = IPMLTP$offNumSamp, offNumPars = vals$offNumPars,
#                              offSizeSamp = IPMLTP$offSizeSamp, breaks = breaks,
#                              offSizePars = vals$offSizePars, Schild=vals$Schild,
#                              sizes=sizes, oneSex = oneSex)
#   
#   tmp<-lapply(1:simmies,function(i) projectStateSpace(sampleStateIPM_red, stateSpaceSampArgs, previousState, yearing-1))
#   counter<-tmp[[1]]
#   for(i in 2:simmies) counter<-counter+tmp[[i]]
#   poptrack%<>%rbind(data.frame(year=0:(yearing-1), pop=colSums(round(counter/simmies)),shift=shsh))
# }
# ggplot(poptrack,aes(year,pop,group=shift))+geom_line(aes(colour=shift))+
#   geom_line(data=data.frame(year=0:(yearing-1),pop=colSums(COUNTS)),aes(year,pop),colour="red")
#################################################################################



# cumy<-cumsum(abs(eigvec$prob)); cumy<-cumy/max(cumy)

# # Setup the size classes based upon the eigenvector
# if(regbinspace){
#   # Let's make sure that the smallest bin contains 1/nbrks of the population
#   minbin<-eigvec$size[which.min(abs(cumy-1/nbks))]
#   lSHEEP$breaks<-seq(minbin,lSHEEP$U-1e-5,length.out=nbks)
# } else {
#   lSHEEP$breaks<-sapply(0:(nbks-1)/(nbks-1),function(quant) eigvec$size[which.min(abs(quant-cumy))])
# }
# # Calculate the ideal shift in the size classes for the state space projection
# if(!manshift) {shift<-CalcShift_Kernel(x0 = vals, IPMLTP = IPMLTP,nbks = nbks,
#                                        breaks = lSHEEP$breaks,
#                                        halfpop = oneSex,L = lSHEEP$L,U = lSHEEP$U)
# } else shift<-0.5
# # Now setup the ideal class sizes
# lSHEEP$sizes <- lSHEEP$breaks[-(nbks)] + shift*diff(lSHEEP$breaks)
# 
# lSHEEP$COUNTS<-sapply(1:yearing, function(t) vectorToCounts(sample(eigvec$size,
#                                                                    poptot+round(poptot*(popinc-1))*(t-1),
#                                                                    prob = abs(eigvec$prob),replace = T),lSHEEP$breaks))
# lSHEEP$cNames<-c(vapply(paste0("yr_",as.character(1:yearing),"__"), 
#                      function(yname) paste0(yname,paste0("sz_",as.character(signif(lSHEEP$sizes,4)))),
#                      character(length(lSHEEP$sizes))))
# 
# lSHEEP$priorProbs<-rowSums(lSHEEP$COUNTS)/sum(lSHEEP$COUNTS)
# saveRDS(list(vals=vals,lSHEEP=lSHEEP),paste0("Results/SimulatedData_",namer,".Rdata"))






# sample the sheep population t times, each time with a total:
# poptot=poptot+round(poptot*(1-popinc))*(t-1)


# Separate real from simulated datasets in Main.R











# lSHEEP$COUNTS<-matrix(NA,nrow = (nbks-1), ncol = yearing)
# for(t in 1:yearing){
#   
#   Weight<-simmedData$size[sample(which(simmedData$census.number==t),
#                                  round(sum(simmedData$census.number==t)*lSHEEP$obsProbTime[t]),replace = F)]
#   lSHEEP$COUNTS[,t]<-vectorToCounts(Weight, lSHEEP$breaks)
#   
# }

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