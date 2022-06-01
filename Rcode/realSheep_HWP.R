# directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
directory<-paste0(getwd(),"/")

list.of.packages <- c("xtable","magrittr","doParallel","Rfast","mc2d", 
                      "abind")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(new.packages)
if(length(new.packages)>0) install.packages(new.packages)

# if(length(list.of.packages[!("mvtnorm" %in% installed.packages()[,"Package"])])){devtools::install_github("rCarto/osrm")}

source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))
source(paste0(directory,'Rcode/SoaySheepData.R'))
library(dissPackage3)
library(xtable)
# library(mcmcse)
library(tidyverse)
library(magrittr)

simulation<-T
poptot<-100
yearing<-10
# Is the population counted one sex or two?
oneSex<-T
# Is the observation probability an empirically-based fixed value or sampled as a R.V.?
fixedObsProb<-F
# Read in the Soay sheep data
lSHEEP<-GetSoaySheep(directory,oneSex=oneSex)
# Number of MCMC simulations
itermax <- 10000
# Number of in-chain parallelised cores
ncores<-1
# Define the number of size class bins
nbks<-10
# Particle filter initialisation function
muModel<-'poisson' # 'multinomial'
# Observation Model
obsModel<-'binomial' #'binomial' # 'multinomial' #'poisson'
manshift<-T
# For the individual and offspring growth function - normal or truncated normal, or otherwise?
normsampler<-"sampleDTN" # 'sampleNorm'

namer<-paste0(ifelse(simulation,paste0("SIM_pop",poptot,"_yr",yearing),"REAL"),"_GSF_",ifelse(fixedObsProb,"fixed","beta"),"_",muModel,"Mu_",obsModel,"Obs_GLMx0_",itermax,"_",nbks,"brks_",normsampler,"_",ifelse(manshift,"manshift","autoshift"))

# Skeleton frame for the parameterisation vector
skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA
)


# Link functions to be used
returnSelf <- function(x) x
linkNum <- function(x) exp(x)+1
if(oneSex) {Schilder <- function(x) 0.5*plogis(x)
} else {Schilder <- function(x) plogis(x)}
links<-c(
  'returnSelf', # Survival Logistic Regression Intercept
  'returnSelf', # Survival Logistic Regression Gradient
  'returnSelf', # Growth Linear Regression Intercept
  'returnSelf', # Growth Linear Regression Gradient
  'exp', # Growth Linear Regression Dispersion (Std. Dev.)
  'returnSelf', # Reproduction Logistic Regression Intercept
  'returnSelf', # Reproduction Logistic Regression Gradient
  'linkNum', # Offspring Number per Birth
  'returnSelf', # Offspring Size Linear Regression Intercept
  'returnSelf', # Offspring Size Linear Regression Gradient
  'exp', # Offspring Size Linear Regression Dispersion (Std. Dev.)
  'Schilder' # Offspring Survival Probability
)

# The start of the list of functions, parameters and formatting for the logTarget
IPMLTP <- list(
  skeleton = skeleton,
  links = links,
  survFunc = match.fun('linLogit'),
  growthSamp = match.fun(normsampler),
  reprFunc = match.fun('linLogit'),
  offNumSamp = match.fun('PoisNum'),
  offSizeSamp = match.fun(normsampler)
)

#if(normsampler=="sampleDTN") {
  IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- doublyTruncatedNormal
#}else growthFunc <- offSizeFunc <- normal

L<-min(c(lSHEEP$solveDF$prev.size, lSHEEP$solveDF$size),na.rm = T)
U<-max(c(lSHEEP$solveDF$prev.size, lSHEEP$solveDF$size),na.rm = T)

# Make the inital values using the piecemeal GLM approach
x0<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb)))
# x0<-readRDS(paste0(directory,"/Results/x0_GSF_fixed_multMu_multObs_GLMx0"))
# x0<-readRDS("./Results/x0_GSF_fixed_multMu_multObs_GLMx0_nbks15_autoshift")
# x0<-c(tx0,x0[(length(x0)-1):length(x0)])

# Number of parameters
Np<-length(x0)

# Import covariance matrix:
# propCOV<-readRDS("./Results/propCOV_GSF_fixed_multMu_multObs_GLMx0_nbrks15_autoshift")
# propCOV<-matrix(0,nrow = length(x0),ncol = length(x0))
# diag(propCOV)<-length(x0)/1000
# propCOV[1:nrow(tpropCOV),1:nrow(tpropCOV)]<-tpropCOV
propCOV<-diag(Np)/60
# diag(propCOV)[c(4,13,14)]<-0.001
# diag(propCOV)[13:14]<-c(1e-3,1e-3)
# propCOV<-readRDS(paste0(directory,"/Results/propCOV_GSF_fixed_multMu_multObs_GLMx0"))
# propCOV<-readRDS(paste0(directory,"/Results/propCOV_fixedObsP_MH_GSF_poisO_Mult"))

# Observed Probability Beta Shape Param 1 & 2
if(!fixedObsProb) IPMLTP$links%<>%c('exp','exp')
# Make sure that the skeleton frame also includes this
if(!fixedObsProb) IPMLTP$skeleton %<>% c(list(obsProbPar = rep(NA,2)))
x0%<>%relist(skeleton = IPMLTP$skeleton)
# Calculate the shift factor that offsets the size class mid-point
if(!manshift) {shift<-CalcShift_Kernel(x0,IPMLTP,nbks,oneSex,L,U)
} else shift<-0.5
print(paste0("Grid shift = ",shift, " for ",nbks," number of breaks." ))
# Get the sheep counts and sizes
lSHEEP<-GetSoaySheep_binned(lSHEEP,shift=shift,oneSex=T,nbks=nbks)

if(simulation) {
    
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
  
  if(!manshift) {shift<-CalcShift_Kernel(vals,IPMLTP,nbks,oneSex,L,U)
  } else shift<-0.5

  for (i in 1:length(IPMLTP$links))  vals[i] <- match.fun(IPMLTP$links[i])(vals[i])
  vals%<>%relist(skeleton=IPMLTP$skeleton)
  
  if(normsampler=="sampleDTN") {
    vals$growthPars%<>%c(L,U)
    vals$offSizePars%<>%c(L,U)
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
  #   Weight<-sample(seq(L,U,length.out=1000),size=(poptot*popinc^t),replace = T,prob = rowsums(abs(popdyn$vectors)))  
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
  
}

# Housekeeping
IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- NULL; x0%<>%unlist()
lSHEEP$solveDF<-NULL

# Selecting dependent upon configuration
if(obsModel=='multinomial'){
  # Multinomial observational model
  if(fixedObsProb) obsfun<-match.fun('fixedMuObs')
  else obsfun<-match.fun('beta_mnomObs')
  # Number of SMC particles
  SMC_parts<-6700  # 6700 # 200
  
} else if(obsModel=='poisson'){
  # Poisson observational model
  if(fixedObsProb) obsfun<-match.fun('fixedPoisObs')
  else obsfun<-match.fun('beta_poisObs')
  # Number of SMC particles
  SMC_parts<-400  # 6700 # 200
  
} else if(obsModel=='binomial'){
  # Binomial observational model
  if(fixedObsProb) obsfun<-match.fun('fixedBinObs')
  else obsfun<-match.fun('beta_binObs')
  # Number of SMC particles
  SMC_parts<-1500  # 6700 # 200
  
} else stop("Error in Mu Model input, we recommend muModel<-'poissonMu' ")

# Particle filter initialisation parameters
if(muModel=='multinomial'){
  # Multinomial Mu initialisation function
  muPar<-list(popSize=sum(lSHEEP$COUNTS[,1]), n=SMC_parts, probs=lSHEEP$priorProbs)
  if(fixedObsProb) mufun<-match.fun("multnomMu")
  else mufun<-match.fun("beta_mnMu")
} else if(muModel=='poisson'){
  # Poisson Mu initialisation function
  muPar<-list(popSize=lSHEEP$COUNTS[,1],n=SMC_parts)
  if(fixedObsProb) mufun<-match.fun("poissonMu")
  else mufun<-match.fun("beta_poisMu")
} else stop("Error in Mu Model input, we recommend muModel<-'poissonMu' ")

# If we want to avoid parameterising the probability of being observed, we can extract it empirically straight from the data
if(fixedObsProb){
  obsProbTime<-lSHEEP$obsProbTime
  muPar$pobs<-obsProbTime[1]
}

########################### SPECIFYING THE PRIORS ##############################

priorName<-"uniform"

priorsIPM<-switch(priorName,uniform=rep("dunif", Np),cauchy=rep("dcauchy", Np),normal=rep("dnorm", Np),stop("Prior distribution not recognised"))

# Setup the prior names:
listy<-switch(priorName,uniform=list(min=-50, max=50),cauchy=list(location=0,scale=0.15),normal=list(mean=10,sd=10),stop("Prior distribution not recognised"))
# Setup the prior hyper par values:
flatPriors <- list(
  # Survival function intercept:
  sFuncInt=listy,
  # Survival function gradient:
  sFuncGra=listy,
  # Growth function intercept:
  gFuncInt=listy,
  # Growth function gradient:
  gFuncGra=listy,
  # log(Growth function sd):
  gFuncSD=listy,
  # Reproduction function intercept:
  rFuncInt=listy,
  # Reproduction function gradient:
  rFuncGra=listy,
  # Number of offspring per litter;
  offNum=listy,
  # Childsize intercept:
  osFuncInt=listy,
  # Childsize gradient:
  osFuncGra=listy,
  # Childsize sd:
  osFuncSD=listy,
  # qlogis(Child survival):
  Schild=listy
)
if(!fixedObsProb) flatPriors%<>%
  c(list(  
    # log(Shape 1 of beta dist for prob of detection of an animal):
    obsProbSh1=listy,
    # log(Shape 2 of beta dist for prob of detection of an animal):
    obsProbSh2=listy
  ))

#################### CREATE THE LOGTARGETPARAMETERS OBJECT #####################
# readRDS(paste0(directory,"RDSobjects/IPMLTP"))

IPMLTP %<>% c(list(
  priorFunc = match.fun('evalPriors'),
  priors = priorsIPM,
  priorPars = flatPriors,
  oneSex = oneSex,
  mu = mufun, # mu = match.fun('multnomMu'), 
  muPar = muPar, # muPar = list(popSize=sum(lSHEEP$COUNTS[,1]), n=SMC_parts, probs=priorProbs),
  b = SMC_parts, 
  Y = lSHEEP$COUNTS,
  obsProb = obsfun, # obsProb = match.fun('detectionNumObs'), # obsProb = match.fun('poissonObs'),
  fixedObsProb=fixedObsProb,
  breaks = lSHEEP$breaks,
  sizes = lSHEEP$sizes
))
if(fixedObsProb) IPMLTP %<>% c(list(obsProbPar = obsProbTime))
#if(normsampler=="sampleDTN") IPMLTP %<>% c(list(DTN = c(L,U)))
IPMLTP %<>% c(list(DTN = c(L,U)))

print("Initial Values Log-Likelihood=")
print(logTargetIPM(x0, logTargetPars = IPMLTP, returnNeg = T, printProp = F))

########################### SIMULATED ANNEALING ################################
# set.seed(101010)
# # get the hessian:
# SANN <-  optim(startValues, logTargetIPM, logTargetPars = IPMLTP, returnNeg = T,
#                method = "SANN", control = list(maxit = 10), printProp = T,
#                hessian = T)
# saveRDS(SANN, "../RDSobjects/HWP/Uniform/SANNuniform")
# # get the covariance of the proposal:
# propCov <- solve(SANN$hessian)
# 
# set.seed(101010)
# IPMLTP2 <- IPMLTP
# IPMLTP2$b <- 1000
# # get the MLEs by using more iterations:
# MLEs <- optim(startValues, logTargetIPM, logTargetPars = IPMLTP2, returnNeg = T,
#               method = "SANN", control = list(maxit = 1000), printProp = T)
# saveRDS(MLEs, "../RDSobjects/HWP/Uniform/MLEsuniform")
# # Check some log lik values
# logTargetIPM(startValues, logTargetPars = IPMLTP, returnNeg = T, printProp = F)
# logTargetIPM(MLEs$par, logTargetPars = IPMLTP, returnNeg = T, printProp = F)
# logTargetIPM(SVjitter, logTargetPars = IPMLTP, returnNeg = T, printProp = F)

######################### RUN THE MCMC CHAINS IN PARALLEL ######################

# high resolution plots:
plotNames <- c("LL", "SurvI", "SurvG", "GrowI", "GrowG", "GrowS", "ReprI", "ReprG", "OffNum",
               "OffSurvI", "OffSurvG", "OffSurvS", "SurvChild", "ObsProbSh1","ObsProbSh2")

# Objects which the cluster needs to know about to run the sampling:
clNeeds = c('logTargetIPM', '%<>%', '%>%', 'sampleNorm', 'returnConstant', 'PoisNum',
            'evalPriors', 'sampleStateIPM_red', 'particleFilter','poissonObs',
            'multnomMu', "poissonMu", 'linLogit', 'vectorisedSamplerIPM', 'rmvnorm',
            'abind','colprods', 'detectionNumObs', 'SplitMNom','dmnom','multinomialObs',
            'countsToProbs','beta_mnomObs','beta_poisObs','linkNum','returnSelf','dbinom',
            'Schilder','beta_mnMu','beta_poisMu','fixedPoisObs','fixedMuObs')

# An example of how chains are run, using a stored proposal covariance structure
# and then stored for iteration:
print("And so it begins...")
ptm <- proc.time()
Sheepies <- pM_GaA(propCOV = propCOV,
                         lTarg = logTargetIPM, lTargPars = IPMLTP,
                         cores = ncores,
                         x0 = x0, # x0 = unlist(startValues),
                         itermax=itermax,
                         Params=list(GreedyStart=100,Pstar=0.234, gamzy0=1, epsilon=2, minVar=1e-6),
                         clNeeds = clNeeds, packages = "dissPackage3", prntPars=TRUE)

# Sheepies <- Nested_pMH(proposal = multvarNormProp, uFunc=NULL, #uFunc = multvarPropUpdate,
#                        propPars = diag(length(SVjitter))/100,
#                        #propPars = readRDS("ObsNumSigma4.1"),
#                        # propPars = propCov,
#                        lTarg = logTargetIPM, lTargPars = IPMLTP,
#                        cores = ncores,
#                        x0 = SVjitter, #x0 = unlist(startValues),
#                        itermax=itermax,
#                        clNeeds = clNeeds, packages = "dissPackage3", prntPars=TRUE)

# Sheepies <- MH(proposal = multvarNormProp, uFunc = NULL,
#                       propPars = diag(length(SVjitter))/100,
#                       #propPars = readRDS("ObsNumSigma4.1"),
#                       #propPars = propCov,
#                       lTarg = logTargetIPM, lTargPars = IPMLTP,
#                       x0 = SVjitter, itermax=itermax, prntPars = T)
ptm_fin<-proc.time() - ptm;
print(paste0("ncpus= ",ncores," : ",ptm_fin))

tag<-paste0(namer,"_",priorName,"_its",itermax,"_",
            gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = ""))
saveRDS(Sheepies, paste0(directory,"Results/",tag))


# 
# MLE<-readRDS("./Results/MLE_svjitter")
# Sheepies<-list(
#   readRDS("./Results/uniform_its18000_2020-04-18")[1:18001,],
#   readRDS("./Results/uniform_its18000_2020-05-08_565"),
#   readRDS("./Results/uniform_its18000_2020-05-10_565"),
#   readRDS("./Results/uniform_its18000_2020-05-11_565")
# )
# 
# 
# chosenChain<-Sheepies
# class(gro)<- 'pMCMCoutput'
# class(off)<- 'pMCMCoutput'
# chosenChain<-list(gro,off)
# 
# multvarNormProp(unlist(startValues),multvarPropUpdate(off))
# 
# getAcceptance(chosenChain)
# 
# # make plots for the first three parameters as figures for the writeup:
# chainFirst3 <- rep(NA, 3) %>% list
# chainFirst3[[1]] <- chosenChain[[1]][,4:6]
# class(chainFirst3) <- 'pMCMCoutput'
# F3plotNames <- plotNames[4:6]
# 
# plot(chainFirst3, cols = 1:3, width=8.3*3*(0.3), height=11.7*1.5*(0.3),
#      cex=1, names=F3plotNames, filePath="chainFirst3_")
# 
# # plot chain before thinning:
# plot(chosenChain[[1]], cols=2:Np, width=8.3*3, height=11.7*3,
#      cex=1, names=plotNames, filePath="")
# 
# # # plot chain after thinning:
# # thinnedChosen <- thinMCMC(chosenChain, alpha = 0.1)
# # plot(thinnedChosen, cols=2:Np, width=8.3*3, height=11.7*3,
# #      cex=2, names=plotNames, filePath="")
# # 
# # # sub-plot for the write-up:
# # plot(thinnedChosen, cols=c(2,4,7,10,13,14), width=8, height=8,
# #      cex=2, names=plotNames[c(1,2,4,7,10,13,14)], filePath="")
# 
# # effective sample sizes for the chains:
# multiESS((chosenChain[[1]][,4:6]))
# multiESS((chosenChain[[2]][,10:12]))
# # combine the thinned chain samples together:
# combinedChosen <- thinnedChosen[[1]]
# for (i in 2:length(chosenChain)) combinedChosen <- rbind(combinedChosen, thinnedChosen[[i]])
# 
# # apply link functions:
# returnSelf <- function(x) x
# links <- rep('returnSelf', 13)
# links[c(5,11)] <- 'exp'
# links[c(8,12,13)] <- 'plogis'
# for (i in 2:Np) combinedChosen[,i] <- match.fun(links[i-1])(combinedChosen[,i])
# 
# # get the MAP parameters:
# ind <- combinedChosen[,1] %>% which.max
# MAP <-  combinedChosen[ind, -1]
# 
# # get the median, mean and 95% CIs:
# means <- apply(combinedChosen[,-1], 2, mean)
# medis <- apply(combinedChosen[,-1], 2, median)
# lower <- apply(combinedChosen[,-1], 2, quantile, probs = 0.025)
# upper <- apply(combinedChosen[,-1], 2, quantile, probs = 0.975)
# 
# # make a data.frame, and then tables for the write-up:
# simulated <- c(-9.65, 3.77, 1.41, 0.56, 0.08, -7.23,
#                2.6, 1, 0.36, 0.71, 0.16, 0.873, 1)
# 
# tMLE<-MLE$par
# # SVjitter[5]<-exp(SVjitter[5]); SVjitter[11]<-exp(SVjitter[11]);
# # SVjitter[12]<-plogis(SVjitter[12]); SVjitter[8]<-plogis(SVjitter[8]); SVjitter[13]<-plogis(SVjitter[13]);
# tMLE[5]<-exp(tMLE[5]); tMLE[11]<-exp(tMLE[11]);
# tMLE[12]<-plogis(tMLE[12]); tMLE[8]<-plogis(tMLE[8]); tMLE[13]<-plogis(tMLE[13]);
# # means[5]<-exp(means[5]); means[11]<-exp(means[11]);
# # means[12]<-plogis(means[12]); means[8]<-plogis(means[8]); means[13]<-plogis(means[13]);
# # medis[5]<-exp(medis[5]); medis[11]<-exp(medis[11]);
# # medis[12]<-plogis(medis[12]); medis[8]<-plogis(medis[8]); medis[13]<-plogis(medis[13]);
# # MAP[5]<-exp(MAP[5]); MAP[11]<-exp(MAP[11]);
# # MAP[12]<-plogis(MAP[12]); MAP[8]<-plogis(MAP[8]); MAP[13]<-plogis(MAP[13]);
# # lower[5]<-exp(lower[5]); lower[11]<-exp(lower[11]);
# # lower[12]<-plogis(lower[12]); lower[8]<-plogis(lower[8]); lower[13]<-plogis(lower[13]);
# # upper[5]<-exp(upper[5]); upper[11]<-exp(upper[11]);
# # upper[12]<-plogis(upper[12]); upper[8]<-plogis(upper[8]); upper[13]<-plogis(upper[13]);
# 
# # summaries <- data.frame(simulated = simulated, SVjitter=SVjitter, MLE=tMLE, MAP = MAP, mean = means,
# #                         median = medis, lower = lower, upper = upper)
# summaries <- data.frame(simulated = simulated, MAP = MAP, mean = means,
#                         median = medis, lower = lower, upper = upper)
# 
# rownames(summaries) <- c("survival.i", "survival.g", "growth.i", "growth.g",
#                          "growth.sd", "repr.i", "repr.g", "litter.size",
#                          "off.size.i", "off.size.g", "off.size.sd",
#                          "child.survival.p", "observation.p")
# 
# xtable(t(summaries)[,c(1:5, 13)])
# xtable(t(summaries)[,-c(1:5, 13)])

# 
# # posterior growth rate distn from the samples:
# cluster <- makeCluster(4)
# CI <- posteriorGrowthRate(thinnedChosen,
#                           IPMLTP, growthFunc = "normal", printRate = 10,
#                           offSizeFunc = "normal", L = 1.5, U = 3.55, m = 500,
#                           cluster = cluster, clNeeds = clNeeds)
# CI