directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"

source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))
library(dissPackage3)
library(xtable)
library(mcmcse)

itermax <- 18000
# itermax <- 15
SMC_parts<-7000
# nchains<-4
# ncores<-detectCores()
# nchains<-4
ncores<-4
# if(ncores%%nchains!=0) {
#   print("Warning: #cores/#nchains not divisible, adjusting #chains=1")
#   nchains<-1
# }

########################### SPECIFYING THE PRIORS ##############################

priorName<-"uniform"

priorsIPM<-switch(priorName,uniform=rep("dunif", 14),cauchy=rep("dcauchy", 14),normal=rep("dnorm", 14),stop("Prior distribution not recognised"))
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
  offNum=list(min=0, max=3),
  # Childsize intercept:
  osFuncInt=listy,
  # Childsize gradient:
  osFuncGra=listy,
  # Childsize sd:
  osFuncSD=listy,
  # qlogis(Child survival):
  Schild=list(min=0, max=50),
  # qlogis(probability of detection of an animal):
  obsProb=list(min=1, max=50)
)

skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA,
  obsProbPar = NA
)

########## CALCULATE SOME VALUES BEFORE WE CAN MAKE THE IPMLTP OBJECT ##########

# simPars <- list(n=100, t=500,
#                 # set survival details:
#                 survFunc = linLogit, survPars = c(-9.65, 3.77),
#                 # set growth details:
#                 growthSamp = sampleDTN,
#                 growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
#                 # set reproduction details:
#                 reprFunc = linLogit, reprPars = c(-7.23, 2.6),
#                 # set offspring number and size distribution details:
#                 offNum=1, offSizeSamp = sampleDTN,
#                 offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
#                 # Child survival probability:
#                 Schild=0.873,
#                 # set other miscelaneous parameters:
#                 Start = 2.7, thresh=700, OneGend = TRUE,
#                 popPrint = F, verbose=F)

# alternative, narrower state-space parameters: 
#     survPars = c(-20.65,7.77)
#     reprPars = c(-11.23, 3.6)
#     growthPars = c(1.41, 0.56, log(0.04), 1.5, 3.55)

# Define the breaks for the 4 group specification:
breaks <- seq(1.5, 3.55, l=6)[-2]

# # Get the prior initial size distribution:
# priorLoopN <- 100
# priorCounts <- rep(0, length(breaks)-1)
# priorStartSize <- rep(NA, priorLoopN)
# 
# set.seed(201020)
# 
# for (i in 1:priorLoopN){
#   # Generate the data:
#   priorData <- do.call(simulateIBM, simPars)  
#   max.prior <- priorData$census.number %>% max
#   # Extract the counts in each size class for the first census:
#   loop.counts <- subset(priorData, priorData$census.number==(max.prior-5)) %>%
#     `$`(size) %>% na.omit %>% vectorToCounts(breaks = breaks)
#   # Update the previous vector:
#   priorCounts %<>% `+`(loop.counts)
#   priorStartSize[i] <- sum(loop.counts)
# }
# 
# # Get the parameters for the initial size distribution:
# priorProbs <- countsToProbs(priorCounts)
# priorStartSize <- (sum(priorCounts)/priorLoopN) %>% ceiling
# priorStartSizeSD <- sd(priorStartSize)
# ############################ GENERATE THE DATA #################################
# 
# # Get the observations:
# set.seed(102010)
# simmedData <- do.call(simulateIBM, simPars)
simmedData<-readRDS(paste0(directory,"RDSobjects/simmedData"))
max.cens <- simmedData$census.number %>% max
# Y is [size_distribution , year]
Y <- getSizeDistns(simmedData, breaks)[,(max.cens-6):max.cens]

#################### CREATE THE LOGTARGETPARAMETERS OBJECT #####################
# readRDS(paste0(directory,"RDSobjects/IPMLTP"))
priorProbs<-c(0.03374183, 0.16805116, 0.39873340, 0.39947362)

IPMLTP <- list(
  priorFunc = match.fun('evalPriors'),
  priors = priorsIPM,
  priorPars = flatPriors,
  skeleton = skeleton,
  survFunc = match.fun('linLogit'),
  growthSamp = match.fun('sampleNorm'),
  reprFunc = match.fun('linLogit'),
  offNumSamp = match.fun('returnConstant'),
  offSizeSamp = match.fun('sampleNorm'),
  oneSex = TRUE,
  mu = match.fun('multnomMu'),
  muPar = c(sum(Y[,1]), list(priorProbs)),
  b = SMC_parts, 
  Y = Y,
  obsProb = match.fun('detectionNumObs'),
  shift = qlogis(0.49),
  breaks = breaks
)

# saveRDS(IPMLTP, "IPMLTP")

########################### MAKE THE START VALUES ##############################

startValues <- list(
  survPars = c(-9.65, 3.77),
  growthPars = c(1.41, 0.56, log(0.08)),
  reprPars = c(-7.23, 2.6),
  offNumPars = 1,
  offSizePars = c(0.36, 0.71, log(0.16)),
  Schild = qlogis(0.873), 
  obsProbPar = 10 # not too close to 50 since this will hurt the chain
)

# startValues %<>% unlist
# SVjitter <- startValues + rnorm(length(startValues), 0, 0.5)
# saveRDS(SVjitter, paste0(directory,"RDSobjects/HWP/startvals_jitter"))
SVjitter<-readRDS(paste0(directory,"RDSobjects/startvals_jitter"))

print("Initial Values Log-Likelihood=")
# print(logTargetIPM(SVjitter, logTargetPars = IPMLTP, returnNeg = T, printProp = F))
print(1563.144)
print(" ")

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

################################################################################

# SANN<-readRDS("../RDSobjects/HWP/Normal/SANNnormal")
# MLEs<-readRDS("../RDSobjects/HWP/Normal/MLEsnormal")

######################### RUN THE MCMC CHAINS IN PARALLEL ######################

# high resolution plots:
plotNames <- c("log(pi)", "Si", "Sg", "Gi", "Gg", "log(Gs)", "Ri", "Rg", "ON",
               "OSi", "OSg", "log(OSs)", "qlogis(Schild)", "qlogis(ObsP)")

# Objects which the cluster needs to know about to run the sampling:
clNeeds = c('logTargetIPM', '%<>%', '%>%', 'sampleNorm', 'returnConstant',
            'detectionNumObs', 'evalPriors', 'sampleStateIPM_red', 'particleFilter',
            'multnomMu', 'linLogit', 'vectorisedSamplerIPM', 'rmvnorm') 


# proposal <- multvarNormProp; uFunc <- NULL;
# propPars <- diag(length(SVjitter))/100;
# #propPars <- readRDS("ObsNumSigma4.1");
# #propPars <- propCov;
# lTarg <- logTargetIPM; lTargPars <- IPMLTP;
# cores <- 4; nChains <- 4;
# x0 <- SVjitter; itermax<-itermax; prntPars <- T;
# clNeeds <- clNeeds; packages <- "dissPackage3"


# An example of how chains are run, using a stored proposal covariance structure
# and then stored for iteration:
ptm <- proc.time()
Sheepies <- Nested_pMH(proposal = multvarNormProp, uFunc = NULL,
                         propPars = diag(length(SVjitter))/100,
                         #propPars = readRDS("ObsNumSigma4.1"),
                         #propPars = propCov,
                         lTarg = logTargetIPM, lTargPars = IPMLTP,
                         cores = 4, x0 = SVjitter, itermax=itermax,
                         clNeeds = clNeeds, packages = "dissPackage3", prntPars=TRUE)
ptm_fin<-proc.time() - ptm;
print(ptm_fin)
# user  system elapsed 
# 0.930   0.743 508.646 
print(" ")

# Sheepies <- MH(proposal = multvarNormProp, uFunc = NULL,
#                 propPars = diag(length(SVjitter))/100,
#                 #propPars = readRDS("ObsNumSigma4.1"),
#                 #propPars = propCov,
#                 lTarg = logTargetIPM, lTargPars = IPMLTP,
#                 x0 = SVjitter, itermax=itermax, prntPars = T)

tag<-paste0(priorName,"_its",itermax,"_",Sys.Date())
saveRDS(Sheepies, paste0(directory,"RDSobjects/",tag))

