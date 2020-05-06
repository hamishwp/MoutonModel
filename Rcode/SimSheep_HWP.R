directory<-"/home/patten/Documents/Coding/MoutonModel/"

list.of.packages <- c("mcmcse","xtable","mvtnorm","magrittr","doParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# if(length(list.of.packages[!("mvtnorm" %in% installed.packages()[,"Package"])])){devtools::install_github("rCarto/osrm")}

source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))
library(dissPackage3)
library(xtable)
library(mcmcse)

itermax <- 18000
#itermax <- 33
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
                         cores = ncores, x0 = SVjitter, itermax=itermax,
                         clNeeds = clNeeds, packages = "dissPackage3", prntPars=TRUE)
ptm_fin<-proc.time() - ptm;
print(paste0("ncpus= 7 : ",ptm_fin))

# SANN <-  optim(SVjitter, logTargetIPM, logTargetPars = IPMLTP, returnNeg = T,
#                method = "SANN", control = list(maxit = 1000), printProp = T,
#                hessian = T)
# 
# Sheepies <- MH(proposal = multvarNormProp, uFunc = NULL,
#                 propPars = diag(length(SVjitter))/100,
#                 #propPars = readRDS("ObsNumSigma4.1"),
#                 #propPars = propCov,
#                 lTarg = logTargetIPM, lTargPars = IPMLTP,
#                 x0 = SVjitter, itermax=itermax, prntPars = T)

tag<-paste0(priorName,"_its",itermax,"_",Sys.Date(),"_",ceiling(runif(1)*1000))
saveRDS(Sheepies, paste0(directory,"Results/",tag))

# MLE<-readRDS("./Results/MLE_svjitter")
# Sheepies<-readRDS("./Results/uniform_its18000_2020-04-18")
# 
# chosenChain<-list(Sheepies)
# 
# 
# 
# 
# 
# getAcceptance(chosenChain)
# 
# # make plots for the first three parameters as figures for the writeup:
# chainFirst3 <- rep(NA, 3) %>% list
# chainFirst3[[1]] <- chosenChain[[1]][,2:4] 
# class(chainFirst3) <- 'pMCMCoutput'
# F3plotNames <- plotNames[1:3]
# 
# plot(chainFirst3, cols = 2:3, width=8.3*3*(0.3), height=11.7*1.5*(0.3),
#      cex=1, names=F3plotNames, filePath="chainFirst3_")
# 
# # plot chain before thinning:
# plot(chosenChain[[1]], cols=2:14, width=8.3*3, height=11.7*3,
#      cex=2, names=plotNames, filePath="")
# 
# # plot chain after thinning:
# thinnedChosen <- thinMCMC(chosenChain, alpha = 0.1)
# plot(thinnedChosen, cols=2:14, width=8.3*3, height=11.7*3,
#      cex=2, names=plotNames, filePath="")
# 
# # sub-plot for the write-up:
# plot(thinnedChosen, cols=c(2,4,7,10,13,14), width=8, height=8,
#      cex=2, names=plotNames[c(1,2,4,7,10,13,14)], filePath="")
# 
# # effective sample sizes for the chains:
# multiESS(thinnedChosen[[1]][,-1])
# # combine the thinned chain samples together:
# combinedChosen <- thinnedChosen[[1]]
# for (i in 2:4) combinedChosen <- rbind(combinedChosen, thinnedChosen[[i]])
# 
# # apply link functions:
# returnSelf <- function(x) x
# links <- rep('returnSelf', 13)
# links[c(5,11)] <- 'exp'
# links[12:13] <- 'plogis'
# for (i in 2:14) combinedChosen[,i] <- match.fun(links[i-1])(combinedChosen[,i])
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
# MLE<-MLE$par
# 
# summaries <- data.frame(simulated = simulated, MLE=MLE, MAP = MAP, mean = means,
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