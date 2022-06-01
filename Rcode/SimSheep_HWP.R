# directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
directory<-paste0(getwd(),"/")
setwd(directory)

list.of.packages <- c("mcmcse","xtable","magrittr","doParallel","Rfast","mc2d", "abind")
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
max.cens <- simmedData$census.number %>% max
simmedData%<>%filter(census.number>=max(c((max.cens-12L),1L)))
# Y is [size_distribution , year]
Y <- getSizeDistns(simmedData, breaks)[,(max.cens-12):max.cens]

priorProbs<-rowSums(Y)/sum(Y)
# priorProbs<-c(0.03374183, 0.16805116, 0.39873340, 0.39947362)
########################### MAKE THE START VALUES ##############################

x0<-getInitialValues(simmedData,printPars = T)
stop("Initial Values don't correspond to realdataset")
vals
unlist(x0)
sqrt(sum((vals-unlist(x0))*(vals-unlist(x0))))

itermax <- 5000
SMC_parts<-6700
ncores<-8

# x0[["obsProbPar"]]<-qlogis(0.95)
# x0[["offNumPars"]]<-1
# Import x0 value:
# x0<-readRDS(paste0(directory,"./Results/x0_MH_GSF_poissonMu_MultinomialObs"))
# x0[["offNumPars"]]<-1
# x0<-0.99*unlist(startValues)

namer<-"zeroReprod_sims_GSF_poissonMu_multinomObs_GLMx0_Retake"

saveRDS(simmedData, paste0("simmedData_",namer))

# propCov<-readRDS(paste0(directory,"gro")) #+ diag(13)*1/100
mufun<-match.fun('poissonMu') #match.fun('multnomMu')
muPar<-list(popSize=Y[,1],n=SMC_parts) #list(popSize=sum(Y[,1]), n=SMC_parts, probs=priorProbs) 

obsfun<-match.fun('multinomialObs') # match.fun('poissonObs') # match.fun('detectionNumObs')

# diagCOV<-x0/x0/5e5 (# Easy way to keep the variable names)
# diagCOV[["offNumPars"]]<-1e-9
# propCOV<-diag(diagCOV)
# Import covariance matrix:
# propCOV<-readRDS(paste0(directory,"/Results/propCOV_MH_GSF_poissonMu_MultinomialObs"))
propCOV<-diag(length(x0))/100

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
  Schild=listy,
  # qlogis(probability of detection of an animal):
  obsProb=listy
)

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
  'Schilder', # Offspring Survival Probability
  'plogis' # Observed Probability
)

########################### AUTOCALCULATE START VALUES ##############################


########################### CALCULATE EXPECTATION VALUES ##############################
# for (i in 1:length(flatPriors)){
#   
#   
#   
# }

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

#################### CREATE THE LOGTARGETPARAMETERS OBJECT #####################
# readRDS(paste0(directory,"RDSobjects/IPMLTP"))

IPMLTP <- list(
  priorFunc = match.fun('evalPriors'),
  priors = priorsIPM,
  priorPars = flatPriors,
  skeleton = skeleton,
  survFunc = match.fun('linLogit'),
  growthSamp = match.fun('sampleNorm'),
  reprFunc = match.fun('linLogit'),
  offNumSamp = match.fun('PoisNum'),
  offSizeSamp = match.fun('sampleNorm'),
  oneSex = oneSex,
  mu = mufun, # mu = match.fun('multnomMu'), 
  muPar = muPar, # muPar = list(popSize=sum(Y[,1]), probs=priorProbs),
  b = SMC_parts, 
  Y = Y,
  obsProb = obsfun, # obsProb = match.fun('detectionNumObs'), # obsProb = match.fun('poissonObs'),
  fixedObsProb=F,
  breaks = breaks,
  sizes=sizes,
  links = links
)
# obsProb<-match.fun("multinomialObs1D")
# saveRDS(IPMLTP, "IPMLTP")

startValues <- unlist(x0)
# SVjitter <- startValues + rnorm(length(startValues), 0, 0.5)
# saveRDS(SVjitter, paste0(directory,"RDSobjects/HWP/startvals_jitter"))
# SVjitter<-readRDS(paste0(directory,"RDSobjects/startvals_jitter"))

# print("Initial Values Log-Likelihood=")
# print(logTargetIPM(SVjitter, logTargetPars = IPMLTP, returnNeg = T, printProp = F))
# print(1563.144)
# print(" ")

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
plotNames <- c("LL", "SurvI", "SurvG", "GrowI", "GrowG", "GrowS", "ReprI", "ReprG", "OffNum",
               "OffSurvI", "OffSurvG", "OffSurvS", "SurvChild", "ObsProb")

# Objects which the cluster needs to know about to run the sampling:
clNeeds = c('logTargetIPM', '%<>%', '%>%', 'sampleNorm', 'returnConstant', 'PoisNum',
            'evalPriors', 'sampleStateIPM_red', 'particleFilter','poissonObs',
            'multnomMu', "poissonMu", 'linLogit', 'vectorisedSamplerIPM', 'rmvnorm',
            'abind','colprods', 'detectionNumObs', 'SplitMNom','dmnom','multinomialObs',
            'countsToProbs','beta_mnomObs','beta_poisObs','linkNum','returnSelf',
            'Schilder','returnConstant','beta_mnMu','beta_poisMu','fixedPoisObs',
            'fixedMuObs') 


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
print("And so it begins...")
ptm <- proc.time()
Sheepies <- pM_GaA(propCOV = propCOV,
                         lTarg = logTargetIPM, lTargPars = IPMLTP,
                         cores = ncores,
                         x0 = x0, # x0 = unlist(startValues),
                         itermax=itermax,
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
ptm_fin<-proc.time() - ptm;
print(paste0("ncpus= ",ncores," : ",ptm_fin))

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
# plot(chosenChain[[1]], cols=2:14, width=8.3*3, height=11.7*3,
#      cex=1, names=plotNames, filePath="")
# 
# # # plot chain after thinning:
# # thinnedChosen <- thinMCMC(chosenChain, alpha = 0.1)
# # plot(thinnedChosen, cols=2:14, width=8.3*3, height=11.7*3,
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