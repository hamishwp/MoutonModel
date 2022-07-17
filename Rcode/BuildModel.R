# Housekeeping
IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- NULL; x0%<>%unlist(); lSHEEP$solveDF<-NULL

# Selecting dependent upon configuration
if(obsModel=='multinomial'){
  # Multinomial observational model
  if(fixedObsProb) obsfun<-match.fun('fixedMuObs')
  else obsfun<-match.fun('beta_mnomObs')
  # Number of SMC particles
  SMC_parts<-7000  # 6700 # 200
  
} else if(obsModel=='poisson'){
  # Poisson observational model
  if(fixedObsProb) obsfun<-match.fun('fixedPoisObs')
  else obsfun<-match.fun('beta_poisObs')
  # Number of SMC particles
  SMC_parts<-600  # 6700 # 200
  
} else if(obsModel=='binomial'){
  # Binomial observational model
  if(fixedObsProb) obsfun<-match.fun('fixedBinObs')
  else obsfun<-match.fun('beta_binObs')
  # Number of SMC particles
  SMC_parts<-1750  # 6700 # 200
  
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
lSHEEP$breaks[c(1,nbks)]<-c(-Inf,Inf)

# Storage object for the ABC-SIR algorithm
outshell<-data.frame(matrix(nrow = 0,ncol = (1L+length(c(lSHEEP$COUNTS))+1L+length(x0))))
names(outshell)<-c("d",lSHEEP$cNames,"Accepted",names(unlist(IPMLTP$skeleton)))

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
  sizes = lSHEEP$sizes,
  outshell = outshell,
  cNames=lSHEEP$cNames,
  pNames=names(unlist(IPMLTP$skeleton)),
  cores=ncores
))
if(fixedObsProb) IPMLTP %<>% c(list(obsProbPar = obsProbTime))
#if(normsampler=="sampleDTN") IPMLTP %<>% c(list(DTN = c(lSHEEP$L,lSHEEP$U)))
IPMLTP %<>% c(list(DTN = data.frame(L=lSHEEP$L,U=lSHEEP$U)))



