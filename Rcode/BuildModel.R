# Housekeeping
IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- NULL; x0%<>%unlist(); lSHEEP$solveDF<-NULL

if(algorithm=="ABCSIR"){
  if(obsModel=='poisson'){
    funcys<-list(
      NoSurv=function(true,sim,wArgs) {dpois(wArgs$pobs*sim+1L, lambda = true+1L,log = T)},
      NoAlive=function(true,sim,wArgs) {dpois(wArgs$pobs*sim+1L, lambda = true+1L,log = T)},
      NoParents=function(true,sim,wArgs) {dpois(wArgs$pobs*sim+1L, lambda = true+1L,log = T)},
      GrowCounts=function(true,sim,wArgs) {dpois(abs(wArgs$pobs*sim-true)+1L, lambda = 1,log = T)},
      avSurvOff=function(true,sim,wArgs) {rep(dbinom(sim[1],sim[2],prob = true[1]/true[2],log = T),D)},
      NoOff=function(true,sim,wArgs) {dpois(wArgs$pobs*sim+1L, lambda = true+1L,log = T)}
    )
  } else if(obsModel=='binomial'){
    funcys<-list(
      NoSurv=function(true,sim,wArgs) {dbinom(sim, true, wArgs$pobs, log = T)},
      NoAlive=function(true,sim,wArgs) {dbinom(sim, true, wArgs$pobs, log = T)},
      NoParents=function(true,sim,wArgs) {dbinom(sim, true, wArgs$pobs, log = T)},
      GrowCounts=function(true,sim,wArgs) {dpois(abs(wArgs$pobs*sim-true)+1L, lambda = 1,log = T)},
      avSurvOff=function(true,sim,wArgs) {rep(dbinom(sim[1],sim[2],prob = true[1]/true[2],log = T),D)},
      NoOff=function(true,sim,wArgs) {dbinom(sim, true, wArgs$pobs, log = T)}
    )
  }
  if(fixedObsProb){
    
    obsfun<-function(true,sim,funcys,wArgs){
      # first modify the true population from the observed to the latent/true values
      # column-wise calculations of objective function elements (is a matrix)
      vapply(1:length(funcys),
             FUN = function(i) funcys[[i]](true[,i,rep(wArgs$time,wArgs$NoParts)],sim[,i,],wArgs),
             FUN.VALUE = numeric(nrow(true)))
    }
  } else {
    
    obsfun<-function(true,sim,funcys,wArgs){
      
      
      # Simulate the obsProbPar from the beta distribution here
      wArgs$pobs<-rbeta(ncol(true),wArgs$p[1],wArgs$p[2])
      # first modify the true population from the observed to the latent/true values
      # column-wise calculations of objective function elements (is a matrix)
      vapply(1:ncol(true),
             FUN = function(i) funcys[[i]](true[,i],sim[,i],wArgs),
             FUN.VALUE = numeric(nrow(true)))
    }
  }
  
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
    
  } else stop("Error in Obs Model input ")
  
  
  
  
  
  
  
  
  
  
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
} else {
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
stop("Sort this out!")
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




