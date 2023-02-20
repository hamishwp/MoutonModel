# Housekeeping
IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- NULL; x0%<>%unlist(); funcys<-NULL

# Temporarily set the number of mSMC particles at 500:
SMC_parts<-500
# Particle filter initialisation parameters
if(muModel=='multinomial'){
  # Multinomial Mu initialisation function
  muPar<-list(popSize=rowSums(lSHEEP$COUNTS[,1]), n=SMC_parts, probs=lSHEEP$priorProbs)
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
  obsProbPar<-lSHEEP$obsProbTime
  muPar$pobs<-obsProbPar[1]
  IPMLTP %<>% c(list(obsProbPar = obsProbPar))
} else {
  obsProbPar <- proposed$obsProbPar
  muPar$pobs<-obsProbPar
}
# For the integral component we need infinite bounds
lSHEEP$breaks[c(1,nbks)]<-c(-Inf,Inf)
# Add to the log target parameters
IPMLTP %<>% c(list(oneSex = oneSex,
                   mu = mufun, # mu = match.fun('multnomMu'), 
                   Y = lSHEEP$COUNTS,
                   SumStats=lSHEEP$SumStats,
                   obsProb = obsfun, # obsProb = match.fun('detectionNumObs'), # obsProb = match.fun('poissonObs'),
                   fixedObsProb=fixedObsProb,
                   breaks = lSHEEP$breaks,
                   sizes = lSHEEP$sizes))

##################### ESTIMATE REQUIRED SMC PARTICLES ########################
# First check if this specific model has already been run
if(!calcParts & file.exists(paste0(directory,"Results/SMCPARTS_calc.RData"))){SMC_parts<-readRDS(paste0(directory,"Results/SMCPARTS_calc.RData"))$SMC_parts
} else print("Please run the code to estimate the lower bound of number of model-SMC particles required")
  
# Particle initialisation of the population size distributions
if(muModel=='multinomial'){
  # Multinomial Mu initialisation function
  muPar<-list(popSize=sum(lSHEEP$COUNTS[,1]), n=SMC_parts, probs=lSHEEP$priorProbs)
} else if(muModel=='poisson'){
  # Poisson Mu initialisation function
  muPar<-list(popSize=lSHEEP$COUNTS[,1],n=SMC_parts)
} else stop("Error in Mu Model input, we recommend muModel<-'poissonMu' ")
# Make sure to add the observation probability parameters
if(fixedObsProb){
  muPar$pobs<-obsProbPar[1]
} else {
  muPar$pobs<-obsProbPar
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

# Proposal Distribution
PropN<-do.call(getInitialValDists,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb,invlinks=IPMLTP$invlinks,MultiSD=1)))
# Either multivariate skew normal or multivariate normal
if(PropDist=="MVSN") {ProposalDist<-function(initSIR) PropN$proposal(ABCNP*ABCk)
}else ProposalDist<-function(initSIR) multvarNormProp(xt=PropN$x0, propPars=PropN$propCOV, n=ABCNP*ABCk)
# Setup the initial values for the ABSSIR algorithm:
initSIR<-list(ProposalDist=ProposalDist, 
              x0=PropN$x0, propCOV=PropN$propCOV,
              itermax=itermax, stepmax=stepmax,
              Np=ABCNP,
              k=ABCk) 

#################### CREATE THE LOGTARGETPARAMETERS OBJECT #####################
# Storage object for the ABC-SIR algorithm
# stop("Sort this out!")

IPMLTP %<>% c(list(
  priorFunc = match.fun('evalPriors'),
  priors = priorsIPM,
  priorPars = flatPriors,
  cores = ncores,
  DTN = data.frame(L=lSHEEP$L,U=lSHEEP$U),
  muPar = muPar, 
  b = SMC_parts,
  solveDF=lSHEEP$solveDF
))

IPMLTP$priorF<-function(thth) IPMLTP$priorFunc(thth, IPMLTP$priors, IPMLTP$priorPars)

if(NpCheck){
  print("Checking that number of samples per ABC-step is more than the minimum for the adaptive ABC threshold algorithm")
  simil<-data.frame()
  for(partsy in c(30,50,100,300,500,1000)){
    sss<-c()
    for(i in 1:150){
      xNew<-multvarNormProp(xt=x0, propPars=propCOV, n=partsy)
      xPrev<-multvarNormProp(xt=x0, propPars=propCOV, n=partsy)
      sss%<>%c(Supremum(0.95,xPrev,xNew,warny = F)[1])
    }
    simil%<>%rbind(data.frame(Np=partsy,similarity=sss))
  }
  p<-ggplot(simil,aes(factor(Np),similarity,group=factor(Np)))+
    geom_violin(aes(fill=factor(Np)),scale = "width")+
    labs(fill="ABC Samples")+xlab("Number of ABC Samples") + ylab("KLIEP Ratio (True Value = 1)");p
  ggsave("./Np-Supremum_similarity.png",p)
  saveRDS(simil,"Np-Supremum_similarity.RData")
} 




