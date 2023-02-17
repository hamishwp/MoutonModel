# Housekeeping
IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- NULL; x0%<>%unlist(); funcys<-NULL

# Combining the distance function with the observation model - 
if(obsModel=="ProbMink"){
  # multinomial observation probability-based
  if(fixedObsProb){
    obsfun<-function(output,wArgs){
      # Calculate the distances per particle
      MultiMod<-function(i) multinomialObs(wArgs$Sd[,i,rep(wArgs$time,dim(wArgs$Sstar)[3])], 
                                           wArgs$Sstar[,i,], 
                                           wArgs$pobs[wArgs$time], 
                                           logy=T)
      # For each type of summary statistic
      output$sw<-rowSums(sapply(1:3,MultiMod)); output$sw[is.na(output$sw)]<-log(.Machine$double.xmin)
      # Calculate the total distances
      output$d<-output$d+sum(output$sw,na.rm = T)
      # Exponentiate and scale the particle weights to become from 0 to 1
      output$sw<-exp(output$sw-max(output$sw))
      # Census-dependent dataframe
      output$shat[,,wArgs$time]<-apply(wArgs$Sstar,1:2,median,na.rm=T)
      
      return(output)
    }
  } else {
    obsfun<-function(output,wArgs){
      # Simulate the obsProbPar from the beta distribution here
      pobs<-rbeta(dim(wArgs$Sstar)[3],wArgs$pobs[1],wArgs$pobs[2])
      # Calculate the distances per particle
      MultiMod<-function(i) multinomialObs(wArgs$Sd[,i,rep(wArgs$time,dim(wArgs$Sstar)[3])], 
                                           wArgs$Sstar[,i,], 
                                           pobs, 
                                           logy=T)
      # For each type of summary statistic
      output$sw<-rowSums(sapply(1:3,MultiMod)); output$sw[is.na(output$sw)]<-log(.Machine$double.xmin)
      # Calculate the total distances
      output$d<-output$d+sum(output$sw,na.rm = T)
      # Exponentiate and scale the particle weights to become from 0 to 1
      output$sw<-exp(output$sw-max(output$sw))
      # Census-dependent dataframe
      output$shat[,,wArgs$time]<-apply(wArgs$Sstar,1:2,median,na.rm=T)
      
      return(output)
    }
  }
} else {
  if(fixedObsProb){
    # The observed probability is extracted directly from the data
    obsfun<-function(output,wArgs){
      # first modify the true population from the observed to the latent/true values
      Sstar<-apply(wArgs$Sstar*wArgs$pobs[wArgs$time],3,rbind)
      # Calculate the distances
      disties<-Minkowski(sest = Sstar,
                         sobs = c(wArgs$Sd[,,wArgs$time]),
                         dimmie=dim(output$shat)[3])
      # Merge into a single output object
      output$d<-disties$d + output$d
      output$sw<-disties$sw
      # Census-dependent dataframe
      output$shat[,,wArgs$time]<-array(disties$shat,dim(output$shat)[1:2])
      return(output)
    }
  } else {
    # The observed probability is sampled, per particle, from a beta distribution
    obsfun<-function(output,wArgs){
      # Simulate the obsProbPar from the beta distribution here
      pobs<-rbeta(dim(wArgs$Sstar)[3],wArgs$pobs[1],wArgs$pobs[2])
      # first modify the true population from the observed to the latent/true values
      Sstar<-apply(wArgs$Sstar*pobs,3,rbind)
      # Calculate the distances
      disties<-Minkowski(sest = Sstar, 
                         sobs = c(wArgs$Sd[,,wArgs$time]),
                         dimmie=dim(output$shat)[3])
      # Merge into a single output object
      output$d<-disties$d + output$d
      output$sw<-disties$sw
      # Census-dependent dataframe
      output$shat[,,wArgs$time]<-array(disties$shat,dim(output$shat)[1:2])
      return(output)
    }
  }
}

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
} else {
  print("Recalculating the number of particles required in the particle filter of the model, this will take some time...")
  # Run particleFilter with varying NoParts
  # First prepare everything for sampling, including initial values
  proposed<-Sample2Physical(x0,IPMLTP)
  stateSpaceSampArgs <- list(survFunc = IPMLTP$survFunc, survPars = proposed$survPars,
                             growthSamp = IPMLTP$growthSamp, growthPars = proposed$growthPars,
                             reprFunc = IPMLTP$reprFunc, reprPars = proposed$reprPars, 
                             offNumSamp = IPMLTP$offNumSamp, offNumPars = proposed$offNumPars,
                             offSizeSamp = IPMLTP$offSizeSamp, breaks = IPMLTP$breaks,
                             offSizePars = proposed$offSizePars, Schild=proposed$Schild,
                             sizes=IPMLTP$sizes, oneSex = IPMLTP$oneSex)
  
  # Now we calculate the standard deviation of the samples, when using specific numbers of particles - NoParts
  tmp<-c(unlist(mclapply(1:2000, function(i) {particleFilter(Sd=lSHEEP$SumStats, mu=mufun, muPar=muPar, obsProb = obsfun,
                                                    sampleState = vectorisedSamplerIPM_ABCSIR,
                                                    sampleStatePar = stateSpaceSampArgs,
                                                    obsProbPar = obsProbPar, 
                                                    fixedObsProb=fixedObsProb,
                                                    NoParts = 500)$d},mc.cores = ncores)))
  # We can see how many samples we need to achieve a stable standard deviation by plotting the cumulative SD (around 1800)
  plot(vapply(seq_along(tmp), function(i) sd(tmp[1:i]), 1),xlab = "Number of Samples",ylab = "Cumulative S.D.",main = "How many samples required for S.D. ~ 600")
  numsam<-2000
  
  parties<-(1:30)*500
  coster<-sapply(parties,function(pp) {
    muPar$n<-pp
    c(unlist(mclapply(1:numsam, function(i) {particleFilter(Sd=lSHEEP$SumStats, mu=mufun, muPar=muPar, obsProb = obsfun,
                      sampleState = vectorisedSamplerIPM_ABCSIR,
                      sampleStatePar = stateSpaceSampArgs,
                      obsProbPar = obsProbPar, 
                      fixedObsProb=fixedObsProb,
                      NoParts = pp)$d},mc.cores=ncores)))
    })
  # Calculate the running (k=3) standard deviation of the distance by number of particles
  runSD<-data.frame(parties=parties[3:length(parties)],runSD=vapply(3:length(parties), function(i) sd(log(-apply(coster,2,mean))[(i-3):i]), 1))
  # To find the ideal number of particles, fit various linear models trained with 
  # the data using sequentially less of the lower no. particle values
  # and predict when the gradient of the linear model equals zero.
  # Note that this method most likely over-estimates the required number of particles.
  tmp<-data.frame(parties=parties[1:(nrow(runSD)-1)],grad=vapply(1:(nrow(runSD)-1),function(i) (lm(runSD~parties,runSD[i:nrow(runSD),]))$coefficients[2],1))
  # Predict the required number of mSMC particles!
  SMC_parts<-round(unname(predict(lm(parties~poly(grad,4),tmp),newdata = data.frame(grad=c(0)))))

  p<-runSD%>%ggplot(aes(parties,runSD))+geom_point()+geom_vline(xintercept = SMC_parts)+
    xlab("Number of Particles in mSMC")+ylab("Running S.D (k=3)")+
    ggtitle("Convergence in No. Particles by S.D")+theme(plot.title = element_text(hjust = 0.5))
  q<-data.frame(meany=apply(coster,2,mean),parties=parties)%>%ggplot(aes(parties,meany))+geom_point()+geom_vline(xintercept = SMC_parts)+
    xlab("Number of Particles in mSMC")+ylab("Average Distance")+
    ggtitle("Convergence in No. Particles w.r.t Distance")+theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(directory,"Plots/Hamish/Convergence/NoPart_mSMC.png"),
         gridExtra::grid.arrange(q,p,nrow=1),width = 10,height = 4); gridExtra::grid.arrange(q,p,nrow=1)
  
  saveRDS(list(SMC_parts=SMC_parts,parties=parties,coster=coster,runSD=runSD),paste0(directory,"Results/SMCPARTS_calc.RData"))
}

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




