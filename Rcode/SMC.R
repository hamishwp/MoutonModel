library(dissPackage3)
library(mc2d)
library(abind)
# library(Rfast,include.only = "colprods")
####################### INITIAL DISTRIBUTION FUNCTIONS #########################

unifMu <- function(popSize, traitClasses, n=1){
  # purpose : returns a vector of individual counts which is uniformly 
  #           distributed across the dimensions of the state space
  # inputs  : popSize      - The integer total population size (to be spread
  #                          uniformly amongst the train classes) 
  #           traitClasses - The integer number of bins that the trait is split
  #                          into for the purpose of the analysis.
  return(rmultinom(n, popSize, rep(1/traitClasses, traitClasses)))
}

multnomMu <- function(popSize, probs, n=1, pobs){
  # purpose : returns a vector of individual counts which is uniformly 
  #           distributed across the dimensions of the state space
  # inputs  : popSize - The integer total population size (to be spread
  #                     uniformly amongst the train classes) 
  #           probs   - The probabilities of beingin each size class
  #           pobs   - probability of being observed
  
  return(rmultinom(n, round(popSize/pobs), probs))
  # return(rmultinom(n, popSize, probs))
}

poissonMu <- function(popSize, n=1, pobs){
  # purpose : returns a vector of 'real' population per size class
  #           calculated from the observed population and observed probability
  # inputs  : popSize - The integer population per size class
  #           pobs   - probability of being observed
  return(popSize+t(sapply((popSize*(1/pobs-1)),rpois,n=n)))
}

beta_mnMu <- function(popSize, probs, n=1, pobs){
  # purpose : returns a vector of individual counts which is uniformly 
  #           distributed across the dimensions of the state space
  # inputs  : popSize - The integer total population size (to be spread
  #                     uniformly amongst the train classes) 
  #           probs   - The probabilities of beingin each size class
  #           pobs   - probability of being observed
  p<-1/rbeta(n,pobs[1],pobs[2])
  return(sapply(p, function(p) rmultinom(1, round(popSize*p), probs) ))
}

beta_poisMu <- function(popSize, n=1, pobs){
  # purpose : returns a vector of 'real' population per size class
  #           calculated from the observed population and observed probability
  # inputs  : popSize - The integer population per size class
  #           pobs   - probability of being observed
  p<-rbeta(n,pobs[1],pobs[2])
  return(popSize+sapply(p,function(p) sapply((popSize*(1/p-1)),rpois,n=1)))
}

####################### IPM STATE SPACE SAMPLER ################################
sampleStateIPM <- function(previousState, survFunc, survPars,
                           growthSamp, growthPars, reprFunc, reprPars, 
                           offNumSamp, offNumPars, offSizeSamp, offSizePars,
                           Schild, breaks, oneSex = TRUE, checks = FALSE,
                           verbose = FALSE, sizes){
  # purpose : A barebones version of simulateIBM which takes as input a size 
  #           distribution at a given time, and given the parameters of the 
  #           vital rate functions, produces a sample of the size distribution
  #           of the individuals at the next time step (census).
  # inputs  : previousState - X: A vector of the number of individuals in each
  #                           size class. NOT OBSERVED BUT ACTUAL. Can be any length.
  #           survFunc      - The function which evaluates the probability of 
  #                           survival. The first argument is the size of the
  #                           individual.
  #           survPars      - The parameters for the survival function.
  #           growthSamp    - The function which samples the new size of
  #                           individuals, given their current size. The first
  #                           argument is the size.
  #           growthPars    - The parameters of the growth function.
  #           reprFunc      - The function which calculates the probability of
  #                           survival. The first argument is the size.
  #           reprPars      - The parameters for the reproduction function.
  #           offNumSamp    - The function which samples from the distribution
  #                           of the number of offspring. 
  #           offNumPars    - The parameters for the offspring number sampling
  #                           distribution.
  #           offSizeSamp   - The function which samples from the distribution
  #                           of offspring sizes given parent sizes. Parent size
  #                           is the first argument.
  #           offSizePars   - The parameters of the child size distribution.
  #           Schild        - The probability that a child survives to the 
  #                           first census.
  #           L             - The lower limit of the size distribution.
  #           U             - The upper limit of the size distribution.
  #           oneSex        - If TRUE, the survey is tracking only one sex, and
  #                           so we simulate the sex of the children with a 
  #                           Bernoulli(0.5) distribution before including them
  #                           in the survey.
  #           breaks        - The breakpoints of the size classes. Should be a 
  #                           vector of length D+1, where D is the number of
  #                           size classes. The extreme entries should be L and
  #                           U, but these will be replaced by -Inf and Inf
  #                           when producing the new size distribution, since
  #                           some non-truncated samplers for the growth and 
  #                           offspring size distributions will produce values
  #                           slightly outside of [L, U]
  #           checks        - If TRUE, performs some input checks
  #           verbose       - If TRUE prints out survival rate, average growth,
  #                           reproductive rate, birth rate, and child censusing
  #                           rate.
  # output  : A vector of counts of the same length as previousState
  # note : Could add functionality for customised breaks, maybe useful if 
  #        population is very clustered in a certain range.
  # note : All function arguments can be passed as the character name of the
  #        function instead
  
  # To ensure function names work too:
  survFunc %<>% match.fun ; growthSamp %<>% match.fun ; reprFunc %<>% match.fun
  offNumSamp %<>% match.fun ; offSizeSamp %<>% match.fun
  
  # Some input checks:
  if (checks){
    if (class(previousState)!='integer') stop('Invalid previous state')
    if (any(previousState%%1!=0)) stop('State space can contain only integers')
    if (!class(survFunc) %in% c('function','character')) stop('Invalid input')
    if (!class(growthSamp) %in% c('function','character')) stop('Invalid input')
    if (!class(reprFunc) %in% c('function','character')) stop('Invalid input')
    if (!class(offNumSamp) %in% c('function','character')) stop('Invalid input')
    if (!class(offSizeSamp) %in% c('function','character'))stop('Invalid input')
  }
  
  # Get midpoint of each user defined interval:
  D <- length(previousState)
  
  # Get the size distribution of individuals that survive:
  newSizes <- survFunc(sizes, survPars) %>% rbinom(n=D, size=previousState) %>%
    rep(x=sizes) %>% growthSamp(growthPars)
  if (length(newSizes)==0) return(rep(0, D))
  
  # Get the size distribution of newborns:
  reprProbs <- reprFunc(newSizes, reprPars)
  reprSizes <- weightedSelection(newSizes, reprProbs)
  if (length(reprSizes)==0) offSizes <- c()
  else{
    # Determine the sex of the offspring:
    probGend <- ifelse(isTRUE(oneSex), 0.5, 1)
    # Modify the number of reproducing sheep 
    reprSizes <- reprSizes[which(as.logical(rbinom(length(reprSizes), 1, probGend)))]
    
    offSizesTemp <- offNumSamp(reprSizes, offNumPars) %>% rep(x=reprSizes) %>%
      offSizeSamp(offSizePars)
    offSizes <- weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp)))
  }
  
  # Print out summaries of the simulation if desired:
  if (verbose){
    cat("survival rate is", length(newSizes)/sum(previousState),"\n")
    oldSizes <- rep(sizes, previousState) %>% mean
    cat("average growth is", mean(newSizes) - oldSizes, "\n")
    cat("reproduction rate is", length(reprSizes)/sum(previousState), "\n")
    cat("rate of births is", length(offSizesTemp)/sum(previousState), "\n")
    cat("child censusing rate is", length(offSizes)/length(offSizesTemp),"\n")
    cat("growth rate is", length(c(offSizes, newSizes))/sum(previousState),"\n")
    cat("\n")
  }
  
  # Get the new distribution of sizes:
  # breaks[c(1, D+1)] <- c(-Inf, Inf)
  vectorToCounts(c(offSizes, newSizes), breaks) %>% return  
  
}


####################### IPM STATE SPACE SAMPLER ################################
sampleStateIPM_red <- function(previousState, survFunc, survPars,
                               growthSamp, growthPars, reprFunc, reprPars, 
                               offNumSamp, offNumPars, offSizeSamp, offSizePars,
                               Schild, breaks, oneSex = TRUE, checks = FALSE,
                               verbose = FALSE, sizes){
  
  # Get the size distribution of individuals that survive:
  D<-length(previousState)
  # Generate the survival probability [0,1] for each 'size' bin
  # Using previousState (X - total #sheep per bin), predict number of OBSERVED sheep per bin
  # Generate a sample population at the bin midpoint (at 'sizes' values) but repeating the values
  newSizes <- rep(rbinom(n=D, size=previousState, prob = survFunc(sizes, survPars)),x=sizes)  
  if (length(newSizes)==0) return(rep(0, D))
  # What would the growth of such a population be in one year?
  newSizes <- growthSamp(newSizes, growthPars)
  
  # Get the probability that each class reproduces based on number of sheep (not observed but actual)
  # reprProbs <- reprFunc(newSizes, reprPars)
  reprProbs <- reprFunc(sizes, reprPars)
  stop("what is reprProbs meant to do?")
  # Get the updated size distribution of adults:
  reprSizes <- weightedSelection(newSizes, reprProbs)
  if (length(reprSizes)==0){ return(vectorToCounts(c(newSizes), breaks))
  }else{
    # Modify the number of reproducing sheep to retain females only
    if(oneSex) {
      reprSizes <- reprSizes[sample(1:length(reprSizes),ceiling(length(reprSizes)*0.5),replace = F)] # as.logical(rbinom(length(reprSizes), 1, 0.5))
      # if (length(reprSizes)==0) return(vectorToCounts(c(newSizes), breaks))
    }
    # Sample the size of the offspring
    offSizesTemp <- offSizeSamp(rep(reprSizes,offNumSamp(length(reprSizes), offNumPars)),offSizePars) 
    # Get the new distribution of sizes:
    offSizes <- weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp)))
    
  }
  
  return(vectorToCounts(c(offSizes, newSizes), breaks))
  
}

sampleStateIPM_ABCSIR <- function(previousState, survFunc, survPars,
                               growthSamp, growthPars, reprFunc, reprPars, 
                               offNumSamp, offNumPars, offSizeSamp, offSizePars,
                               Schild, breaks, oneSex = TRUE, checks = FALSE,
                               verbose = FALSE, sizes){
  # Initialise output data.frame
  outer<-data.frame()
  # Get the size distribution of individuals that survive:
  D<-length(previousState); reppie<-rep(0,D)
  # Generate the survival probability [0,1] for each 'size' bin
  # Using previousState (X - total #sheep per bin), predict number of OBSERVED sheep per bin
  # Generate a sample population at the bin midpoint (at 'sizes' values) but repeating the values
  newSizesI <- rep(rbinom(n=D, size=previousState, prob = survFunc(sizes, survPars)),x=sizes)  
  # If none survive, return empty set
  if (length(newSizesI)==0) return(
    array(c(reppie,reppie,reppie),
          dim = c(length(sizes),3),
          dimnames = list(round(sizes,digits = 2),
                          c("NoSurv","NoParents","NoOff"))))
  # What would the growth of such a population be in one year?
  newSizes <- growthSamp(newSizesI, growthPars)
  # Bin these by breaks
  newCounts<-vectorToCounts(newSizes, breaks)
  # Get the probability that each class reproduces based on number of sheep (not observed but actual)
  # stop("Something is wrong here, what is reprFunc meant to do and why did it have z=previousState")
  reprProbs <- reprFunc(sizes, reprPars)
  # Parent counts
  # stop("what is reprProbs meant to do? weightedSelection?")
  reprCounts <- rbinom(D, newCounts, reprProbs)
  # Check if any births occurred
  if (sum(reprCounts)==0) {
    return(array(c(newCounts,
                   reprCounts,
                   reppie),
          dim = c(length(sizes),3),
          dimnames = list(round(sizes,digits = 2),
                          c("NoSurv","NoParents","NoOff"))))
  } else {
    # Modify the number of reproducing sheep to retain females only
    if(oneSex) {
      # Are all sheep included in the study female?
      offCounts<-rbinom(D, reprCounts, 0.5)
    } else offCounts<-reprCounts
    # Number of offspring born per PARENT size class
    offCounts[offCounts!=0]<-offNumSamp(offCounts[offCounts!=0],offNumPars)
    # Save the total number born
    bornCount<-sum(offCounts)
    # How many survive to the next year?
    offCounts <- rbinom(D, offCounts, Schild)
    # Convert to parent sizes
    reprSizes <- rep(sizes,offCounts)
    # Sample the size of the offspring
    offSizes <- offSizeSamp(reprSizes,offSizePars) 
  }
  # Return the ABCSIR object with all the count data, for this year only. This will be combined with rbind later
  # Note that this data.frame has D rows (number of bins/breaks)
  # stop("Change reprFunc and weightedSelection in sampleSpaceIPM_red")
  # stop("Does this explain why simulated data comparisons never worked when true parameters were provided?")
  return(array(c(vectorToCounts(c(newSizes),breaks), # survived population based on previous census size
                 reprCounts, # parents
                 vectorToCounts(c(offSizes),breaks)), # offspring counts
               dim = c(length(sizes),3),
               dimnames = list(round(sizes,digits = 2),
                               c("NoSurv","NoParents","NoOff"))))
}

################################# UTILS ########################################
vectorisedSamplerIPM_ABCSIR <- function(initialStates, SamplerArgs){
  # purpose : Uses the sampleStateIPM function to project a n different state
  #           spaces forwards one step, given they all have the same parameters
  # inputs  : samplerArgs       - The list of NAMED arguments for the state
  #                               space sampler. Any argument which is a
  #                               function should be passed as a character 
  #                               instead. For example 'mean' instead of mean.
  #           initialStates     - The matrix containing the intital states.
  #                               Each column is a separate initial state. 
  # output  : A Dxt matrix, where D is the dimension of the state space and t
  #           is the number of time steps we project forwards for.
  
  helper <- function(vec) do.call(sampleStateIPM_ABCSIR, c(list(vec), SamplerArgs))
  
  array(apply(initialStates,2, helper,simplify = T),
        dim=c(length(SamplerArgs$sizes),3,ncol(initialStates)))
  
}

vectorisedSamplerIPM <- function(initialStates, SamplerArgs){
  # purpose : Uses the sampleStateIPM function to project a n different state
  #           spaces forwards one step, given they all have the same parameters
  # inputs  : samplerArgs       - The list of NAMED arguments for the state
  #                               space sampler. Any argument which is a
  #                               function should be passed as a character 
  #                               instead. For example 'mean' instead of mean.
  #           initialStates     - The matrix containing the intital states.
  #                               Each column is a separate initial state. 
  # output  : A Dxt matrix, where D is the dimension of the state space and t
  #           is the number of time steps we project forwards for.
  
  helper <- function(vec) do.call(sampleStateIPM_red, c(list(vec), SamplerArgs))
  return(apply(initialStates,2, helper))
}

# vectorisedSamplerIPM <- function(initialStates, SamplerArgs){
#   nmz<-paste0("p",1:SamplerArgs$D)
#   parsy<-str_flatten(paste0(nmz,"=initialStates[",1:SamplerArgs$D,",]"),collapse = ",")
#   tnm<-str_flatten(nmz,collapse = ",")
#   
#   eval(parse(text = paste0("tfun<-function(",tnm,"){do.call(sampleStateIPM_red,c(list(c(",tnm,")), SamplerArgs))}")))
#   VIPM<-Vectorize(tfun,vectorize.args = nmz)
#   return(eval(parse(text = paste0("VIPM(",parsy,")"))))
# }

projectStateSpace <- function(stateSpaceSampler, SamplerArgs, initialState, t,
                               ...){
  # purpose : Uses the sampleStateIPM function to project a current state space
  #           for multiple time steps. 
  # inputs  : stateSpaceSampler - The function which samples the state space at
  #                               the next time step, given the current state.
  #                               it is expected that the first argument of this
  #                               function is the vector representing the value
  #                               of the previous state
  #           samplerArgs       - The list of NAMED arguments for the state
  #                               space sampler. Any argument which is a
  #                               function should be passed as a character 
  #                               instead. For example 'mean' instead of mean.
  #           initialState      - The vector containing the intital state
  #           t                 - The number of time steps to simulate for
  # output  : A Dxt matrix, where D is the dimension of the state space and t
  #           is the number of time steps we project forwards for.
  
  D <- length(initialState)
  output <- matrix(NA, nrow=D, ncol=t+1)
  output[,1] <- initialState
  
  for (i in 2:(t+1)){
    tArgs <- c(list(output[,i-1]), SamplerArgs, ...)
    output[,i] <- do.call(stateSpaceSampler, tArgs)
  }
  
  return(output)
}

getSizeForCensus <- function(DF, censusNum, breaks){
  # purpose : Takes as input a data.frame of simulate data made by simulateIBM
  #           and returns the distribution of sizes at the given census number.
  # inputs  : DF        - The data.frame of simulated data
  #           censusNum - The  integer number of the census we want
  #           breaks    - The breaks for the size classes
  DF %<>% subset(!is.na(DF$size) & DF$census.number==censusNum)
  vectorToCounts(DF$size, breaks) %>% return
}

getSizeDistns <- function(DF, breaks){
  # purpose : Produce a matrix of oberved size distributions at each time step,
  #           given a data.frame of simulated data produced by the simulateIBM
  #           function.
  # inputs  : DF - The data.frame of simulated data
  #           s  - The number of size classes used for the discretisation.

  # In case user specified breaks are bad:
  # breaks %<>% sort
  # breaks[c(1, length(breaks))] <- c(-Inf, Inf)
  
  # Loop through the censuses and get the size distn:
  return(sapply(unique(DF$census.number), getSizeForCensus, DF=DF, breaks=breaks))
}

makeFakeTrueSizeDist <- function(probs, sizes){
  # purpose : creates a quick and dirty fake vector of true sizes,
  #           multinomially distributing observations across size classes
  # inputs  : probs - The probability of being in each size class
  #           sizes - The vector of the size of the population at each time 
  #                   step
  output <- matrix(NA, nrow=length(probs), ncol=length(sizes))
  for (i in 1:length(sizes)) output[,i] <- rmultinom(1, sizes[i], probs)
  return(output)
}
# Min-max scaling for the particle filter sample function
minmaxScale<-function(sw) (sw-min(sw,na.rm = T))/(max(sw,na.rm = T)-min(sw,na.rm = T))
eminmaxScale<-function(sw) exp((sw-min(sw,na.rm = T))/(max(sw,na.rm = T)-min(sw,na.rm = T)))/exp(1)
################### PARTICLE FILTER AND PMCMC IMPLEMENTATION ###################

# IPM Particle Filter function
particleFilter <- function(Sd, mu, muPar, sampleState, sampleStatePar, obsProb,
                           obsProbPar, NoParts, fixedObsProb=F){
  # purpose : Produces an estimate of the log likelihood of the model, given
  #           NoParts particles projected through the state space
  # inputs  : Sd             - Matrix of summary statistics of the observations.
  #                            Columns are time steps.
  #           mu             - The function which produces a sample from the 
  #                            initial distribution of the state space.
  #                          - multinomial with probability for each size class detection
  #           muPar          - The parameters for the function which samples
  #                            the initial distribution of the state space. It
  #                            is expected to have final argument 'n', the
  #                            number of samples we require (the argument name
  #                            can in fact be anything). Must be a list.
  #           sampleState    - The function which given the state space at time
  #                            t-1, produces a sample of the state space at time
  #                            t.
  #           sampleStatePar - The parameters for the function which samples
  #                            the state space.
  #           obsProb        - Function which evaluates the probability of the
  #                            observed out of true counts at time t, given the state at time t.
  #                          - binomial with observed probability p
  #           obsProbPars    - The parameters of the conditional likelihood of
  #                            the observations. Must be a list.
  #           NoParts        - The number of particles to produce.
  #           fixedObsProb   - Binary indicator of whether the observed probability 
  #                            is taken directly from the data (TRUE) or 
  #                            sampled from a beta distribution
  # output  : A matrix or array with n rows, t columns and the third dimension
  #           equal to the dimension of the state space.
  # In case of a single observation:
  t <- dim(Sd)[3]
  # Setup initial states
  prevStates <- do.call(mu, muPar)
  # Setup weight matrix and standardised weight matrix:
  output <- list(distance=0, sw=rep(1,NoParts),
                 shat=array(NA, dim=c(nrow(Sd),3,t)))
  # Update weights for first time step:
  wArgs <- list(Sd=Sd, pobs=obsProbPar, NoParts=NoParts)
  
  for (time in 1:t){
    wArgs$time<-time
    # Importance resample from the previous states 
    particleIndices <- sample(1L:NoParts, NoParts, replace = T, 
                              prob = output$sw)
    sampledStates <-  prevStates[, particleIndices]
    # IPM push forward
    wArgs$Sstar <- sampleState(sampledStates, sampleStatePar)
    # prevStates is the total population, by size bin
    prevStates<-wArgs$Sstar[,1,]+wArgs$Sstar[,3,]
    # Convert to observed from latent space & calculate the objective function vector
    output <- obsProb(output,wArgs)
  }
  output$sw<-NULL
  return(output)
  
}










