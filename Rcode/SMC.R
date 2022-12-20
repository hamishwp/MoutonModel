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

################### OBSERVATION CONDITIONAL LIKELIHOODS ########################

binomialObs <- function(Y, X, p){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability, so that the number of observed
  #           animals in each size class is binomially distributed with number
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # For a single time step:
  helper <- function(yx, p){
    n <- length(yx)
    y <- yx[1:n]
    x <- yx[-(1:n)]
    dbinom(y, x, p) %>% log %>% sum %>% `-` %>% return
  }
  
  # Loop through all time steps and return:
  rbind(Y, X) %>% apply(2, helper, p=p) %>% return
}

detectionNumObs1D <- function(Y, X, p, logy=F){
  # purpose : A simple model which ignores the probability of the distribution
  #           of sizes, and says that the probability of observing what we
  #           observe is simply the probability of detecting the number of 
  #           animals we did, given the true distribution.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  trueSizes <- sum(X)
  obsSizes <- sum(Y)
  
  return(dbinom(obsSizes, trueSizes, p, log=logy))
}

detectionNumObs <- function(Y, X, p, logy=T){
  # purpose : A simple model which ignores the probability of the distribution
  #           of sizes, and says that the probability of observing what we
  #           observe is simply the probability of detecting the number of 
  #           animals we did, given the true distribution.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  return(dbinom(colSums(Y), colSums(X), p, log=logy))
  
}


poissonObs <- function(Y, X, p, logy=T){
  # purpose : A model which assumes the number of individuals in each size class
  #           is poisson distributed with mean equal to the number of true
  #           individuals in the classx the detection probability.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  results <- dpois(Y, (X*p), log=logy)
  # results %<>% replace(is.infinite(results), log(.Machine$double.xmin))
  if (logy) return(colSums(results,na.rm = T))
  else return(colprods(results,na.rm = T))
}

prob1Obs <- function(Y, X, p, logy=T){
  # purpose : returns a uniform probability so that we can test the particle 
  #           filter more easily
  # inputs  : Y   - The observation matrix (not used)
  #           X   - The True state matrix (not used)
  #           p   - The probability of detection (not used)
  #           log - if TRUE, returns the log probability instead
  Y %<>% as.matrix %>% ncol
  if (logy) return(runif(Y) %>% log)
  else return(runif(Y))
}

SplitMNom<-function(vec,logy) dmultinom(x = vec[,1],prob=vec[,2],log=logy)

dmnom<-function(Y,multiProbs,logy){
  out1<-apply(abind(Y,multiProbs,along=3),2,SplitMNom,log=logy)
}

multinomialObs <- function(Y, X, p, logy=T){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability. Then, given the total number of 
  #           observed animals, we assumed that they are multinomially
  #           distributed according to the same size distribution as the true
  #           population
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # Sometimes the particle filter requires all arguments to be a list, since
  # p is only a single parameter, this can cause it to be passed as a list of 
  # one element, this tries to elegantly cope with the edge case:
  p %<>% as.numeric
  
  # Some required quantities:
  truePopSizes <- colSums(X)
  obsPopSizes <- colSums(Y)
  multiProbs <- apply(X, 2, countsToProbs)
  # multiProbs[is.na(multiProbs)]<-1
    
  # out1<- dmultinomial(t(Y), prob = t(multiProbs), log = logy)
  indy<-truePopSizes!=0
  out1<-out2<-rep(ifelse(logy,-Inf,0),length(truePopSizes))
  out1[indy]<-dmnom(Y[,indy],multiProbs[,indy],logy)
  out2[indy]<-dbinom(obsPopSizes[indy], truePopSizes[indy], p, log = logy)

  if (logy) return(out1+out2)
  else return(out1*out2)
    
  # if (logy) {
  #   output<-out1+out2
  #   output %>% replace(is.infinite(output), log(.Machine$double.xmin)) %>% return()
  # }
  # else {
  #   output<-out1*out2
  #   output %>% replace(is.infinite(output), .Machine$double.xmin) %>% return()
  # }
  
}

beta_binObs <- function(Y, X, p, logy=T, time, NoParts){
  # purpose : A model which assumes the number of individuals in each size class
  #           is poisson distributed with mean equal to the number of true
  #           individuals in the classx the detection probability.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  pobs<-rbeta(ncol(Y),p[1],p[2])
  
  if (logy) return(dbinom(colSums(X), colSums(Y), pobs, log=logy))
  else return(dbinom(colSums(X), colSums(Y), pobs, log=logy))
}

beta_poisObs <- function(Y, X, p, logy=T, time, NoParts){
  # purpose : A model which assumes the number of individuals in each size class
  #           is poisson distributed with mean equal to the number of true
  #           individuals in the classx the detection probability.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  pobs<-rbeta(ncol(Y),p[1],p[2])
  results <- dpois(Y, (X*pobs), log=logy)
  # results %<>% replace(is.infinite(results), log(.Machine$double.xmin))
  if (logy) return(colSums(results))
  else return(colprods(results))
}

beta_mnomObs<-function(Y, X, p, logy=T, time, NoParts){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability. Then, given the total number of 
  #           observed animals, we assumed that they are multinomially
  #           distributed according to the same size distribution as the true
  #           population, using a beta distribution to generate the observed probabilities
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # Sometimes the particle filter requires all arguments to be a list, since
  # p is only a single parameter, this can cause it to be passed as a list of 
  # one element, this tries to elegantly cope with the edge case:
  
  # Some required quantities:
  truePopSizes <- colSums(X)
  obsPopSizes <- colSums(Y)
  multiProbs <- apply(X, 2, countsToProbs)
  indy<-truePopSizes!=0
  
  pobs<-rbeta(length(truePopSizes),p[1],p[2])
  
  # out1<- dmultinomial(t(Y), prob = t(multiProbs), log = logy)
  # Make sure that if any X distributions have zero individuals, the likelihood equals zero
  out1<-out2<-rep(ifelse(logy,-Inf,0),length(truePopSizes))
  out1[indy]<-dmnom(Y[,indy],multiProbs[,indy],logy)
  out2[indy]<-dbinom(obsPopSizes[indy], truePopSizes[indy], pobs[indy], log = logy)
  
  if (logy) return(out1+out2)
  else return(out1*out2)
  
}

fixedMuObs<-function(Y, X, p, logy=T, time, NoParts){
  return(multinomialObs(Y[,rep(time,NoParts)], X, p[time], logy=logy))
}

fixedPoisObs<-function(Y, X, p, logy=T, time, NoParts){
  return(poissonObs(Y[,rep(time,NoParts)], X, p[time], logy=logy))
}

fixedBinObs<-function(Y, X, p, logy=T, time, NoParts){
  return(detectionNumObs(Y[,rep(time,NoParts)], X, p[time], logy=logy))
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
    array(c(reppie,reppie,reppie,reppie,reppie),
          dim = c(length(sizes),5),
          dimnames = list(round(sizes,digits = 2),
                          c("NoSurv","NoAlive","NoParents","avSurvOff","NoOff"))))
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
    return(array(c(vectorToCounts(c(newSizesI),breaks),
                   newCounts,
                   reprCounts,
                   newCounts,
                   reppie),
          dim = c(length(sizes),5),
          dimnames = list(round(sizes,digits = 2),
                          c("NoSurv","NoAlive","NoParents","avSurvOff","NoOff"))))
  } else {
    # Modify the number of reproducing sheep to retain females only
    if(oneSex) {
      # Are all sheep included in the study female?
      offCounts<-rbinom(D, reprCounts, 0.5)
    }
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
  return(array(c(vectorToCounts(c(newSizesI),breaks),
                 vectorToCounts(c(offSizes, newSizes), breaks),
                 reprCounts,
                 newCounts,
                 vectorToCounts(c(offSizes),breaks)),
               dim = c(length(sizes),5),
               dimnames = list(round(sizes,digits = 2),
                               c("NoSurv","NoAlive","NoParents","avSurvOff","NoOff"))))
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
        dim=c(length(SamplerArgs$sizes),5,ncol(initialStates)))
  
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
################### PARTICLE FILTER AND PMCMC IMPLEMENTATION ###################
# Minkowski distance function
# Distances were taken from https://www.biorxiv.org/content/10.1101/2021.07.29.454327v1.full.pdf
# Note that if p=1 we are using the L1 distances - our default
Minkowski<-function(sest,sobs,dimmie,p=1){
  # Median of columns
  meds<-apply(sest,1,median,na.rm=T)
  # Median Absolute Deviation (MAD) to the sample
  MAD<-abs(sest-meds)
  # Median Absolute Deviation to the Observation (MADO)
  MADO<-abs(sest-sobs)
  # Don't punish the summary statistics that deviate more from obs data initially than others:
  if(sum(MADO>2*MAD,na.rm = T)/length(sobs)<1/3) PCMAD<-MAD+MADO else PCMAD<-MAD
  # Calculate the Minkowski distance per summary statistic
  d_i<--abs(sest-sobs)/PCMAD
  # Find number of infinite and NaN values to artificially add to the LL:  
  infies<-length(d_i[is.infinite(d_i) | is.na(d_i)])
  infiesSW<-apply(d_i,2,function(dd) length(dd[is.infinite(dd) | is.na(dd)]))
  # output total distance
  return(list(shat=meds,
              d=pracma::nthroot(sum(d_i[!is.infinite(d_i)]^p,na.rm = T),p)*infies,
              sw=apply(d_i,2,function(dd) pracma::nthroot(sum(dd[!is.infinite(dd)]^p,na.rm = T),p))*infiesSW))
}
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
  t <- ncol(Sd)
  # Setup initial states
  prevStates <- do.call(mu, muPar)
  # Setup weight matrix and standardised weight matrix:
  output <- list(d=1, sw=rep(1.,NoParts),
                 shat=array(NA, dim=c(nrow(Sd),5,t)))
  # Update weights for first time step:
  wArgs <- list(Sd=Sd, pobs=obsProbPar, NoParts=NoParts)
  
  for (time in 1:t){
    wArgs$time<-time
    # Importance resample from the previous states 
    particleIndices <- sample(1L:NoParts, NoParts, replace = T, 
                              prob = (output$sw)/min(output$sw,na.rm = T))
    sampledStates <-  prevStates[, particleIndices]
    # IPM push forward
    wArgs$Sstar <- sampleState(sampledStates, sampleStatePar)
    # Convert to observed from latent space & calculate the objective function vector
    output <- obsProb(output,wArgs)
  }
  output$sw<-NULL
  return(output)
  
}

##################################################
### TESTING ZONE - PLEASE DELETE WHEN POSSIBLE ###
##################################################
# 
# directory<-paste0(getwd(),"/")
# source(paste0(directory,"Rcode/SimulateData.R"))
# source(paste0(directory,"Rcode/piecemealFunctions.R"))
# IPMLTP<-readRDS(paste0(directory,"RDSobjects/IPMLTP"))
# library(ggplot2)
# #
# # # Y is [size_distribution , year]
# #
# # # list(
# # #   priorFunc = 'evalPriors',
# # #   priors = priorsIPM,
# # #   priorPars = flatPriors,
# # #   skeleton = skeleton,
# # #   survFunc = 'linLogit',
# # #   growthSamp = 'sampleNorm',
# # #   reprFunc = 'linLogit',
# # #   offNumSamp = 'returnConstant',
# # #   offSizeSamp = 'sampleNorm',
# # #   oneSex = TRUE,
# # #   mu = 'multnomMu',
# # #   muPar = c(sum(Y[,1]), list(priorProbs)),
# # #   NoParts = 250,
# # #   Y = Y,
# # #   obsProb = 'detectionNumObs',
# # #   shift = qlogis(0.49),
# # #   breaks = breaks
# # # )
# #
# #
# #
# #
# #
# # ##################################################################################
# breaks <- seq(1.5, 3.55, l=6)[-2]
# shift = qlogis(0.49)
# 
# priorProbs<-c(0.03374183, 0.16805116, 0.39873340, 0.39947362)
# 
# startValues <- list(
#   survPars = c(-9.65, 3.77),
#   growthPars = c(1.41, 0.56, log(0.08)),
#   reprPars = c(-7.23, 2.6),
#   offNumPars = 1,
#   offSizePars = c(0.36, 0.71, log(0.16)),
#   Schild = qlogis(0.873),
#   obsProbPar = 10 # not too close to 50 since this will hurt the chain
# ) %>% unlist()
# 
# SVjitter<-readRDS(paste0(directory,"RDSobjects/startvals_jitter"))
# 
# simmedData<-readRDS(paste0(directory,"RDSobjects/simmedData"))
# max.cens <- simmedData$census.number %>% max
# Y <- getSizeDistns(simmedData, breaks)[,(max.cens-12):max.cens]
# 
# # mu<-match.fun('poissonMu')
# # muPar <- list(popSize=Y[,1])
# mu<-match.fun('multnomMu')
# muPar <- list(popSize=sum(Y[,1]), probs=priorProbs)
# 
# obsProb<-match.fun("detectionNumObs")
# # obsProb<-match.fun("multinomialObs")
# # obsProb<-match.fun("poissonObs")
# obsProbPar <- 10
# 
# # NoParts<-IPMLTP$NoParts
# NoParts <- 1600
# 
# survFunc <- linLogit; survPars <- c(-9.65, 3.77)
# growthSamp <- sampleNorm; growthPars <- c(1.41, 0.56, log(0.08), 1.5, 3.55)
# reprFunc <- linLogit; reprPars <- c(-7.23, 2.6);
# offNum<-1;
# offNumSamp <- match.fun('returnConstant'); offNumPars <- offNum;
# offSizeSamp <- sampleNorm; offSizePars <- c(0.36, 0.71, log(0.16), 1.5, 3.55);
# Schild<-qlogis(0.873); oneSex<-T;
# 
# D <- length(breaks)-1
# sizes <- breaks[-(D+1)] + plogis(shift)*diff(breaks)
# 
# sampleStatePar <- list(survFunc = survFunc, survPars = survPars,
#                        growthSamp = growthSamp, growthPars = growthPars,
#                        reprFunc = reprFunc, reprPars = reprPars,
#                        offNumSamp = offNumSamp, offNumPars = offNumPars,
#                        offSizeSamp = offSizeSamp, breaks = breaks,
#                        offSizePars = offSizePars, Schild=Schild,
#                        oneSex = oneSex, shift = shift, D=D,sizes=sizes)
# 
# # tmp<-Sheep[1500:2500,-1] %>% apply(2, mean,na.rm=T) %>% relist(skeleton=IPMLTP$skeleton)
# tmp<-SVjitter %>% relist(skeleton=IPMLTP$skeleton)
# # tmp<-startValues %>% relist(skeleton=IPMLTP$skeleton)
# obsProbPar <- tmp$obsProbPar
# tmp$obsProbPar<-NULL
# sampleStatePar <- append(list(survFunc = survFunc,
#                        growthSamp = growthSamp,
#                        reprFunc = reprFunc,
#                        offNumSamp = offNumSamp,
#                        offSizeSamp = offSizeSamp, breaks = breaks,
#                        oneSex = oneSex, shift = shift, D=D,sizes=sizes),tmp)
# 
# ll <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, NoParts = NoParts, returnW = T)
# # 
# # 
# # 
# # 
# sampleStatePar$Schild %<>% as.numeric %>% plogis
# obsProbPar %<>% as.numeric %>% plogis
# 
# # If half the children are of a sex not being tracked:
# if (isTRUE(sampleStatePar$oneSex)) sampleStatePar$Schild %<>% `/`(2)
# 
# breaks<-sampleStatePar$breaks
# sampleStatePar$D<-length(breaks)-1
# sampleStatePar$sizes <- breaks[-(sampleStatePar$D+1)] + plogis(as.numeric(sampleStatePar$shift))*diff(breaks)
# sampleStatePar$shift<-NULL
# sampleStatePar$breaks[c(1, sampleStatePar$D+1)] <- c(-Inf, Inf)
# 
# sampleStatePar$growthPars[3]%<>%exp()
# sampleStatePar$offSizePars[3]%<>%exp()
# 
# muPar$n<-NoParts; muPar$pobs<-obsProbPar
# prevStates <- do.call(mu, muPar)
# 
# previousState<-prevStates[,1]
# survPars<-sampleStatePar$survPars
# growthPars<-sampleStatePar$growthPars
# reprPars<- sampleStatePar$reprPars
# offNumPars<-sampleStatePar$offNumPars
# offSizePars<-sampleStatePar$offSizePars
# Schild<-sampleStatePar$Schild
# 
# ###################################
# ########## IPM Algorithm ##########
# ###################################
# sizesamp<-seq(min(sizes),max(sizes),length.out=100)
# 
# # Probability of surviving for each size class
# plot(sizesamp,survFunc(sizesamp,survPars))
# survFunc(sizes, survPars)
# # Number that make it through per class
# survFunc(sizes, survPars) %>% rbinom(n=D, size=previousState)
# # Weights of surviving sheep that survive (number of unique weights = number of size classes)
# newSizes <- survFunc(sizes, survPars) %>% rbinom(n=D, size=previousState) %>%
#   rep(x=sizes)
# # What weight do they grow to?
# newSizes %<>% growthSamp(growthPars)
# plot(sizesamp,growthSamp(sizesamp,growthPars))
# # Are these sheep likely to reproduce?
# reprProbs <- reprFunc(newSizes, reprPars)
# plot(sizesamp,reprFunc(sizesamp,reprPars))
# # What size are the newborns?
# reprSizes <- weightedSelection(newSizes, reprProbs)
# ### The rest, have a play!
# if (length(reprSizes)==0){ return(vectorToCounts(c(newSizes), breaks))
# }else{
#   offSizesTemp <- offNumSamp(reprSizes, offNumPars) %>% rep(x=reprSizes) %>%
#     offSizeSamp(offSizePars)
#   plot(sizesamp,offNumSamp(sizesamp, offNumPars))
#   # Get the new distribution of sizes:
#   # weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp))) %>%
#   #   c(newSizes) %>% vectorToCounts(breaks =  breaks) %>% return
#   offSizes <- weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp)))
# }
# vectorToCounts(c(offSizes, newSizes), breaks) %>% return
# ###################################
# ###################################
# ###################################
# 
# 
# 
# 
# # 
# # ll %>% apply(2,countsToProbs) %>% `^`(2) %>% apply(2, sum) %>% `^`(-1) %>% `[`(-1) %>% mean
# # muPar$n<-NoParts; muPar$pobs<-obsProbPar
# sampledStates<-do.call(mu, list(popSize=muPar$popSize, probs=muPar$probs,n=NoParts,pobs=obsProbPar))
# nextState<-vectorisedSamplerIPM(prevStates,sampleStatePar)

# ll<-numeric(length = 50)
# ptm <- proc.time()
# for (i in 1:50){
#     ll[i] <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, NoParts = NoParts, returnW = F)
# }
# ptm_fin<-proc.time() - ptm; ptm_fin/50

########## pfMCMC : ###########
# Original ~55 seconds

# user   system  elapsed 
# 38.83666  0.01482 38.90688 

#   user   system  elapsed 
# 23.32004  0.02008 28.83988 

# user   system  elapsed 
# 26.17082  0.01480 26.28622


# ptm <- proc.time()
# ll1 <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, NoParts = 500, returnW = T)
# ptm_fin<-proc.time() - ptm; ptm_fin
# print("Round 1")



# ptm <- proc.time()
# ll2 <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, NoParts = 2000, returnW = T)
# ptm_fin<-proc.time() - ptm; ptm_fin
# print("Round 2")
# ptm <- proc.time()
# ll3 <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, NoParts = 10000, returnW = T)
# ptm_fin<-proc.time() - ptm; ptm_fin
# print("Round 3")
# 
# for (i in 1:6){
#   plot(ll1$sw[,i])
#   points(ll2$sw[,i],col="red")
#   points(ll3$sw[,i],col="green")
# }
# 
  
##################################################################################
# 
# mu<-match.fun('multnomMu')
# muPar <- list(popSize=sum(Y[,1]), probs=priorProbs)
# 
# lister<-seq.int(from = 50,to = 1000,by = 50)
# ##### Binomial Observation Model #####
# obsProb<-match.fun("detectionNumObs")
# 
# i<-1
# mnll1<-c()
# for (bb in lister){
#   mnll1 <- c(mnll1,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                             sampleState = vectorisedSamplerIPM,
#                             sampleStatePar = sampleStatePar,
#                             obsProbPar = obsProbPar, NoParts = bb, returnW = F))
#   i<-i+1
# }
# 
# ##### Poisson Observation Model #####
# obsProb<-match.fun("poissonObs")
# 
# i<-1
# mnll2<-c()
# for (bb in lister){
#   mnll2 <- c(mnll2,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                               sampleState = vectorisedSamplerIPM,
#                               sampleStatePar = sampleStatePar,
#                               obsProbPar = obsProbPar, NoParts = bb, returnW = F))
#   i<-i+1
# }
# 
# ##### Multinomial Observation Model #####
# obsProb<-match.fun("multinomialObs1D")
# 
# i<-1
# mnll3<-c()
# for (bb in lister){
#   mnll3 <- c(mnll3,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                               sampleState = vectorisedSamplerIPM,
#                               sampleStatePar = sampleStatePar,
#                               obsProbPar = obsProbPar, NoParts = bb, returnW = F))
#   i<-i+1
# }
# 
# mu<-match.fun('poissonMu')
# muPar <- list(popSize=Y[,1])
# 
# lister<-seq.int(from = 50,to = 1000,by = 50)
# ##### Binomial Observation Model #####
# obsProb<-match.fun("detectionNumObs")
# 
# i<-1
# pnll1<-c()
# for (bb in lister){
#   pnll1 <- c(pnll1,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                               sampleState = vectorisedSamplerIPM,
#                               sampleStatePar = sampleStatePar,
#                               obsProbPar = obsProbPar, NoParts = bb, returnW = F))
#   i<-i+1
# }
# 
# ##### Poisson Observation Model #####
# obsProb<-match.fun("poissonObs")
# 
# i<-1
# pnll2<-c()
# for (bb in lister){
#   pnll2 <- c(pnll2,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                               sampleState = vectorisedSamplerIPM,
#                               sampleStatePar = sampleStatePar,
#                               obsProbPar = obsProbPar, NoParts = bb, returnW = F))
#   i<-i+1
# }
# 
# ##### Multinomial Observation Model #####
# obsProb<-match.fun("multinomialObs1D")
# 
# i<-1
# pnll3<-c()
# for (bb in lister){
#   pnll3 <- c(pnll3,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                               sampleState = vectorisedSamplerIPM,
#                               sampleStatePar = sampleStatePar,
#                               obsProbPar = obsProbPar, NoParts = bb, returnW = F))
#   i<-i+1
# }
# 
# ##################################################################################
# ##### CHECK VARIANCE CONVERGENCE #####
# mulist<-c('multnomMu','poissonMu')
# muParlist<-list(list(popSize=sum(Y[,1]), probs=priorProbs),list(popSize=Y[,1]))
#   
# DFf<-data.frame()
# i<-1
# for (k in 1:2){
#   mu<-match.fun(mulist[k])
#   muPar<-muParlist[[k]]
#   
#   for(obsProb in c("detectionNumObs","poissonObs","multinomialObs1D")){
#     # for(obsProb in c("poissonObs","multinomialObs1D")){
#     ll<-tmp2<-c()
#     op<-match.fun(obsProb)
#     # bb<- seq.int(from = 50,to = 1000,by = 50)
#     bb<-c(10,50,200,500,2000,5000,8000)
#     # bb<-c(20,100)
#     for (NoParts in bb){
#       for (j in 1:50){
#         ll <- c(ll,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = op,
#                                   sampleState = vectorisedSamplerIPM,
#                                   sampleStatePar = sampleStatePar,
#                                   obsProbPar = obsProbPar, NoParts = NoParts, returnW = F))
#         tmp2<-c(tmp2,NoParts)
#       }
#     }
#     
#     DF<-data.frame(ll=ll,bb=tmp2,obsProb=rep(obsProb,length(tmp2)),mu=rep(mulist[k],length(tmp2)))
#     # DFsd<-DF%>%group_by(bb)%>%summarise(lsd=sd(ll))
#     DFf<-rbind(DFf,DF)
#     
#     i<-i+1
#   }
# }
# 
# p<-ggplot(DFf,aes(bb,-ll,group=mu))+geom_point(aes(colour=obsProb,shape=mu)) + scale_y_log10() +
#   xlab("#Particles") + ylab("Negative Log-Likelihood")
# ggsave("pf_conv_nll.png", plot=p,path = paste0(directory,'Plots/NEW_WORK/'),width = 8,height = 5)
# 
# DFsd<-DFf%>%group_by(mu,obsProb,bb)%>%summarise(lsd=sd(ll))
# p<-ggplot(DFsd,aes(bb,lsd,group=mu))+geom_point(aes(colour=obsProb,shape=mu))  + scale_y_log10() + scale_x_log10() +
#   xlab("#Particles") + ylab("Standard Deviation of LL")
# ggsave("pf_conv_sd.png", plot=p,path = paste0(directory,'Plots/NEW_WORK/'),width = 8,height = 5)
# ##################################################################################
# ##### CONVERGENCE OF NLL & SD(LL) #####
# 
# convDF<-DFf%>%filter(bb==8000)%>%group_by(mu,obsProb)%>%summarise(mean=mean(ll))
# 
# predy<-DFsd%>%filter(mu=="multnomMu" & obsProb=="detectionNumObs")%>%lm(formula = log(bb) ~ log(lsd), data=.) 
# print(exp(predict(object = predy,data.frame(lsd=c(1.2)))))
# # 1422.193
# print(convDF$mean[1])
# # -22.60212
# 
# predy<-DFsd%>%filter(mu=="multnomMu" & obsProb=="poissonObs")%>%lm(formula = log(bb) ~ log(lsd), data=.) 
# print(exp(predict(object = predy,data.frame(lsd=c(1.2)))))
# # 156.0012
# print(convDF$mean[2])
# # -101.4621
# 
# predy<-DFsd%>%filter(mu=="multnomMu" & obsProb=="multinomialObs1D")%>%lm(formula = log(bb) ~ log(lsd), data=.) 
# print(exp(predict(object = predy,data.frame(lsd=c(1.2)))))
# # 6104.29
# print(convDF$mean[3])
# # -95.32274
# 
# predy<-DFsd%>%filter(mu=="poissonMu" & obsProb=="detectionNumObs")%>%lm(formula = log(bb) ~ log(lsd), data=.) 
# print(exp(predict(object = predy,data.frame(lsd=c(1.2)))))
# # 1370.334 
# print(convDF$mean[4])
# # -22.93806
# 
# predy<-DFsd%>%filter(mu=="poissonMu" & obsProb=="poissonObs")%>%lm(formula = log(bb) ~ log(lsd), data=.) 
# print(exp(predict(object = predy,data.frame(lsd=c(1.2)))))
# # 155.0295
# print(convDF$mean[5])
# # -99.82453
# 
# predy<-DFsd%>%filter(mu=="poissonMu" & obsProb=="multinomialObs1D")%>%lm(formula = log(bb) ~ log(lsd), data=.) 
# print(exp(predict(object = predy,data.frame(lsd=c(1.2)))))
# # 6364.99
# print(convDF$mean[6])
# # -92.7514
  
##################################################################################

# predict(data.frame(lsd=c(log(1.2))),data=.)%>%exp()
# library(mgcv)
# model <- gam(lsd ~ s(bb,k = 4), data = DFsd)
# # model <- glm(lsd ~poly(bb,4), data = DFsd)
# summary(model)
# predy <- model %>% predict(data.frame(bb=bb))
# pd<-data.frame(bb=bb,pred=exp(predy))
# varLL<-ggplot(DFsd,aes(x=bb,y=lsd)) + geom_point() + geom_line(data = pd,aes(x=bb,y=pred),colour="red") +
#   geom_hline(yintercept = 1.2)+ geom_text(x=7000,y=log10(4),label=paste0("var=1.2 at ",bb[which.min(abs(exp(predy)-1.2))]," particles")) +
#   scale_y_log10()+xlab("Number of pfMCMC Particles") + ylab("Variance in log-likelihood (n=150)")
# NLL<-ggplot(DFold,aes(x=bb,y=ll))+geom_point()+ylim(c(23,50))+ylab("Negative Log-Likelihood") +
#   xlab("Number of pfMCMC Particles")+geom_smooth(method="auto",se = FALSE)
# gridExtra::grid.arrange(NLL,varLL,nrow=1)

# print(paste0("Number of pfMCMC particles = ",bb[which.min(abs(exp(predy)-1.2))]))

##################################################################################






# 
# sampleStatePar <- list(survFunc = survFunc, survPars = survPars,
#                            growthSamp = growthSamp, growthPars = growthPars,
#                            reprFunc = reprFunc, reprPars = reprPars,
#                            offNumSamp = offNumSamp, offNumPars = offNumPars,
#                            offSizeSamp = offSizeSamp, breaks = breaks,
#                            offSizePars = offSizePars, Schild=(plogis(Schild)*0.5),
#                            oneSex = oneSex, shift = plogis(shift), D=D,sizes=sizes)
# breaks<-sampleStatePar$breaks
# sampleStatePar$D<-length(breaks)-1
# sampleStatePar$sizes <- breaks[-(sampleStatePar$D+1)] + sampleStatePar$shift*diff(breaks)
# sampleStatePar$breaks[c(1, D+1)] <- c(-Inf, Inf)
# breaks[c(1, D+1)] <- c(-Inf, Inf)
# sampleStatePar$growthPars[3]%<>%exp()
# growthPars<-sampleStatePar$growthPars
# sampleStatePar$offSizePars[3]%<>%exp()
# offSizePars<-sampleStatePar$offSizePars
# 
# sampledStates<-do.call(mu, c(muPar, list(NoParts)))
# subsamp<-sampledStates[,50]
# 
# sampleStateIPM_red(subsamp, survFunc, survPars,
#                growthSamp, growthPars, reprFunc, reprPars,
#                offNumSamp, offNumPars, offSizeSamp, offSizePars,
#                (plogis(Schild)*0.5), breaks, oneSex = TRUE, checks = FALSE,
#                verbose = FALSE, shift = plogis(shift), D=D,sizes=sizes)
# 
# 
# ptm <- proc.time()
# for (i in 1:50){
#   prevStates<-vectorisedSamplerIPM(sampledStates,sampleStatePar)  
# }
# ptm_fin<-proc.time() - ptm; ptm_fin/50


# user   system  elapsed 
# 13.70282  0.05876 14.08398 

# user   system  elapsed 
# 10.21486  0.01654 10.48674 

# user  system elapsed 
# 6.00668 0.00976 6.05336 

# user  system elapsed 
# 6.42244 0.00230 6.44178 

# user  system elapsed 
# 5.63376 0.00232 5.64658

# user  system elapsed 
# 4.04028 0.00084 4.04284

# user  system elapsed 
# 5.38704 0.00464 5.41094 


