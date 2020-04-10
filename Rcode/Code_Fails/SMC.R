library(dissPackage3) #for the countsToProbs function

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

multnomMu <- function(popSize, probs, n=1){
  # purpose : returns a vector of individual counts which is uniformly 
  #           distributed across the dimensions of the state space
  # inputs  : popSize - The integer total population size (to be spread
  #                     uniformly amongst the train classes) 
  #           probs   - The probabilities of beingin each size class
  return(rmultinom(n, popSize, probs))
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

detectionNumObs1D <- function(Y, X, p, log=F){
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
  
  return(dbinom(obsSizes, trueSizes, p, log=log))
}

detectionNumObs <- function(Y, X, p, log=F){
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
  
  return(dbinom(colSums(Y), colSums(X), p, log=log))
  
}

poissonObs <- function(Y, X, p, log=T){
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
  
  # In case Y and X are for a single time step:
  Y %<>% as.matrix ; X %<>% as.matrix
  
  helper <- function(yx, log=T){
    n <- length(yx)
    y <- yx[1:(n/2)]
    x <- yx[(n/2+1):n]*p
    results <- dpois(y, x, log=log)
    if (log) return(sum(results))
    else return(prod(results))
  }
  
  rbind(Y, X) %>% apply(2, helper) %>% return
}

prob1Obs <- function(Y, X, p, log=F){
  # purpose : returns a uniform probability so that we can test the particle 
  #           filter more easily
  # inputs  : Y   - The observation matrix (not used)
  #           X   - The True state matrix (not used)
  #           p   - The probability of detection (not used)
  #           log - if TRUE, returns the log probability instead
  Y %<>% as.matrix %>% ncol
  if (log) return(runif(Y) %>% log)
  else return(runif(Y))
}

multinomialObs <- function(Y, X, p, log=F){
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
  if (class(p)=='list') p %<>% as.numeric
  
  helper <- function(yp){
    n <- length(yp)/2
    y <- yp[1:n]
    pr <- yp[-(1:n)]
    dmultinom(y, prob = pr, log = log)
  }
  
  if (isTRUE(log)) operator <- `+`
  else operator <- `*`
  
  # In case Y and X are for a single time step:
  Y %<>% as.matrix ; X %<>% as.matrix
  
  # Some required quantities:
  truePopSizes <- apply(X, 2, sum)
  obsPopSizes <- apply(Y, 2, sum)
  multiProbs <- apply(X, 2, countsToProbs)
  
  # Calculate the prob of seeing the distribution we saw given the number we saw
  # and the probability of seeing the number we saw
  output <- rbind(Y, multiProbs) %>% apply(2, helper) %>%
    operator(dbinom(obsPopSizes, truePopSizes, p, log = log))
  
  # To avoid returning -Inf log probabilities:
  output %<>% replace(is.infinite(output), log(.Machine$double.xmin))
  return(output)
}


####################### IPM STATE SPACE SAMPLER ################################
sampleStateIPM <- function(previousState, survFunc, survPars,
                                growthSamp, growthPars, reprFunc, reprPars, 
                                offNumSamp, offNumPars, offSizeSamp, offSizePars,
                                Schild, breaks, oneSex = TRUE, checks = FALSE,
                                verbose = FALSE, shift = qlogis(0.5)){
  # purpose : A barebones version of simulateIBM which takes as input a size 
  #           distribution at a given time, and given the parameters of the 
  #           vital rate functions, produces a sample of the size distribution
  #           of the individuals at the next time step (census).
  # inputs  : previousState - A vector of the number of individuals in each
  #                           size class. Can be any length.
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
  #           shift         - In many cases, if there is a correlation between
  #                           size and survival, assuming that individuals are
  #                           at the midpoint of the interval can cause bias,
  #                           leading to large shifts in population size over
  #                           time. Shift corrects for this, by putting the
  #                           assumed size of the individual right at the start
  #                           of the interval (value 0), right at the end (value
  #                           1), or anything in between.
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
  
  # Some link functions:
  Schild %<>% plogis
  shift %<>% plogis
  
  # If half the children are of a sex not being tracked:
  if (isTRUE(oneSex)) Schild %<>% `/`(2)
  
  # Get midpoint of each user defined interval:
  D <- length(previousState)
  sizes <- breaks[-(D+1)] + shift*diff(breaks)
  
  # Get the size distribution of individuals that survive:
  newSizes <- survFunc(sizes, survPars) %>% rbinom(n=D, size=previousState) %>%
    rep(x=sizes) %>% growthSamp(growthPars)
  if (length(newSizes)==0) return(rep(0, D))
  
  # Get the size distribution of newborns:
  reprProbs <- reprFunc(newSizes, reprPars)
  reprSizes <- weightedSelection(newSizes, reprProbs)
  if (length(reprSizes)==0) offSizes <- c()
  else{
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
  breaks[c(1, D+1)] <- c(-Inf, Inf)
  vectorToCounts(c(offSizes, newSizes), breaks) %>% return  
  
}

################################ REDUCED ########################################
sampleStateIPM_red <- function(previousState, survFunc, survPars,
                               growthSamp, growthPars, reprFunc, reprPars, 
                               offNumSamp, offNumPars, offSizeSamp, offSizePars,
                               Schild, breaks, oneSex = TRUE, checks = FALSE,
                               verbose = FALSE, shift = qlogis(0.5), D=4, sizes){
  
  # Get the size distribution of individuals that survive:
  newSizes <- survFunc(sizes, survPars) %>% rbinom(n=D, size=previousState) %>%
    rep(x=sizes) %>% growthSamp(c(growthPars,D))
  if (length(newSizes)==0) return(rep(0, D))
  
  # Get the size distribution of newborns:
  reprProbs <- reprFunc(newSizes, reprPars)
  reprSizes <- weightedSelection(newSizes, reprProbs)
  if (length(reprSizes)==0){ vectorToCounts(c(offSizes,newSizes), breaks) %>% return
  }else{
    offSizesTemp <- offNumSamp(reprSizes, offNumPars) %>% rep(x=reprSizes) %>%
      offSizeSamp(offSizePars)
    # Get the new distribution of sizes:
    offSizes <- weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp)))
    # weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp))) %>% 
    #   c(newSizes) %>% vectorToCounts(breaks =  breaks) %>% return    
  }
  
  vectorToCounts(c(offSizes, newSizes), breaks) %>% return
  
}

################################# UTILS ########################################

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
  initialStates %>% apply(2, helper) %>% return
}

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
  breaks %<>% sort
  breaks[c(1, length(breaks))] <- c(-Inf, Inf)
  
  # Loop through the censuses and get the size distn:
  maxCens <- DF$census.number %>% max
  sapply(1:maxCens, getSizeForCensus, DF=DF, breaks=breaks) %>% return
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

particleFilter <- function(Y, mu, muPar, sampleState, sampleStatePar, obsProb,
                           obsProbPar, b, checks=F, returnW=F, l=T, cap=-300){
  # purpose : Produces an estimate of the log likelihood of the model, given
  #           b particles projected through the state space
  # inputs  : Y              - The matrix of observations. Columns are time
  #                            steps.
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
  #           b              - The number of particles to produce.
  #           t              - The number of time steps to simulate each 
  #                            particle for.
  #           checks         - If TRUE, does some type checking on inputs.
  #           returnW        - if TRUE, returns the weights instead of the 
  #                            del Moral estimator of the likelihood
  #           l              - If TRUE, returns the log lik
  #           cap            - Any log likelihood smaller than cap for a
  #                            single observation will be rounded up to cap,
  #                            to help limit numerical errors
  # output  : A matrix or array with n rows, t columns and the third dimension
  #           equal to the dimension of the state space.
  
  # Some link functions:
  sampleStatePar$Schild %<>% plogis
  sampleStatePar$shift %<>% plogis
  obsProbPar %<>% as.numeric %>% plogis
  
  # If half the children are of a sex not being tracked:
  if (isTRUE(sampleStatePar$oneSex)) sampleStatePar$Schild %<>% `/`(2)
  
  breaks<-sampleStatePar$breaks
  sampleStatePar$D<-length(breaks)-1
  sampleStatePar$sizes <- breaks[-(sampleStatePar$D+1)] + sampleStatePar$shift*diff(breaks)
  sampleStatePar$breaks[c(1, D+1)] <- c(-Inf, Inf)
  
  sampleStatePar$growthPars<-c(sampleStatePar$growthPars[1:2],exp(sampleStatePar$growthPars[3]),sampleStatePar$D,sampleStatePar$growthPars[4:5])
  sampleStatePar$offSizePars<-c(sampleStatePar$offSizePars[1:2],exp(sampleStatePar$offSizePars[3]),sampleStatePar$D,sampleStatePar$offSizePars[4:5])
  
  # In case of a single observation:
  Y %<>% as.matrix
  
  # So function names can be passed:
  mu <- match.fun(mu) ; sampleState <- match.fun(sampleState)
  obsProb <-  match.fun(obsProb)
  
  # Setup initial state, weight matrix and standardised weight matrix:
  t <-  ncol(Y)
  initDists <- do.call(mu, c(muPar, list(b)))
  logw <- matrix(NA, nrow=b, ncol=t)
  sw <- matrix(NA, nrow=b, ncol=t)
  
  # Setup matrix to store the states for the previous time step:
  prevStates <- initDists # each particle is a column
  
  # Update weights for first time step:
  initArgs <- list(Y = Y[, rep(1, b)], X = prevStates, obsProbPar, log = T)
  # This is dbinom(observedvalues,totalpopulation,capturelikelihood)
  # Represents the likelihood of having sum(Y) observed out of X total when captureprob=p 
  logw[,1] <- do.call(obsProb, initArgs)
  
  # Remove the max for numerical stability, and then get standardised weights:
  M <- max(logw[,1]) ; residw <- logw[,1] - M
  # Sometimes removing M makes -Inf values when we exp again, so we need to 
  # reset them back to something that doesn't break the maths:
  residw %<>% replace(residw < cap, cap)
  sw[,1] <- exp(residw)
  
  # Create the LL object:
  output <- M + log(mean(sw[,1]))
  
  for (time in 2:t){
    # Sample from the previous states:
    # Using importance resampling (weighted & with replacement)
    particleIndices <- sample(1:b, b, replace = T, prob = sw[, time - 1])
    sampledStates <-  prevStates[, particleIndices]
    prevStates <- sampleState(sampledStates, sampleStatePar)
    
    # Evaluate the conditional likelihood of the data given the states:
    #@@@@@@@@@ Get rid of this unnecessary creation and combination of matrices @@@@@@@
    wArgs <- list(Y = Y[, rep(time, b)], X = prevStates, obsProbPar, log = T)
    logw[, time] <- do.call(obsProb, wArgs)
    # Note: log(.Machine$double.xmin) = -708.39
    logw[, time] %<>% `<`(cap) %>% replace(x = logw[, time], values = cap)
    
    # Update the weight matrices:
    M <- max(logw[, time])
    #@@@@@@@@@ Get rid of residw by using piping options @@@@@@@
    residw <- logw[, time] - M
    # residw %<>% replace(residw < -100, -100)
    sw[, time] <- exp(residw)
    
    # Update the LL:
    output <- output + M + log(mean(sw[, time]))
  }
  
  if (returnW) return(sw)
  else return(output)
}

##################################################
### TESTING ZONE - PLEASE DELETE WHEN POSSIBLE ###
##################################################

directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
source(paste0(directory,"Rcode/SimulateData.R"))
source(paste0(directory,"Rcode/piecemealFunctions.R"))
IPMLTP<-readRDS(paste0(directory,"RDSobjects/IPMLTP"))
library(ggplot2)

# Y is [size_distribution , year]

# list(
#   priorFunc = 'evalPriors',
#   priors = priorsIPM,
#   priorPars = flatPriors,
#   skeleton = skeleton,
#   survFunc = 'linLogit',
#   growthSamp = 'sampleNorm',
#   reprFunc = 'linLogit',
#   offNumSamp = 'returnConstant',
#   offSizeSamp = 'sampleNorm',
#   oneSex = TRUE,
#   mu = 'multnomMu',
#   muPar = c(sum(Y[,1]), list(priorProbs)),
#   b = 250,
#   Y = Y,
#   obsProb = 'detectionNumObs',
#   shift = qlogis(0.49),
#   breaks = breaks
# )





##################################################################################
Y<-IPMLTP$Y ; mu<-IPMLTP$mu ; muPar<-IPMLTP$muPar ; obsProb<-IPMLTP$obsProb ;
breaks<-IPMLTP$breaks ; shift<-IPMLTP$shift

rm(IPMLTP)

# b<-IPMLTP$b
b <- 7500

survFunc <- linLogit; survPars <- c(-9.65, 3.77)
growthSamp <- sampleNorm; growthPars <- c(1.41, 0.56, log(0.08), 1.5, 3.55)
reprFunc <- linLogit; reprPars <- c(-7.23, 2.6);
offNum<-1;
offNumSamp = match.fun('returnConstant'); offNumPars = offNum;
offSizeSamp <- sampleNorm; offSizePars <- c(0.36, 0.71, log(0.16), 1.5, 3.55);
Schild<-0.873; oneSex<-T;
obsProbPar = 10
D <- length(breaks)-1
sizes <- breaks[-(D+1)] + shift*diff(breaks)

sampleStatePar <- list(survFunc = survFunc, survPars = survPars,
                           growthSamp = growthSamp, growthPars = growthPars,
                           reprFunc = reprFunc, reprPars = reprPars,
                           offNumSamp = offNumSamp, offNumPars = offNumPars,
                           offSizeSamp = offSizeSamp, breaks = breaks,
                           offSizePars = offSizePars, Schild=Schild,
                           oneSex = oneSex, shift = shift, D=D,sizes=sizes)

# sampleStateIPM(previousState, survFunc, survPars,
#               growthSamp, growthPars, reprFunc, reprPars,
#               offNumSamp, offNumPars, offSizeSamp, offSizePars,
#               Schild, breaks, oneSex = TRUE, checks = FALSE,
#               verbose = FALSE, shift = qlogis(0.5))

# sampledStates<-do.call(mu, c(muPar, list(b)))
# prevStates<-vectorisedSamplerIPM(sampledStates,sampleStatePar)

ptm <- proc.time()
for (i in 1:50){
    ll <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
                     sampleState = vectorisedSamplerIPM,
                     sampleStatePar = sampleStatePar,
                     obsProbPar = obsProbPar, b = b, returnW = F)
}
ptm_fin<-proc.time() - ptm; ptm_fin/50

ll <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
                     sampleState = vectorisedSamplerIPM,
                     sampleStatePar = sampleStatePar,
                     obsProbPar = obsProbPar, b = b, returnW = T)
ll %>% apply(2,countsToProbs) %>% `^`(2) %>% apply(2, sum) %>% `^`(-1) %>% `[`(-1) %>% mean

########## pfMCMC : ###########
# user   system  elapsed 
# 38.83666  0.01482 38.90688 

#   user   system  elapsed 
# 23.32004  0.02008 28.83988 


# ptm <- proc.time()
# ll1 <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, b = 500, returnW = T)
# ptm_fin<-proc.time() - ptm; ptm_fin
# print("Round 1")



# ptm <- proc.time()
# ll2 <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, b = 2000, returnW = T)
# ptm_fin<-proc.time() - ptm; ptm_fin
# print("Round 2")
# ptm <- proc.time()
# ll3 <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                      sampleState = vectorisedSamplerIPM,
#                      sampleStatePar = sampleStatePar,
#                      obsProbPar = obsProbPar, b = 10000, returnW = T)
# ptm_fin<-proc.time() - ptm; ptm_fin
# print("Round 3")
# 
# for (i in 1:6){
#   plot(ll1$sw[,i])
#   points(ll2$sw[,i],col="red")
#   points(ll3$sw[,i],col="green")
# }
# 
# i<-1
# ll<-c()
# for (bb in seq.int(from = 500,to = 10000,by = 50)){
# ll <- c(ll,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                       sampleState = vectorisedSamplerIPM,
#                       sampleStatePar = sampleStatePar,
#                       obsProbPar = obsProbPar, b = bb, returnW = F))
#   i<-i+1
# }
# 
# lenny<-length(ll)
# meany<-sdy<-c()
# for(i in 6:lenny){
#   meany<-c(meany,mean(ll[(i-5):i]))
#   sdy<-c(sdy,sd(ll[(i-5):i]))
# }
# x<-tmp[6:lenny]
# plot(x,meany,ylim=range(c(meany-sdy,meany+sdy)),pch=19)
# arrows(x, meany-sdy, x, meany+sdy, length=0.05, angle=90, code=3)
# 
# ll2<-tmp2<-c()
# for (bb in c(250,500,750,1000,2500,5000,7500,10000)){
#   for (j in 1:150){
#   ll2 <- c(ll2,particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
#                             sampleState = vectorisedSamplerIPM,
#                             sampleStatePar = sampleStatePar,
#                             obsProbPar = obsProbPar, b = bb, returnW = F))
#   tmp2<-c(tmp2,bb)
#   }
# }
# 
# DF<-data.frame(ll=ll2,bb=tmp2)
# DFsd<-DF%>%group_by(bb)%>%summarise(sd=sd(ll))
# bb<- seq.int(from = 200,to = 10000,by = 100)
# 
# library(mgcv)
# model <- gam(lsd ~ s(bb,k = 4), data = DFsd)
# # model <- glm(lsd ~poly(bb,4), data = DFsd)
# summary(model)
# predy <- model %>% predict(data.frame(bb=bb))
# pd<-data.frame(bb=bb,pred=exp(predy))
# varLL<-ggplot(DFsd,aes(x=bb,y=sd)) + geom_point() + geom_line(data = pd,aes(x=bb,y=pred),colour="red") + 
#   geom_hline(yintercept = 1.2)+ geom_text(x=7000,y=log10(4),label=paste0("var=1.2 at ",bb[which.min(abs(exp(predy)-1.2))]," particles")) +
#   scale_y_log10()+xlab("Number of pfMCMC Particles") + ylab("Variance in log-likelihood (n=150)")
# # NLL<-ggplot(DFold,aes(x=bb,y=ll))+geom_point()+ylim(c(23,50))+ylab("Negative Log-Likelihood") + 
# #   xlab("Number of pfMCMC Particles")+geom_smooth(method="auto",se = FALSE)
# # gridExtra::grid.arrange(NLL,varLL,nrow=1)
# 
# print(paste0("Number of pfMCMC particles = ",bb[which.min(abs(exp(predy)-1.2))]))


##################################################################################







sampleStatePar <- list(survFunc = survFunc, survPars = survPars,
                           growthSamp = growthSamp, growthPars = growthPars,
                           reprFunc = reprFunc, reprPars = reprPars,
                           offNumSamp = offNumSamp, offNumPars = offNumPars,
                           offSizeSamp = offSizeSamp, breaks = breaks,
                           offSizePars = offSizePars, Schild=(plogis(Schild)*0.5),
                           oneSex = oneSex, shift = plogis(shift), D=D,sizes=sizes)
breaks<-sampleStatePar$breaks
sampleStatePar$D<-length(breaks)-1
sampleStatePar$sizes <- breaks[-(sampleStatePar$D+1)] + sampleStatePar$shift*diff(breaks)
sampleStatePar$breaks[c(1, D+1)] <- c(-Inf, Inf)
sampleStatePar$growthPars<-c(sampleStatePar$growthPars[1:2],exp(sampleStatePar$growthPars[3]),sampleStatePar$D,sampleStatePar$growthPars[4:5])
sampleStatePar$offSizePars<-c(sampleStatePar$offSizePars[1:2],exp(sampleStatePar$offSizePars[3]),sampleStatePar$D,sampleStatePar$offSizePars[4:5])


sampledStates<-do.call(mu, c(muPar, list(b)))
subsamp<-sampledStates[,50]
sampleStateIPM_red(previousState = subsamp,sampleStatePar)
prevStates<-vectorisedSamplerIPM(sampledStates,sampleStatePar)

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
