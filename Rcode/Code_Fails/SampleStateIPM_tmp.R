
sampleNorm <- function(x, pars){
  # purpose : produces samples as DTN without the double truncation.
  intercept<-pars[1];gradient<-pars[2];sd<-pars[3];n<-pars[4]
  shitz<-intercept+gradient*x
  val<-rnorm(n = n, mean = shitz , sd = sd)
  if(is.na(val)) stop(paste0("error: par = ",unlist(pars),", x =",x))
  return(val)
  # return(rnorm(n = pars[4], mean = (x*pars[2] + pars[1]) , sd = pars[3]))
}

linLogit <- function(z, par){
  # purpose : Calculates the probability an individual of size x survives to 
  #           the next time step. Uses a logistic form.
  # inputs  : z         - The size of the individual (continuous)
  #           par       - A vector with entries 'intercept' and 'gradient'
  # output  : The density of survival
  
  intercept <- par[1]
  gradient <- par[2]
  shit<-intercept+gradient*z
  return(plogis(shit))
  # calculate probability of survival:
  # (par[1] + par[2]*z) %>% plogis() %>% return
}

####################### IPM STATE SPACE SAMPLER ################################
sampleStateIPM_red <- function(previousState, survFunc, survPars,
                           growthSamp, growthPars, reprFunc, reprPars, 
                           offNumSamp, offNumPars, offSizeSamp, offSizePars,
                           Schild, breaks, oneSex = TRUE, checks = FALSE,
                           verbose = FALSE, shift = qlogis(0.5), D=4, sizes){
  
  # To insure function names work too:
  # survFunc %<>% match.fun ; growthSamp %<>% match.fun ; reprFunc %<>% match.fun
  # offNumSamp %<>% match.fun ; offSizeSamp %<>% match.fun
  
  # # Get midpoint of each user defined interval:
  # D <- length(previousState)
  # sizes <- breaks[-(D+1)] + shift*diff(breaks)
  
  # Get the size distribution of individuals that survive:
  newSizes <- survFunc(sizes, survPars) %>% rbinom(n=D, size=previousState) %>%
    rep(x=sizes) %>% growthSamp(growthPars)
  if (length(newSizes)==0) return(rep(0, D))
  
  # Get the size distribution of newborns:
  reprProbs <- reprFunc(newSizes, reprPars)
  reprSizes <- weightedSelection(newSizes, reprProbs)
  if (length(reprSizes)==0){ return(vectorToCounts(c(newSizes), breaks))
  }else{
    offSizesTemp <- offNumSamp(reprSizes, offNumPars) %>% rep(x=reprSizes) %>%
      offSizeSamp(offSizePars)
    # Get the new distribution of sizes:
    # weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp))) %>% 
    #   c(newSizes) %>% vectorToCounts(breaks =  breaks) %>% return
    offSizes <- weightedSelection(offSizesTemp,rep(Schild,length(offSizesTemp)))
  }
  
  vectorToCounts(c(offSizes, newSizes), breaks) %>% return
  
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
  initialStates %>% apply(2, helper) %>% return
}