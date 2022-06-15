
# Skeleton frame for the parameterisation vector
skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA
)
if(oneSex) {sexprob<-0.5} else {sexprob<-1}
# Parameter support ranges
supports<-data.frame(
  lower=c(-Inf, # Survival Logistic Regression Intercept
          -Inf, # Survival Logistic Regression Gradient
          -Inf, # Growth Linear Regression Intercept
          -Inf, # Growth Linear Regression Gradient
          0,    # Growth Linear Regression Dispersion (Std. Dev.)
          -Inf, # Reproduction Logistic Regression Intercept
          -Inf, # Reproduction Logistic Regression Gradient
          1,    # Offspring Number per Birth
          -Inf, # Offspring Size Linear Regression Intercept
          -Inf, # Offspring Size Linear Regression Gradient
          0,    # Offspring Size Linear Regression Dispersion (Std. Dev.)
          0),    # Offspring Survival Probability
  upper=c(Inf, # Survival Logistic Regression Intercept
          Inf, # Survival Logistic Regression Gradient
          Inf, # Growth Linear Regression Intercept
          Inf, # Growth Linear Regression Gradient
          Inf,    # Growth Linear Regression Dispersion (Std. Dev.)
          Inf, # Reproduction Logistic Regression Intercept
          Inf, # Reproduction Logistic Regression Gradient
          Inf,    # Offspring Number per Birth
          Inf, # Offspring Size Linear Regression Intercept
          Inf, # Offspring Size Linear Regression Gradient
          Inf,    # Offspring Size Linear Regression Dispersion (Std. Dev.)
          sexprob),    # Offspring Survival Probability
)

# When using the data obervation values directly rather than assuming a distribution of observation probabilities
if(!fixedObsProb) {
  # Observed Probability Beta Shape Param 1 & 2
  IPMLTP$links%<>%c('exp','exp')
  # Make sure that the skeleton frame also includes this
  IPMLTP$skeleton %<>% c(list(obsProbPar = rep(NA,2)))
  # Add to the supports:
  supports%<>%rbind(data.frame(lower=c(0,0),upper=c(Inf,Inf)))
}
# Relist the initial values with respect to the skeleton template
x0%<>%relist(skeleton = IPMLTP$skeleton)

# Some basic link functions (to be used in MCMC calculations)
returnSelf <- function(x) x
linkNum <- function(x) exp(x)+1
invlinkNum <- function(x) log(x-1)
if(oneSex) {Schilder <- function(x) 0.5*plogis(x)
} else {Schilder <- function(x) plogis(x)}
# Generate the appropriate functions to be passed into the MCMC algorithm
links<-invlinks<-acceptTrans<-c()
for (i in 1:nrow(supports)){
  # Real numbers
  if(is.infinite(supports$lower[i]) & is.infinite(supports$upper[i])){
    links%<>%c(returnSelf)
    invlinks%<>%c(returnSelf)
    acceptTrans%<>%c(function(xold,xnew){return(1)})
    
  # x>a
  } else if(!is.infinite(supports$lower[i]) & is.infinite(supports$upper[i])){
    links%<>%c(function(x){exp(x)+supports$lower[i]})
    invlinks%<>%c(function(x){log(x-supports$lower[i])})
    acceptTrans%<>%c(function(xold,xnew){return((xnew-supports$lower[i])/(xold-supports$lower[i]))})
    
  # x<b
  } else if(is.infinite(supports$lower[i]) & !is.infinite(supports$upper[i])){
    links%<>%c(function(x){supports$upper[i]-exp(x)})
    invlinks%<>%c(function(x){log(supports$upper[i]-x)})
    acceptTrans%<>%c(function(xold,xnew){return((supports$upper[i]-xnew)/(supports$upper[i]-xold))})

  # a<x<b
  } else if(supports$lower[i]==0 & supports$upper[i]==1){
    links%<>%c(function(x){(supports$upper[i]*exp(x)+supports$upper[i])/(exp(x)+1)})
    invlinks%<>%c(function(x){log()})
    acceptTrans%<>%c(function(xold,xnew){return((supports$upper[i]-xnew)/(supports$upper[i]-xold))})
    
  } else stop("Check your model parameter support ranges: link function not defined (CodeSkeleton.R)")
  
}

# The start of the list of functions, parameters and formatting for the logTarget
IPMLTP <- list(
  skeleton = skeleton,
  links = links,
  survFunc = match.fun('linLogit'),
  growthSamp = match.fun(normsampler),
  reprFunc = match.fun('linLogit'),
  offNumSamp = match.fun('PoisNum'),
  offSizeSamp = match.fun(normsampler)
)

# Function used to sample the growth & offspring growth distribution
if(normsampler=="sampleDTN") {
  IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- doublyTruncatedNormal
}else IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- normal


