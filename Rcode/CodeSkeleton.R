
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
          sexprob)    # Offspring Survival Probability
)
# When using the data obervation values directly rather than assuming a distribution of observation probabilities
if(!fixedObsProb) {
  # Make sure that the skeleton frame also includes this
  skeleton %<>% c(list(obsProbPar = rep(NA,2)))
  # Add to the supports:
  supports%<>%rbind(data.frame(lower=c(0,0),upper=c(Inf,Inf)))
}

# Some basic link functions (to be used in MCMC calculations)
returnSelf <- function(x) x
linkNum <- function(x) exp(x)+1
invlinkNum <- function(x) log(x-1)
if(oneSex) {Schilder <- function(x) 0.5*plogis(x)
} else {Schilder <- function(x) plogis(x)}

# Generate the appropriate functions to be passed into the MCMC algorithm
links<-invlinks<-acceptTrans<-skew<-initPrior<-c()
for (i in 1:nrow(supports)){
  # Real numbers
  if(is.infinite(supports$lower[i]) & is.infinite(supports$upper[i])){
    links%<>%c(returnSelf)
    invlinks%<>%c(returnSelf)
    acceptTrans%<>%c(function(xold,xnew){return(1)})
    skew%<>%c(F)
    initPrior%<>%c(function(n,a,b) {rnorm(n,a,b)})
    
  # WARNING: R USES GLOBAL ENV FOR VARIABLES DEFINED IN FUNCTION
  # TO AVOID PROBLEMS, WE DEFINE supports$lower[i] & supports$upper[i] IN LOCAL ENV
  # x>0
  } else if(supports$lower[i]==0 & is.infinite(supports$upper[i])){
    links%<>%c(function(x){exp(x)+0})
    invlinks%<>%c(function(x){log(x-0)})
    acceptTrans%<>%c(function(xold,xnew){return((xnew-0)/(xold-0))})
    skew%<>%c(T)
  # x>1
  } else if(supports$lower[i]==1 & is.infinite(supports$upper[i])){
    links%<>%c(function(x){exp(x)+1})
    invlinks%<>%c(function(x){log(x-1)})
    acceptTrans%<>%c(function(xold,xnew){return((xnew-1)/(xold-1))})
    skew%<>%c(T)
  # x<b
  # } else if(is.infinite(supports$lower[i]) & !is.infinite(supports$upper[i])){
  #   links%<>%c(function(x){supports$upper[i]-exp(x)})
  #   invlinks%<>%c(function(x){log(supports$upper[i]-x)})
  #   acceptTrans%<>%c(function(xold,xnew){return((supports$upper[i]-xnew)/(supports$upper[i]-xold))})

  # 0<x<1
  } else if(supports$lower[i]==0 & supports$upper[i]==1){
    links%<>%c(function(x){(1*exp(x)+0)/(exp(x)+1)})
    invlinks%<>%c(function(x){log(x-0) - log(1-x)})
    acceptTrans%<>%c(function(xold,xnew){return((1-xnew)/(1-xold)*(xnew-0)/(xold-0))})
    skew%<>%c(T)
    
  # 0<x<0.5
  } else if(supports$lower[i]==0 & supports$upper[i]==0.5){
    links%<>%c(function(x){(0.5*exp(x)+0)/(exp(x)+1)})
    invlinks%<>%c(function(x){log(x-0) - log(0.5-x)})
    acceptTrans%<>%c(function(xold,xnew){return((0.5-xnew)/(0.5-xold)*(xnew-0)/(xold-0))})
    skew%<>%c(T)
    
  } else stop("Check your model parameter support ranges: link function not defined (CodeSkeleton.R)")
  
}

# Calculate the multiplicative factor required in the acceptance probability calculation
modifyAcc<-function(xold,xnew) prod(sapply(1:length(xold),function(i){do.call(acceptTrans[[i]],list(xold=xold[i],xnew=xnew[i]))}))

# The start of the list of functions, parameters and formatting for the logTarget
IPMLTP <- list(
  skeleton = skeleton,
  links = links,
  invlinks = invlinks,
  skew = skew,
  initPrior = initPrior,
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

if(algorithm=="ABCSIR"){
  IPMLTP$sampleState <-vectorisedSamplerIPM_ABCSIR
} else if(algorithm=="AMCMC"){
  IPMLTP$sampleState <-vectorisedSamplerIPM
} else stop("Please choose an algorithm to parameterise the model, default is 'ABCSIR'")

##################### CONVERT VALUES WITH LINK FUNCTIONS  ######################
Sample2Physical<-function(x0,IPMLTP){
  x0%<>%unlist()
  for (i in 1:length(IPMLTP$links))  x0[i] <- IPMLTP$links[[i]](x0[i])
  x0%<>%relist(skeleton=IPMLTP$skeleton)
  if(!is.null(IPMLTP$DTN)) {
    x0$growthPars%<>%c(IPMLTP$DTN$L,IPMLTP$DTN$U)
    x0$offSizePars%<>%c(IPMLTP$DTN$L,IPMLTP$DTN$U)
  }
  return(x0)
}

Physical2Sample<-function(x0,IPMLTP){
  x0%<>%unlist()
  for (i in 1:length(IPMLTP$invlinks))  x0[i] <- IPMLTP$invlinks[[i]](x0[i])
  x0%<>%relist(skeleton=IPMLTP$skeleton)
  if(!is.null(IPMLTP$DTN)) {
    x0$growthPars%<>%c(IPMLTP$DTN$L,IPMLTP$DTN$U)
    x0$offSizePars%<>%c(IPMLTP$DTN$L,IPMLTP$DTN$U)
  }
  return(x0)
}


