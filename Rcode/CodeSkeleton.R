
# Skeleton frame for the parameterisation vector
skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA
)

# Link functions to be used
returnSelf <- function(x) x
linkNum <- function(x) exp(x)+1
if(oneSex) {Schilder <- function(x) 0.5*plogis(x)
} else {Schilder <- function(x) plogis(x)}
links<-c(
  'returnSelf', # Survival Logistic Regression Intercept
  'returnSelf', # Survival Logistic Regression Gradient
  'returnSelf', # Growth Linear Regression Intercept
  'returnSelf', # Growth Linear Regression Gradient
  'exp', # Growth Linear Regression Dispersion (Std. Dev.)
  'returnSelf', # Reproduction Logistic Regression Intercept
  'returnSelf', # Reproduction Logistic Regression Gradient
  'linkNum', # Offspring Number per Birth
  'returnSelf', # Offspring Size Linear Regression Intercept
  'returnSelf', # Offspring Size Linear Regression Gradient
  'exp', # Offspring Size Linear Regression Dispersion (Std. Dev.)
  'Schilder' # Offspring Survival Probability
)

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
}else growthFunc <- offSizeFunc <- normal

# Observed Probability Beta Shape Param 1 & 2
if(!fixedObsProb) IPMLTP$links%<>%c('exp','exp')
# Make sure that the skeleton frame also includes this
if(!fixedObsProb) IPMLTP$skeleton %<>% c(list(obsProbPar = rep(NA,2)))
x0%<>%relist(skeleton = IPMLTP$skeleton)

