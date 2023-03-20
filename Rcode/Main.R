#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@                                                                           @@@#
#@@@                      Soay Sheep Bayesian-IPM Model                        @@@#
#@@@                                                                           @@@#
#@@@                             Hamish Patten                                 @@@#
#@@@                      Postdoc - Biodemography Group                        @@@#
#@@@                        Department of Statistics                           @@@#
#@@@                          University of Oxford                             @@@#
#@@@                                                                           @@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
directory<-paste0(getwd(),"/")
# Load the simulation parameters required for this run
source(paste0(directory,'Rcode/LoadParameters.R'))
# Load the packages required for the code, and install them if they don't already exist
source(paste0(directory,'Rcode/GetPackages.R'))
# Load the IPM object skeleton
source(paste0(directory,'Rcode/CodeSkeleton.R'))
# Either we read in the real data or generate our own!
if(simulation) {source(paste0(directory,'Rcode/SimulateSheepData.R'))
} else source(paste0(directory,'Rcode/SoaySheepData.R'))
# Build up the model based on the parameters and sheep data
source(paste0(directory,'Rcode/BuildModel.R'))
# Check that the log-likelihood values make sense before running the full parallelised code
print("Initial Values Log-Likelihood=")
ptm <- proc.time()[3]
tmp<-logTargetIPM(initSIR$x0, logTargetPars = IPMLTP, returnNeg = F, printProp = F)
print(paste0("Distance at GLM: ",tmp$distance))
print(paste0("Median SS: ",median(tmp$shat)))
ptm_fin<-(proc.time()[3] - ptm); timeouter<-ptm_fin*3; print(paste0("Timeout = ",timeouter))
# Ensure that the minimum delta distance ensures that all-zero summary stats are not included
initSIR%<>%CalcMinDelta(IPMLTP)
# Make sure this timeout is integrated into the simulations
initSIR$timeouter<-timeouter
# In case user wants to start from a previous simulation
prevSim<-NULL
# prevSim<-"./Results/output_SIM_pop100_yr10_ABCSIR_pert_GlobSkewCov_fixed_poissonMu_FudgerObs_60000_10brks_regbinspaceFALSE_sampleDTN_autoshift_rand142000"
# if(!is.null(prevSim) & file.exists(paste0("output_",namer))) {
#   output<-readRDS(prevSim)
#   exty<-str_split(prevSim,"output_")[[1]][2]
#   inexty<-grep(exty,list.files("./Results/"),value = T)
#   inpy<-readRDS(paste0("./Results/",inexty))
#   IPMLTP<-inpy$IPMLTP
#   initSIR<-inpy$initSIR
#   initSIR$output<-initSIR$output[[length(initSIR$output)]]
# }
# print out the initial proposal distribution parameters
print(initSIR$x0)
print(diag(initSIR$propCOV))
# Save everything we need to replicate this run:
earlytag<-paste0(namer,"_",priorName,"_its",itermax,"_",gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = ""))
saveRDS(list(
  initSIR=initSIR,
  IPMLTP=IPMLTP
), paste0(directory,"Results/INPUTS/INPUT_",earlytag)); IPMLTP$solveDF<-IPMLTP$detectedNum<-NULL
###################################################################################
####################### Parameterise the model using ABCSMC #######################
###################################################################################
print("And so it begins...")
ptm <- proc.time()[3]

Sheepies<-ABCSIR(initSIR, lTarg = logTargetIPM, lTargPars = IPMLTP)

ptm_fin<-proc.time()[3] - ptm;
print(paste0("ncpus= ",ncores," : ",ptm_fin))
###################################################################################
################################ Save the output! #################################
###################################################################################
# Create a unique tag identifier for this specific simulation
tag<-paste0(namer,"_",priorName,"_its",itermax,"_",
            gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = ""))
# Save chains object
saveRDS(Sheepies, paste0(directory,"Results/",tag))
###################################################################################
###################################################################################
###################################################################################

# TO DO MARCH:
# 1) Start with MultiSD=3
# 2) Try the different distance measures
# 3) Extend to a larger number of years
# 4) Use Supremum instead of QuantESS



# Short-term to do
# - Sort out sampling from the priors?
# - Add different resample & perturbation functions
# - what is NoParts doing in all the old obsProb models?
# - Sample from priors for initialisation: 
#   if no prior values are provided, setup the priors distributions from the GLM
# - Do we try using as many threshold values as there are summary statistics?
#   by doing the local-distance weighting in NN MVN models we could use this as 
#   a way of justifying not having to use k-epsilon thresholds which is more computationally expensive
# - Encode the ESS and include a threshold of N/2 lower than which the code spits out a warning (or stops entirely)
#   this ESS threshold comes from Moral, et al, 2011 -  10.1007/s11222-011-9271-y
# - Implement Euclidean distance function which is scaled in parameter space (by the range of values or variance?)
#   that will be used to find the M nearest neighbours


# A tutorial introduction to Bayesian inference for stochastic epidemic
# models using Approximate Bayesian Computation.
# T. Kypraios, 2016, mathematical biosciences
# 10.1016/j.mbs.2016.07.001

# ABCSMC with high summary statistic dimensionality (Kypraios, 2016):
# 'The problem is that simulations which produce good matches
# of all summaries simultaneously become increasingly unlikely as p
# grows. A formal treatment of the issue is given by Blum [9], Barber
# et al. [3] and Biau et al. [8] amongst others.'

# Adaptive ABC threshold schemes (Kypraios, 2016):
# 'We found that adaptive tolerances sometimes perform poorly
# for discrete summary statistics, which several of our applications
# use. The problem occurs when most of the distances from accepted
# simulations exactly equal the tolerance epsilon_t.' 
# What happens is that epsilon_t = epsilon_t-1 and the algorithm gets stuck.
# This could be navigated by having monotonically decreasing thresholds

# Union metric - one single ABC threshold that requires a distance metric of the summary statistics 
#                in order to accept or reject the particle
# Intersection metric - either individual thresholds are produced, per summary statistic,
#                       or every summary statistic needs to pass a single threshold individually
#                       making the threshold max(S_h)<epsilon instead
# We could create and use three types here: union, max-intersection and full-intersection
# whereby the full-intersection method would be created by quantiles of each individual summary statistic
# just as the union method currently does for us!

# Approximate Bayesian Computation and Simulation-Based Inference for Complex Stochastic Epidemic Models
# T.J. McKinley, 2018, 10.1214/17-STS618
# 'Combining metrics in a sensible manner is sometimes challenging,
# especially if they are defined on different scales (see 
# e.g., Conlan et al., 2012). Union metrics can sometimes
# lead to simulations being regularly accepted when they
# do not fit certain outputs very well at all, whereas in-
#   tersection metrics can penalise misfitting simulations
# more, but at a cost of reduced acceptance rates.'

# PERTURBATION KERNELS
# need to distinguish between guided-centered normal distributions and non-guided
# - Local-distance mean- & covariance-weighted kernel using multivariate normal distribution centered and scaled by M-nearest neighbours
# - Local-distance mean-weighted only kernel using multivariate normal distribution centered and scaled by M-nearest neighbours
# - Global-distance weighted kernel using multivariate normal distribution centered and scaled by M-nearest neighbours
# - Un-weighted kernel using multivariate normal distribution centered and scaled by M-nearest neighbours
# - Global multivariate normal distrbution
# - Global univariate (component-wise) normal distribution 
# - Multivariate normal kernel with optimal local covariance matrix
# - Guided (by weighted - by distance - average of theta* and mean(theta_t-1)) multivariate normal kernel with optimal local covariance matrix


# PLAN OF TO-DO
# 1) modify logTarget function to output objective function, per ABC particle
# 2) develop the perturbation function (fullcond algorithm)
# 3) particle weight recalculation based on the objective function summary stats
#   DONE!!!!


# Issues to resolve:
# - Initial acceptance threshold (DONE)
# - Adaptive threshold method (DONE)
# - Distance metric (DONE)
# - Resampling & perturbation method (NEEDS IMPLEMENTING)
# - Stopping criteria (DONE)

# Code up a new particle_filter, including providing outputting the summary statistics of each particle at each time step
# Finish the ABCSMC algorithm!





# proposed<-Sample2Physical(x0,IPMLTP)
# stateSpaceSampArgs <- list(survFunc = IPMLTP$survFunc, survPars = proposed$survPars,
#                            growthSamp = IPMLTP$growthSamp, growthPars = proposed$growthPars,
#                            reprFunc = IPMLTP$reprFunc, reprPars = proposed$reprPars,
#                            offNumSamp = IPMLTP$offNumSamp, offNumPars = proposed$offNumPars,
#                            offSizeSamp = IPMLTP$offSizeSamp, breaks = IPMLTP$breaks,
#                            offSizePars = proposed$offSizePars, Schild=proposed$Schild,
#                            sizes=IPMLTP$sizes, oneSex = IPMLTP$oneSex)
# sampleState <- vectorisedSamplerIPM_ABCSIR







# NOTES on ABCSMC
# MAX:
# Initialisation
#   setup empty vectors, ABC thresholds, particle weights
#   sample particles from initial prior distributions
#   calculate distances for initial particles
# Enter step-loop
#   calculate new ABC tolerance based on the old one,
#   (using alpha proportion target of particles that should make it through based on desired value ~0.9) 
#   recalculate particle weights based on meeting the ABC threshold
#   weighted resample of particles (if necessary... see line 818 Method.R of ODDRIN)
#   note that Np/Ess = Np/sum(weights_i^2)
#   add the perturbation to resampled particles
#     update propCOV of perturbation function
#     (using weights of the previous-step particles that would also meet current threshold)
#   calculate distances for each perturbed particle
#   calculate MCMC acceptance probability 
#   (remembering to adjust acceptance ratios based on link functions - modifyAcc function)
#   if passes acceptance test, perturbed particle is accepted, otherwise keep original particle


# Initialisation
#   setup empty vectors, ABC thresholds, particle weights
#   sample particles from initial prior distributions
#   calculate distances for initial particles
# Enter step-loop


# Differences with Clement/Zhixiao work
#   we define an adaptive ABC-threshold, based on culling a certain number of particles
#     this is slow and requires many steps
#   they pre-define the number of particles to make it through and reduce the threshold accordingly
#     this leads to restricting parameter space to only leave the particles with smallest distances
#   do they still add a perturbation component?
#   does Max keep going until Np particles are simulated?


# Make sure I am not confusing the quartile-reduction method (using alpha - see Max)
# with something else - should I be continuously running until all Np particles are accepted?
# 




# Sheepies <- pM_GaA(propCOV = propCOV,
#                          lTarg = logTargetIPM, lTargPars = IPMLTP,
#                          cores = ncores,
#                          x0 = x0,
#                          itermax=itermax,
#                          Params=list(GreedyStart=500,Pstar=0.234, gamzy0=NA, epsilon=1, minVar=1e-9),
#                          clNeeds = clNeeds, packages = "dissPackage3", prntPars=TRUE)