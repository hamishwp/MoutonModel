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
# Load the packages required for the code, and install them if they don't already exist
source(paste0(directory,'Rcode/GetPackages.R'))
# Load the simulation parameters required for this run
source(paste0(directory,'Rcode/LoadParameters.R'))
# Load the IPM object skeleton
source(paste0(directory,'Rcode/CodeSkeleton.R'))
# Either we read in the real data or generate our own!
if(simulation) {source(paste0(directory,'Rcode/SimulateSheepData.R'))
} else source(paste0(directory,'Rcode/SoaySheepData.R'))
# Build up the model based on the parameters and sheep data
source(paste0(directory,'Rcode/BuildModel.R'))
# Check that the log-likelihood values make sense before running the full parallelised code
print("Initial Values Log-Likelihood=")
print(logTargetIPM(x0, logTargetPars = IPMLTP, returnNeg = T, printProp = F))
# Save everything we need to replicate this run:
earlytag<-paste0(namer,"_",priorName,"_its",itermax,"_",
                 gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = ""),"_rand",round(runif(1,max = 1000)))
saveRDS(list(
  x0=x0,
  propCOV=propCOV,
  IPMLTP=IPMLTP
), paste0(directory,"Results/INPUT_",earlytag))
###################################################################################
######################## Parameterise the model using MCMC ########################
###################################################################################
print("And so it begins...")
ptm <- proc.time()
Sheepies <- pM_GaA(propCOV = propCOV,
                         lTarg = logTargetIPM, lTargPars = IPMLTP,
                         cores = ncores,
                         x0 = x0,
                         itermax=itermax,
                         Params=list(GreedyStart=500,Pstar=0.234, gamzy0=NA, epsilon=1, minVar=1e-9),
                         clNeeds = clNeeds, packages = "dissPackage3", prntPars=TRUE)
ptm_fin<-proc.time() - ptm;
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

# Short-term to do
# - parallelise properly obsFun when funcys is only one function
# - check Sstar dimensionality is consistent across all possible models
# - instead of mu initialisation for only total number of individuals,
#   shouldn't we generate all summary statistics?
# - check if sw is iteratively-additive/multiplicative
# - check the output dimension order for output$d_i and output$shat
# - code-up algorithm to evaluate number of SMC particles required
# - sort outshell thing
# - make it possible to choose between Minkowski or distribution-based obsFun
# - median summary stats to be output to the ABCSMC algorithm
# - what is NoParts doing in all the old obsProb models?
# - is output$d_i even needed?


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




