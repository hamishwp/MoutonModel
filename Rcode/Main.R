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




# to do:
# - Auto-define link functions & asymmetric acceptance probability based on upper and lower parameter limits
# - auto-calculate the number of particles required in algorithm based on logTarget
# - re-calculate number of particles every 250 iterations?
# - use CI of initial values to build initial covariance function




directory<-paste0(getwd(),"/")
# Load the packages required for the code, and install them if they don't already exist
source(paste0(directory,'Rcode/GetPackages.R'))
# Load the simulation parameters required for this run
source(paste0(directory,'Rcode/LoadParameters.R'))
# Load the IPM object skeleton
source(paste0(directory,'Rcode/CodeSkeleton.R'))
# Load the data and wrangle it appropriately
source(paste0(directory,'Rcode/SoaySheepData.R'))
# If we want simulated data, simulate the data into the lSHEEP object
if(simulation) source(paste0(directory,'Rcode/SimulateSheepData.R'))
# Build up the model based on the parameters and sheep data
source(paste0(directory,'Rcode/BuildModel.R'))
# Check that the log-likelihood values make sense before running the full parallelised code
print("Initial Values Log-Likelihood=")
print(logTargetIPM(x0, logTargetPars = IPMLTP, returnNeg = T, printProp = F))
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
