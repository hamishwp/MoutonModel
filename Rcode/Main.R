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
