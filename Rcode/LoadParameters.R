###################################################################################
################# Load the parameters required for the simulation #################
###################################################################################
# Are we using the real Soay sheep data or are we simulating it?
simulation<-F
# SIMULATION PARAMETERS - UNUSED IF simulation = F
poptot<-100 # The number of years for the simulation (IF SIMULATED)
yearing<-10 # The total population (IF SIMULATED)
# Is the population counted one sex or two?
oneSex<-T
# Is the observation probability an empirically-based fixed value or sampled as a R.V.?
fixedObsProb<-T
# Number of MCMC simulations
itermax <- 30000
# Do we need to calculate the number of particles required for the 
calcParts<-F
# Number of in-chain parallelised cores
ncores<-14
# Define the number of size class bins
nbks<-10
# Bins method:
regbinspace<-F
# Particle filter initialisation function
muModel<-'poisson' #'multinomial'
# Observation Model
obsModel<-'poisson' #'multinomial' #'binomial' 
# Do we automatically calculate the shift in the staggered grid of the size-class bins, based on the IPM kernal?
manshift<-F
# For the individual and offspring growth function - normal or truncated normal, or otherwise?
normsampler<-"sampleDTN"
# What algorithm to use to parameterise the model?
algorithm<-"ABCSIR"
# Which perturbation function to the aSMC resampler?
perturber<-"pert_GlobCov"
# If using a nearest neighbour perturbation, how many neighbours are required?
pNNs<-50
# Use these parameters to create a name for the output file from the simjulation
namer<-paste0(ifelse(simulation,paste0("SIM_pop",poptot,"_yr",yearing),"REAL"),
              "_",algorithm,"_",perturber,"_",ifelse(fixedObsProb,"fixed","beta"),"_",muModel,"Mu_",obsModel,
              "Obs_GLMx0_",itermax,"_",nbks,"brks_","regbinspace",regbinspace,"_",normsampler,"_",ifelse(manshift,"manshift","autoshift"))

###################################################################################
###################################################################################
###################################################################################