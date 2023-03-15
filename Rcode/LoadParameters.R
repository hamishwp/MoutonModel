###################################################################################
################# Load the parameters required for the simulation #################
###################################################################################
# Are we using the real Soay sheep data or are we simulating it?
simulation<-T
if(simulation){
  poptot<-100 # The number of years for the simulation (IF SIMULATED)
  yearing<-10 # The total population (IF SIMULATED)
}
# Is the population counted one sex or two?
oneSex<-T
# Is the observation probability an empirically-based fixed value or sampled as a R.V.?
fixedObsProb<-T
# Number of MCMC simulations
itermax <- 60000
stepmax <- 10
ABCNP<-1000L  # this is the number of particles to pass the ABC threshold
ABCk<-1.2 # this sets the number of particles to trial in ABC as N_trial=k*N (see table 2, U. Simola, et al, Bayesian Analysis (2021) 16, Number 2, Adaptive Approximate Bayesian Computation
# Do we need to calculate the minimum number of particles required for the adaptive-epsilon algorithm?
calcParts<-F
# Number of in-chain parallelised cores
ncores<-50
# Define the number of size class bins
nbks<-10
# Bins method:
regbinspace<-F
# Particle filter initialisation function
muModel<-'poisson' #'multinomial'
# Observation Model
obsModel<-"MADadaptdist"
# Adaptive ABC threshold method
DeltaCalc<-"QuantESS"
# How wide do we want to initialise our proposal distribution?
InitSD<-3
# Do we narrow things down with high level priors?
HLPon<-F
# Do we automatically calculate the shift in the staggered grid of the size-class bins, based on the IPM kernal?
manshift<-F
# For the individual and offspring growth function - normal or truncated normal, or otherwise?
normsampler<-"sampleDTN"
# Proposal distribution
PropDist<-"MVSN"
# What algorithm to use to parameterise the model?
algorithm<-"ABCSIR"
# For the ABCSMC particle weights, do we want the standard or the alternative (Filipi 2012) method?
altWeights<-T
# Which perturbation function to the aSMC resampler?
perturber<-"pert_GlobSkewCov"
# If using a nearest neighbour perturbation, how many neighbours are required?
pNNs<-50
# Check the minimum number of ABCSMC particles for the adaptive epsilon threshold function
NpCheck<-F
# Use these parameters to create a name for the output file from the simjulation
namer<-paste0(ifelse(simulation,paste0("SIM_pop",poptot,"_yr",yearing),"REAL"),
              "_",algorithm,"_",perturber,"_",ifelse(fixedObsProb,"fixed","beta"),"_",muModel,"Mu_",obsModel,
              "Obs_",itermax,"_",nbks,"brks_","regbinspace",regbinspace,"_",normsampler,"_",ifelse(manshift,"manshift","autoshift"),"_rand",round(runif(1,max = 1000)))

###################################################################################
###################################################################################
###################################################################################