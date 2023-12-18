
R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> # directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
> directory<-paste0(getwd(),"/")
> 
> list.of.packages <- c("xtable","magrittr","doParallel","Rfast","mc2d", 
+                       "abind")
> new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
> print(new.packages)
character(0)
> if(length(new.packages)>0) install.packages(new.packages)
> 
> # if(length(list.of.packages[!("mvtnorm" %in% installed.packages()[,"Package"])])){devtools::install_github("rCarto/osrm")}
> 
> source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
Loading required package: foreach
Loading required package: iterators
Loading required package: Rcpp
Loading required package: RcppZiggurat

Attaching package: ‘Rfast’

The following objects are masked from ‘package:mvtnorm’:

    dmvnorm, dmvt, rmvt

> source(paste0(directory,'Rcode/SimulateData.R'))
> source(paste0(directory,'Rcode/ModelSpecIPM.R'))
> source(paste0(directory,'Rcode/piecemealFunctions.R'))
> source(paste0(directory,'Rcode/SMC.R'))

Attaching package: ‘mc2d’

The following objects are masked from ‘package:base’:

    pmax, pmin

> source(paste0(directory,'Rcode/SoaySheepData.R'))
── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
✔ ggplot2 3.3.5     ✔ purrr   0.3.4
✔ tibble  3.1.6     ✔ dplyr   1.0.7
✔ tidyr   1.1.4     ✔ stringr 1.4.0
✔ readr   2.1.1     ✔ forcats 0.5.1
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ purrr::accumulate() masks foreach::accumulate()
✖ tidyr::extract()    masks magrittr::extract()
✖ dplyr::filter()     masks stats::filter()
✖ purrr::is_integer() masks Rfast::is_integer()
✖ dplyr::lag()        masks stats::lag()
✖ dplyr::nth()        masks Rfast::nth()
✖ purrr::set_names()  masks magrittr::set_names()
✖ purrr::transpose()  masks Rfast::transpose()
✖ purrr::when()       masks foreach::when()

Attaching package: ‘pracma’

The following object is masked from ‘package:purrr’:

    cross

The following objects are masked from ‘package:magrittr’:

    and, mod, or

The following objects are masked from ‘package:Rfast’:

    Norm, Rank, squareform

> library(dissPackage3)
> library(xtable)
> # library(mcmcse)
> library(tidyverse)
> library(magrittr)
> 
> namer<-"REAL_GSF_beta_multMu_MultObs_GLMx0_final_5brks_autoshift"
> 
> # Is the population counted one sex or two?
> oneSex<-T
> # Is the observation probability an empirically-based fixed value or sampled as a R.V.?
> fixedObsProb<-T
> # Read in the Soay sheep data
> lSHEEP<-GetSoaySheep(directory,oneSex=oneSex)
> # Number of MCMC simulations
> itermax <- 20000
> # Number of SMC particles
> SMC_parts<-6700  # 6700 # 200
> # Number of in-chain parallelised cores
> ncores<-4
> 
> # For the individual and offspring growth function - normal or truncated normal, or otherwise?
> normsampler<-"sampleDTN"
> 
> # Skeleton frame for the parameterisation vector
> skeleton = list(
+   survPars = rep(NA, 2),
+   growthPars = rep(NA, 3),
+   reprPars = rep(NA, 2),
+   offNumPars = NA,
+   offSizePars = rep(NA, 3),
+   Schild = NA
+ )
> 
> 
> # Link functions to be used
> returnSelf <- function(x) x
> linkNum <- function(x) exp(x)+1
> if(oneSex) {Schilder <- function(x) 0.5*plogis(x)
+ } else {Schilder <- function(x) plogis(x)}
> links<-c(
+   'returnSelf', # Survival Logistic Regression Intercept
+   'returnSelf', # Survival Logistic Regression Gradient
+   'returnSelf', # Growth Linear Regression Intercept
+   'returnSelf', # Growth Linear Regression Gradient
+   'exp', # Growth Linear Regression Dispersion (Std. Dev.)
+   'returnSelf', # Reproduction Logistic Regression Intercept
+   'returnSelf', # Reproduction Logistic Regression Gradient
+   'linkNum', # Offspring Number per Birth
+   'returnSelf', # Offspring Size Linear Regression Intercept
+   'returnSelf', # Offspring Size Linear Regression Gradient
+   'exp', # Offspring Size Linear Regression Dispersion (Std. Dev.)
+   'Schilder' # Offspring Survival Probability
+ )
> 
> # The start of the list of functions, parameters and formatting for the logTarget
> IPMLTP <- list(
+   skeleton = skeleton,
+   links = links,
+   survFunc = match.fun('linLogit'),
+   growthSamp = match.fun(normsampler),
+   reprFunc = match.fun('linLogit'),
+   offNumSamp = match.fun('PoisNum'),
+   offSizeSamp = match.fun(normsampler)
+ )
> 
> if(normsampler=="sampleDTN") {
+   IPMLTP$growthFunc <- IPMLTP$offSizeFunc <- doublyTruncatedNormal
+   L<-min(c(lSHEEP$solveDF$prev.size, lSHEEP$solveDF$size),na.rm = T)
+   U<-max(c(lSHEEP$solveDF$prev.size, lSHEEP$solveDF$size),na.rm = T)
+ }else growthFunc <- offSizeFunc <- normal
> 
> # Make the inital values using the piecemeal GLM approach
> x0<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb)))
[1] "Survival: -7.15806203640638, 2.90891568175715"
[1] "Growth: 1.49111882577943, 0.530693641638976, 0.0880691780935396"
[1] "Reproduction: -4.61944325120285, 1.36969697430432"
[1] "Offspring number: 1.06720430107527"
[1] "Offspring size: 0.688701918947741, 0.593104183012289, 0.207365851083968"
[1] "Offspring survival: 0.997819433517082"
[1] "Mean observed probability: 0.499194384802533"
> # x0<-readRDS(paste0(directory,"/Results/x0_GSF_fixed_multMu_multObs_GLMx0"))
> # x0<-readRDS("./Results/x0_GSF_fixed_multMu_multObs_GLMx0_nbks15_autoshift")
> # x0<-c(tx0,x0[(length(x0)-1):length(x0)])
> 
> # Number of parameters
> Np<-length(x0)
> 
> # Import covariance matrix:
> # propCOV<-readRDS("./Results/propCOV_GSF_fixed_multMu_multObs_GLMx0_nbrks15_autoshift")
> # propCOV<-matrix(0,nrow = length(x0),ncol = length(x0))
> # diag(propCOV)<-length(x0)/1000
> # propCOV[1:nrow(tpropCOV),1:nrow(tpropCOV)]<-tpropCOV
> propCOV<-diag(Np)/60
> # diag(propCOV)[13:14]<-c(1e-3,1e-3)
> # propCOV<-readRDS(paste0(directory,"/Results/propCOV_GSF_fixed_multMu_multObs_GLMx0"))
> # propCOV<-readRDS(paste0(directory,"/Results/propCOV_fixedObsP_MH_GSF_poisO_Mult"))
> 
> # Observed Probability Beta Shape Param 1 & 2
> if(!fixedObsProb) IPMLTP$links%<>%c('exp','exp')
> # Make sure that the skeleton frame also includes this
> if(!fixedObsProb) IPMLTP$skeleton %<>% c(list(obsProbPar = rep(NA,2)))
> x0%<>%relist(skeleton = IPMLTP$skeleton)
> # Define the number of size class bins
> nbks<-5
> # Calculate the shift factor that offsets the size class mid-point
> shift<-CalcShift_Kernel(x0,IPMLTP,nbks,oneSex,L,U)
> stop(shift)
Error: 0.247603159213544
Execution halted
