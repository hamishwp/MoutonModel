# warning("Please run 'sudo apt-get install libgsl-dev' in terminal (Linux, Mac only)")

list.of.packages <- c("xtable","magrittr","doParallel","Rfast", "tidyverse",
                      "abind","R.utils","sn", "boot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)>0) {
  print(new.packages)
  install.packages(new.packages)
}

if(length(list.of.packages[!("densratio" %in% installed.packages()[,"Package"])])){devtools::install_github("hoxo-m/densratio")}

source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))
source(paste0(directory,'Rcode/ObsDistance.R'))

library(dissPackage3)
library(xtable)
# library(mcmcse)
library(tidyverse)
library(magrittr)
library(Rfast)
library(densratio)

dir.create("./Results",showWarnings = F)
