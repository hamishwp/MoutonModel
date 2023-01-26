# warning("Please run 'sudo apt-get install libgsl-dev' in terminal (Linux, Mac only)")

list.of.packages <- c("xtable","magrittr","doParallel","Rfast", "tidyverse",
                      "abind")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(new.packages)
if(length(new.packages)>0) install.packages(new.packages)

# if(length(list.of.packages[!("mvtnorm" %in% installed.packages()[,"Package"])])){devtools::install_github("rCarto/osrm")}

source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))

library(dissPackage3)
library(xtable)
# library(mcmcse)
library(tidyverse)
library(magrittr)
library(Rfast)

dir.create("./Results",showWarnings = F)
