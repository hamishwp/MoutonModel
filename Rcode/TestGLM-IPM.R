# directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
directory<-paste0(getwd(),"/")
setwd(directory)
source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))
library(dissPackage3,lib.loc=directory)
library(xtable)
library(mcmcse)
library(tidyverse)
library(magrittr)

startValues <- list(
  survPars = c(-9.65, 3.77),
  growthPars = c(1.41, 0.56, log(0.08)),
  reprPars = c(-7.23, 2.6),
  offNumPars = 1,
  offSizePars = c(0.36, 0.71, log(0.16)),
  Schild = qlogis(0.873),
  obsProbPar = 0.5 # not too close to 50 since this will hurt the chain
)

# Link functions to be used
# returnSelf <- function(x) x
# linkNum <- function(x) exp(x)+1
# if(oneSex) {Schilder <- function(x) 0.5*plogis(x)
# } else {Schilder <- function(x) plogis(x)}
# links<-c(
#   'returnSelf', # Survival Logistic Regression Intercept
#   'returnSelf', # Survival Logistic Regression Gradient
#   'returnSelf', # Growth Linear Regression Intercept
#   'returnSelf', # Growth Linear Regression Gradient
#   'exp', # Growth Linear Regression Dispersion (Std. Dev.)
#   'returnSelf', # Reproduction Logistic Regression Intercept
#   'returnSelf', # Reproduction Logistic Regression Gradient
#   'returnSelf', # Offspring Number per Birth
#   'returnSelf', # Offspring Size Linear Regression Intercept
#   'returnSelf', # Offspring Size Linear Regression Gradient
#   'exp', # Offspring Size Linear Regression Dispersion (Std. Dev.)
#   'Schilder', # Offspring Survival Probability
#   'plogis'
# )

plotNames <- c("Survival Intercept", "Survival Gradient", "Growth Intercept", "Growth Gradient", "log(Growth Sigma)", "Reproduction Intercept", "Reproduction Gradient", "Offspring Number",
               "Offspring Size Gradient", "Offspring Size Intercept", "log(Offspring Size Sigma)", "qlogis(Survival Offspring)", "Observed Probability")

poptots<-c(10,30,50,100,300,500)

# lenny<-1:length(unlist(startValues))
lenny<-c(1:13)
########## CALCULATE SOME VALUES BEFORE WE CAN MAKE THE IPMLTP OBJECT ##########
# ind<-13
nsimsI<-50
vMat<-data.frame(
  survG=seq(from=-8.5,to=-12,length.out=nsimsI),
  survI=seq(from=3.6,to=5,length.out=nsimsI),
  growG=seq(from=1,to=1.8,length.out=nsimsI),
  growI=seq(from=0.1,to=0.65,length.out=nsimsI),
  growS=seq(from=-3,to=0,length.out=nsimsI),
  reprG=seq(from=-9,to=-5,length.out=nsimsI),
  reprI=seq(from=1.5,to=3.5,length.out=nsimsI),
  offNum=seq(from=1,to=20,length.out=nsimsI),
  osizeG=seq(from=0.2,to=0.5,length.out=nsimsI),
  osizeI=seq(from=0.5,to=0.9,length.out=nsimsI),
  osizeS=seq(from=-2,to=-1,length.out=nsimsI),
  Schild=seq(from=0.3,to=0.95,length.out=nsimsI),
  obsP=seq(from=0.2,to=0.95,length.out=nsimsI)
)

# rMSE<-array(dim=c(nsimsI,2,length(lenny)))


funcMSE<-function(pp){
  rMSE<-data.frame()   
  k<-0
  for (j in lenny){
    k<-k+1
    
    vals<-c(-9.65, 3.77,
            1.41, 0.56, log(0.08),
            -7.23, 2.6,
            1.01,
            0.36, 0.71, log(0.16),
            0.9,
            1.0
    )
    
    for (i in 1:nsimsI){
      
      vals[j]<-vMat[i,j]
      # print(paste0("Parameter:",plotNames[j],", value: ",vMat[i,j]))
      
      popinc<-tryCatch(eigen(kernelOneVar(m = 500, 
                           survPars = vals[1:2],
                           survFunc = match.fun('linLogit'),
                           repFunc = match.fun('linLogit'),
                           growthFunc = doublyTruncatedNormal,
                           offSizeFunc = doublyTruncatedNormal,
                           growthPars = c(vals[3:4], exp(vals[5]), 1.5, 3.55),
                           repPars = vals[6:7],
                           offNum=round(vals[8]),
                           offSizePars = c(vals[9:10], exp(vals[11]), 1.5, 3.55),
                           childSurv=vals[12],
                           L = 1.5, U = 3.55,
                           shift=0.5,halfPop = T)), error = function(e) NA)
      if(is.na(sum(popinc$values))) {print(paste0("fucked it, pp=",pp,", ",plotNames[j],"i=",i)); next}
      startWeight<-sample(seq(1.5,3.55,length.out=500),size=pp,replace = T,prob = rowsums(abs(popinc$vectors)))
      
      ############################ GENERATE THE PARAMETERS #################################
      simPars <- list(n=pp, t=10,
                      # set survival details:
                      survFunc = linLogit, survPars = vals[1:2],
                      # set growth details:
                      growthSamp = sampleDTN,
                      growthPars = c(vals[3:4], exp(vals[5]), 1.5, 3.55),
                      # set reproduction details:
                      reprFunc = linLogit, reprPars = vals[6:7],
                      # set offspring number and size distribution details:
                      offNumSamp=PoisNum, offNumPars=vals[8],
                      offSizeSamp = sampleDTN,
                      offSizePars = c(vals[9:10], exp(vals[11]), 1.5, 3.55),
                      # Child survival probability:
                      Schild=vals[12], obsProb=vals[13],
                      # set other miscelaneous parameters:
                      Start = startWeight, thresh=300000, OneGend = TRUE,
                      popPrint = F, verbose=F)
      
      
      ############################ GENERATE THE DATA #################################
      # Get the observations:
      set.seed(102010)
      simmedData <- do.call(simulateIBM, simPars)
      
      x0<-tryCatch(getInitialValues(simmedData,printPars = F), error = function(e) NA)
      # for (i in 1:length(x0)) x0[i] <- match.fun(links[i])(x0[i])
      
      rMSE%<>%rbind(data.frame(Variable=plotNames[j],Value=vals[j],
                               lRMSE=log(sqrt(sum((vals-x0)*(vals-x0)))),
                               estimate=x0[j],
                               TotalPop=pp))
      
    }
    
  }
  return(rMSE)
}

MSEout<-mclapply(X = poptots,
                        FUN = funcMSE,
                        mc.cores = length(poptots))


p<-ggplot(data=filter(rMSE,TotalPop!=10),aes(x=estimate,y=Value))+geom_point(aes(colour=factor(TotalPop)))+ geom_abline(slope=1)+
  xlab("GLM Estimated Value")+ggtitle(plotNames)+theme(plot.title = element_text(hjust = 0.5))+ylab("True Value")+labs(col="Population Size")
p<-p+facet_wrap( ~ Variable, nrow = 4, scales = "free")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("GLM Piecemeal Fit Error");p
ggsave('GLMErrorPopSize.png', plot=p,path = paste0(directory,'Plots/Hamish'),width = 13,height = 10.)

















################### BIVARIATE CORRELATIONS ###################
ind1<-10
nsimsI<-30
vMat1<-data.frame(
  survI=seq(from=-8.5,to=-12,length.out=nsimsI),
  survG=seq(from=3.6,to=5,length.out=nsimsI),
  growI=seq(from=1,to=1.8,length.out=nsimsI),
  growG=seq(from=0.3,to=0.65,length.out=nsimsI),
  growS=seq(from=-3,to=0,length.out=nsimsI),
  reprI=seq(from=-9,to=-5,length.out=nsimsI),
  reprG=seq(from=1.5,to=3.5,length.out=nsimsI),
  offNum=seq(from=1,to=20,length.out=nsimsI),
  osizeI=seq(from=0.2,to=0.5,length.out=nsimsI),
  osizeG=seq(from=0.5,to=0.9,length.out=nsimsI),
  osizeS=seq(from=-2,to=-1,length.out=nsimsI),
  Schild=seq(from=0.3,to=0.95,length.out=nsimsI),
  obsP=seq(from=0.2,to=0.95,length.out=nsimsI)
)
vMat1<-vMat1[,ind1]

ind2<-12
nsimsJ<-30
vMat2<-data.frame(
  survI=seq(from=-8.5,to=-12,length.out=nsimsJ),
  survG=seq(from=3.6,to=5,length.out=nsimsJ),
  growI=seq(from=1,to=1.8,length.out=nsimsJ),
  growG=seq(from=0.3,to=0.65,length.out=nsimsJ),
  growS=seq(from=-3,to=0,length.out=nsimsJ),
  reprI=seq(from=-9,to=-5,length.out=nsimsJ),
  reprG=seq(from=1.5,to=3.5,length.out=nsimsJ),
  offNum=seq(from=1,to=20,length.out=nsimsJ),
  osizeI=seq(from=0.2,to=0.5,length.out=nsimsJ),
  osizeG=seq(from=0.5,to=0.9,length.out=nsimsJ),
  osizeS=seq(from=-2,to=-1,length.out=nsimsJ),
  Schild=seq(from=0.3,to=0.95,length.out=nsimsJ),
  obsP=seq(from=0.2,to=0.95,length.out=nsimsJ)
)
vMat2<-vMat2[,ind2]

# rMSE2D<-array(dim=c(nsimsI,nsimsJ))
rMSE2D<-data.frame() 

vals<-c(-9.65, 3.77,
        1.41, 0.56, log(0.08),
        -7.23, 2.6,
        1.1,
        0.36, 0.71, log(0.16),
        0.873,
        0.5
)

for (j in 1:nsimsJ){

  vals[ind2]<-vMat2[j]
  
  for (i in 1:nsimsI){
    
    vals[ind1]<-vMat1[i]
    
    ############################ GENERATE THE PARAMETERS #################################
    simPars <- list(n=100, t=5000,
                    # set survival details:
                    survFunc = linLogit, survPars = vals[1:2],
                    # set growth details:
                    growthSamp = sampleDTN,
                    growthPars = c(vals[3:5], 1.5, 3.55),
                    # set reproduction details:
                    reprFunc = linLogit, reprPars = vals[6:7],
                    # set offspring number and size distribution details:
                    offNum=vals[8], offSizeSamp = sampleDTN,
                    offSizePars = c(vals[9:11], 1.5, 3.55),
                    # Child survival probability:
                    Schild=vals[12], obsProb=vals[13],
                    # set other miscelaneous parameters:
                    Start = 2.7, thresh=7000, OneGend = TRUE,
                    popPrint = F, verbose=F)
    
    
    ############################ GENERATE THE DATA #################################
    # Get the observations:
    set.seed(102010)
    simmedData <- do.call(simulateIBM, simPars)
    max.cens <- simmedData$census.number %>% max
    simmedData%<>%filter(census.number>=max(c((max.cens-12L),1L)))
    x0<-tryCatch(getInitialValues(simmedData,printPars = F), error = function(e) NA)

    rMSE2D%<>%rbind(data.frame(Vali=vMat1[i],Valj=vMat2[j],lRMSE=log(sqrt(sum((vals-x0)*(vals-x0))))))
    
  }
  
}

vals<-c(-9.65, 3.77,
        1.41, 0.56, log(0.08),
        -7.23, 2.6,
        1.01,
        0.36, 0.71, log(0.16),
        0.9,
        0.5
)

colnames(rMSE2D)<-c(plotNames[ind1],plotNames[ind2],"log(Root Mean-Squared-Error)")
colnames(rMSE2D)<-c("Var1", "Var2", "lRMSE")
ggplot(rMSE2D,aes(Var1,Var2,lRMSE)) + geom_raster(aes(fill=lRMSE),na.rm = T) +
  xlab(plotNames[ind1]) + ylab(plotNames[ind2]) + ggtitle("log(Root Mean-Squared-Error)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data=data.frame(Var1=vals[ind1],Var2=vals[ind2]),
             mapping=aes(Var1,Var2),colour="red",shape = 4,size=5,stroke = 2)
# image(log(rMSE2))
# 
# b<-apply(rMSE2,1,rev)
# filled.contour(y = 50*(1:100)/nsimsI+1,x= 0.8*(1:100)/nsimsJ,log(b),
#                plot.title = title(main = "log(Root Mean-Squared-Error)",
#                             xlab="Observed Probability",ylab="Number of Offspring"))
