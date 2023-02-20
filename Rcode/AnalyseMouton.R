# directory<-"/home/patten/Documents/Coding/Oxford/MoutonModel/"
directory<-paste0(getwd(),"/")

source(paste0(directory,'Rcode/AdaptivePMCMC.R'))
source(paste0(directory,'Rcode/SimulateData.R'))
source(paste0(directory,'Rcode/ModelSpecIPM.R'))
source(paste0(directory,'Rcode/piecemealFunctions.R'))
source(paste0(directory,'Rcode/SMC.R'))
# library(dissPackage3,lib.loc=directory)
library(xtable)
# library(mcmcse)
library(magrittr)
library(tidyverse)
# install.packages(c("rcorr","corrplot"))
# devtools::install_github('harrelfe/Hmisc')
# library(Hmisc)
# library(corrplot)

x0<-vals<-c(-7.25, 3.77,
            1.49111883, 0.53069364, log(0.08806918),
            -4.619443,  1.369697,
            log(0.067204),
            0.6887019, 0.5931042, log(0.2073659),
            IPMLTP$invlinks[[12]](0.873),
            log(50),
            log(10))
Np<-length(unlist(x0))

# Convert to physical coordinates
vals%<>%Sample2Physical(IPMLTP)

if(fixedObsProb) {
  # Calculate the expected observation probability
  obsMean<-exp(x0[13])/(exp(x0[13])+exp(x0[14]))
  # Remove the last two parameters
  x0<-x0[1:12]; Np<-length(unlist(x0))
  # Temporal observation probability must be kept fixed
  lSHEEP$obsProbTime <- rep(obsMean,yearing+1)
} else lSHEEP$obsProbTime <- rbeta(yearing+1,vals$obsProbPar[1],vals$obsProbPar[2])

# Convert to physical coordinates
names(x0)<-names(unlist(vals))[-c(6,7,14,15)]
x0true<-x0

details = file.info(list.files("./","output_SIM"))
details = details[with(details, order(as.POSIXct(mtime),decreasing = T)), ]
files = rownames(details); rm(details)
# Take the most recent, by default
exty<-str_split(files[1],"output_")[[1]][2]
output<-readRDS(paste0("./output_",exty))
inexty<-grep(exty,list.files("./Results/"),value = T)
inpy<-readRDS(paste0("./Results/",inexty))
IPMLTP<-inpy$IPMLTP
initSIR<-inpy$initSIR

lSHEEP<-IPMLTP
lSHEEP$detectedNum<-lSHEEP$solveDF%>%filter(survived==1)%>%group_by(census.number)%>%
  summarise(detectedNum=length(size),.groups = 'drop_last')%>%pull(detectedNum)%>%unname()

source(paste0(directory,'Rcode/BuildModel.R'))

# TO CHECK THE INITIAL SAMPLE DISTRIBUTION
sheepies<-cbind(data.frame(Distance=runif(1500)),as.data.frame(PropN$proposal(1500)))
names(sheepies)<-c("Distance",names(x0true))
shdens<-t(apply(sheepies[,-1],1,function(x) unname(unlist(Sample2Physical(x,IPMLTP)))[-c(6,7,14,15)]))
shdens<-cbind(sheepies$Distance,shdens)%>%as.data.frame()
colnames(shdens)<-colnames(sheepies)
shdens%<>%reshape2::melt(id.vars=c("Distance"))
abliny<-data.frame(variable=colnames(sheepies)[-1],Z=unname(unlist(Sample2Physical(x0true,IPMLTP)))[-c(6,7,14,15)])
q<-shdens%>%filter(variable!="Distance")%>%ggplot(aes(value))+geom_histogram(aes(colour=variable,fill=variable),alpha=0.5)+
  geom_vline(data = abliny, aes(xintercept = Z),colour="red")
q<-q+facet_wrap(. ~ variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q
# ggsave("InitialProposalDist_VariableDensities.png", plot=q,path = paste0(directory,'Plots/Hamish/'),width = 12,height = 8)

disties<-scoring<-data.frame()
for(i in 1:length(output)){
  # inds<-output[[i]]$distance>output[[i]]$delta[i]
  inds<-rep(T,length(output[[i]]$distance))
  sheepies<-as.data.frame(cbind(log(-output[[i]]$distance[inds]),output[[i]]$theta[inds,]))
  names(sheepies)<-c("distance",names(x0))
  sheepies<-sheepies[!is.na(sheepies$distance) & !is.infinite(sheepies$distance),]  
  disties%<>%rbind(cbind(reshape2::melt(sheepies,id.vars=c("distance")),data.frame(iteration=rep(i,nrow(sheepies)*(ncol(sheepies)-1)))))
  
  ptmp<-RMSE<-c()
  for(j in 2:ncol(sheepies)) {
    minnie<-min(sum(sheepies[,j]<=x0true[j-1]),sum(sheepies[,j]>=x0true[j-1]))
    ptmp%<>%c(minnie/nrow(sheepies))
    RMSE%<>%c(sqrt(sum((sheepies[,j]-x0true[j-1])^2)))
  }
  scoring%<>%rbind(data.frame(pvalue=ptmp,RMSE=RMSE,variable=names(x0),iteration=i))
} 

scoring%>%group_by(iteration)%>%summarise(ptotal=prod(pvalue),RMSEtot=prod(RMSE))
istep<-scoring%>%group_by(iteration)%>%summarise(ptotal=prod(pvalue),RMSEtot=prod(RMSE))%>%pull(RMSEtot)%>%which.min()
# istep<-length(output)
# istep<-1

scoring%>%ggplot()+geom_point(aes(iteration,RMSE,colour=pvalue))

dimmie<-c(min(disties$distance),log(-output[[1]]$delta[1]))
disties%>%ggplot(aes(distance,value,group=iteration))+geom_density(aes(colour=iteration,fill=iteration,y=..density..),alpha=0.3)+
  xlim(dimmie)+facet_grid(rows=vars(iteration))

output[[istep]]$q_thresh
output[[istep]]$delta

inds<-output[[istep]]$distance>output[[istep]]$delta[istep]
# inds<-rep(T,length(output[[istep]]$distance))
sheepies<-as.data.frame(cbind(log(-output[[istep]]$distance[inds]),output[[istep]]$theta[inds,]))
names(sheepies)<-c("Distance",names(x0))
sheepies<-sheepies[!is.na(sheepies$Distance) & !is.infinite(sheepies$Distance),]

newSh<-data.frame()
for(i in 1:length(output)){
  newSh%<>%rbind(data.frame(NewD=log(apply(output[[i]]$shat[output[[i]]$distance>output[[i]]$delta[[i]],],1,function(x) sum(abs(x-c(IPMLTP$SumStats))))),
                  CurD=log(-output[[i]]$distance[output[[i]]$distance>output[[i]]$delta[[i]]]),iteration=i))
}
newSh%>%ggplot()+geom_point(aes(NewD,CurD))+facet_grid(rows=vars(iteration))
newSh%>%ggplot()+geom_histogram(aes(NewD))+scale_y_log10()+facet_grid(rows=vars(iteration))
newSh%>%ggplot()+geom_histogram(aes(CurD))+scale_y_log10()+facet_grid(rows=vars(iteration))

min(newSh$NewD)


tmp<-as.matrix(sheepies)
colnames(tmp)<-c("Distance",names(x0))
means <- apply(tmp[,-1], 2, mean,na.rm=T)
medis <- apply(tmp[,-1], 2, median,na.rm=T)
lower <- apply(tmp[,-1], 2, quantile, probs = 0.025,na.rm=T)
upper <- apply(tmp[,-1], 2, quantile, probs = 0.975,na.rm=T)
wmean <- apply(tmp[,-1],2,weighted.mean,w=exp(tmp[,1]),na.rm = T)
# get the MAP parameters:
ind <- tmp[,1] %>% which.max
MAP <-  tmp[ind, -1]

summaries <- data.frame(True=x0true,median = medis, mean = means, GLM=initSIR$x0, 
                        lower = lower, upper = upper, weighted_mean=wmean) #, MLE=MLE)
# summaries <- data.frame(true = simulated, MAP = MAP, mean = means, GLM=GLM,
#                         median = medis, lower = lower, upper = upper, weighted_mean=wmean) #, MLE=MLE)
# summaries <- data.frame(simulated = simulated, MAP = MAP, mean = means,MLE=MLE)  
summaries%<>%cbind(data.frame(inCI=summaries$GLM<summaries$upper & summaries$GLM>summaries$lower))

print(summaries)
xtable(summaries)

sheepies%>%ggplot()+geom_density(aes(Distance))

qq<-list()
for(i in 2:ncol(sheepies)){
  shtmp<-sheepies[,c(1,i)]; varname<-names(shtmp)[2]; names(shtmp)[2]<-"Value"
  p<-shtmp%>%ggplot(aes(x=Value,y=Distance))+
    geom_density_2d_filled()+ggtitle(varname)+
    ylab("log(-Distance)") + geom_vline(xintercept = x0true[i-1],colour="red")
    # geom_vline(xintercept = x0[i],colour="red")
  qq<-c(qq,list(p))
}
gridExtra::grid.arrange(qq[[1]],qq[[2]],qq[[3]],qq[[4]],qq[[5]],qq[[6]],
                        qq[[7]],qq[[8]],qq[[9]],qq[[10]],qq[[11]],qq[[12]])

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(colour=..density..), geom="tile", contour = FALSE) +
    scale_fill_gradientn(colours=rainbow(100))
  p
}

GGally::ggpairs(sheepies, lower=list(continuous=my_fn))

# sorters<-sort(sheepies$Distance,index.return=T)
# redSheep<-sheepies[sorters$ix[sorters$x<45],]
# GGally::ggpairs(redSheep, lower=list(continuous=my_fn))

sheepies[which.max(sheepies$Distance),]%>%Sample2Physical(IPMLTP)


cov.wt(sheepies[,-1],sheepies[,1])$center%>%Sample2Physical(IPMLTP)
diag(cov.wt(sheepies[,-1],sheepies[,1])$cov)%>%Sample2Physical(IPMLTP)

meany<-colMeans(sheepies[,-1])%>%Sample2Physical(IPMLTP)
medy<-apply(sheepies[,-1],2,median)%>%Sample2Physical(IPMLTP)
UL<-(cov.wt(sheepies[,-1],sheepies[,1])$center+apply(sheepies[,-1],2,sd))%>%Sample2Physical(IPMLTP)
LL<-(cov.wt(sheepies[,-1],sheepies[,1])$center-apply(sheepies[,-1],2,sd))%>%Sample2Physical(IPMLTP)



sheepies%>%ggplot(aes(x=survPars1,y=survPars2))+
  stat_summary_2d(aes(z=Distance),bins = 40)+stat_density_2d(aes(colour=..level..))+
  geom_hline(yintercept = vals$survPars[2],colour="red")+
  geom_vline(xintercept = vals$survPars[1],colour="red")+
  scale_color_gradient2(mid="red",high="yellow")

sheepies%>%ggplot(aes(x=reprPars1,reprPars2,z=Distance))+
  stat_summary_2d(aes(z=Distance),bins = 40)+stat_density_2d(aes(colour=..level..))+
  geom_hline(yintercept = vals$reprPars[2],colour="red")+
  geom_vline(xintercept = vals$reprPars[1],colour="red")+
  scale_color_gradient2(mid="red",high="yellow")


sheepies%>%ggplot(aes(x=growthPars1,growthPars2,z=Distance))+
  stat_summary_2d(aes(z=Distance),bins = 40)+stat_density_2d(aes(colour=..level..))+
  geom_hline(yintercept = vals$growthPars[2],colour="red")+
  geom_vline(xintercept = vals$growthPars[1],colour="red")+
  scale_color_gradient2(mid="red",high="yellow")

sheepies%>%ggplot(aes(x=growthPars1,growthPars3,z=Distance))+
  stat_summary_2d(aes(z=Distance),bins = 40)+stat_density_2d(aes(colour=..level..))+
  geom_hline(yintercept = log(vals$growthPars[3]),colour="red")+
  geom_vline(xintercept = vals$growthPars[1],colour="red")+
  scale_color_gradient2(mid="red",high="yellow")

sheepies%>%ggplot(aes(x=growthPars2,growthPars3,z=Distance))+
  stat_summary_2d(aes(z=Distance),bins = 40)+stat_density_2d(aes(colour=..level..))+
  geom_hline(yintercept = log(vals$growthPars[3]),colour="red")+
  geom_vline(xintercept = vals$growthPars[2],colour="red")+
  scale_color_gradient2(mid="red",high="yellow")

sheepies%>%ggplot(aes(x=offSizePars1,offSizePars2,z=Distance))+
  stat_summary_2d(aes(z=Distance),bins = 40)+stat_density_2d(aes(colour=..level..))+
  geom_hline(yintercept = vals$offSizePars[2],colour="red")+
  geom_vline(xintercept = vals$offSizePars[1],colour="red")+
  scale_color_gradient2(mid="red",high="yellow")

initDist<-logTargetIPM(x0true, logTargetPars = IPMLTP, returnNeg = F, printProp = F, returnW=T)
initDist$d

hist(log(abs(c(initDist$shat)-c(IPMLTP$SumStats))))
hist(log(abs(output[[istep]]$shat[which.max(sheepies$Distance),]-IPMLTP$SumStats)))

sum(abs(initDist$shat-IPMLTP$SumStats))
sum(abs(output[[istep]]$shat[which.max(sheepies$Distance),]-IPMLTP$SumStats))

hist(log(apply(output[[istep]]$shat,1,function(x) sum(abs(x-IPMLTP$SumStats)))))
log(sum(abs(initDist$shat-IPMLTP$SumStats)))


medSS<-apply(output[[1]]$shat,2,median)
plot(log(apply(output[[istep]]$shat,1,function(x) sum(abs(x-IPMLTP$SumStats)))),log(-output[[istep]]$distance))
apply(output[[istep]]$shat,2,median)
apply(output[[1]]$shat,2,median)

length(medSS)
normMSE<-vapply(1:ncol(output[[istep]]$shat),function(i) abs(output[[istep]]$shat[,i]-IPMLTP$SumStats[i])/medSS[i],FUN.VALUE = numeric(length(output[[istep]]$distance)))
hist(normMSE)
normMSE<-vapply(1:ncol(output[[istep]]$shat),function(i) abs(output[[istep]]$shat[,i]-IPMLTP$SumStats[i])/medSS[i],FUN.VALUE = numeric(length(output[[istep]]$distance)))

sum(abs(output[[istep]]$shat[which.max(sheepies$Distance),]-IPMLTP$SumStats)/normMSE)
sum(abs(initDist$shat-IPMLTP$SumStats)/normMSE)



# Is the observation probability a fixed or random effects model?
# fixedObsProb<-F
# Load the IPM object skeleton
source(paste0(directory,'Rcode/CodeSkeleton.R'))

lSHEEP<-GetSoaySheep(directory,oneSex=T)
GLM<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb)))
for (i in 1:length(GLM)) GLM[i] <- links[[i]](GLM[i])

# make a data.frame, and then tables for the write-up:
# simulated <- c(-9.65, 3.77, 1.41, 0.56, 0.08, -7.23,
#                2.6, 1, 0.36, 0.71, 0.16, 0.873, 1)

# MLE<-readRDS("./Results/MLE_svjitter")$par
# for (i in 1:length(MLE)) MLE[i] <- links[[i]](MLE[i])

# SVjitter<-readRDS(paste0(directory,"RDSobjects/startvals_jitter"))
# for (i in 1:length(SVjitter)) SVjitter[i] <- links[[i]](SVjitter[i])

if(fixedObsProb) {x0<-readRDS(paste0(directory,"RDSobjects/GLM_real"))
} else x0<-readRDS(paste0(directory,"RDSobjects/GLM_real_beta"))
GLM<-x0
for (i in 1:length(GLM)) GLM[i] <- links[[i]](x0[i])

# propCOV<-diag(length(GLM))/60

Sheep1<-readRDS("Results/uniform_its18000_2020-05-08_565")
Sheep2<-readRDS("Results/uniform_its18000_2020-05-10_565")
Sheep3<-readRDS("Results/uniform_its18000_2020-05-11_565")

for (i in 1:4){
  tmp<-get(paste0("Sheep",i))
  
  print(paste0("Max LL: ",max(tmp[,1]),", chain number: ",i))
  
  for (j in 2:ncol(tmp)) tmp[,j] <- links[[j-1]](tmp[,j])
  
  means <- apply(tmp[,-1], 2, mean,na.rm=T)
  medis <- apply(tmp[,-1], 2, median,na.rm=T)
  lower <- apply(tmp[,-1], 2, quantile, probs = 0.025,na.rm=T)
  upper <- apply(tmp[,-1], 2, quantile, probs = 0.975,na.rm=T)
  wmean <- apply(tmp[,-1],2,weighted.mean,w=exp(tmp[,1]),na.rm = T)
  # get the MAP parameters:
  ind <- tmp[,1] %>% which.max
  MAP <-  tmp[ind, -1]

  summaries <- data.frame(GLM=GLM, MAP = MAP, mean = means,
                          lower = lower, upper = upper, weighted_mean=wmean) #, MLE=MLE)
  # summaries <- data.frame(true = simulated, MAP = MAP, mean = means, GLM=GLM,
  #                         median = medis, lower = lower, upper = upper, weighted_mean=wmean) #, MLE=MLE)
  # summaries <- data.frame(simulated = simulated, MAP = MAP, mean = means,MLE=MLE)  
  summaries%<>%cbind(data.frame(inCI=summaries$GLM<summaries$upper & summaries$GLM>summaries$lower))
  
  print(summaries)
  xtable(summaries)
  # print(data.frame(GLM=abs(summaries$simulated-summaries$GLM),
  #                  mean=abs(summaries$simulated-summaries$mean),
  #                  MAP=abs(summaries$simulated-summaries$MAP)))
  
  print(paste0("mean:",round(sqrt(sum((summaries$GLM-summaries$mean)^2)),3)))
  print(paste0("MAP:",round(sqrt(sum((summaries$GLM-summaries$MAP)^2)),3)))
  print(paste0("median:",round(sqrt(sum((summaries$GLM-summaries$weighted_mean)^2)),3)))
  # print(paste0("GLM:",round(sqrt(sum(abs(summaries$simulated-summaries$GLM)^2,na.rm=T)),3)))
  # print(paste0("mean:",round(sqrt(sum((summaries$simulated-summaries$mean)^2)),3)))
  # print(paste0("MAP:",round(sqrt(sum((summaries$simulated-summaries$MAP)^2)),3)))
  # print(paste0("median:",round(sqrt(sum((summaries$simulated-summaries$median)^2)),3)))
  
  pairs(tmp,pch=19,labels = c("LL",names(unlist(GLM))))
  
  # MAPsum<- data.frame(simulated = simulated, GLM=GLM, MultinomObs_PoissonMu=MAP)
  # MEANsum<- data.frame(simulated = simulated, GLM=GLM, MultinomObs_PoissonMu=means)
  # UPsum<- data.frame(simulated = simulated, GLM=GLM, MultinomObs_PoissonMu=upper)
  # LOWsum<- data.frame(simulated = simulated, GLM=GLM, MultinomObs_PoissonMu=lower)
  
  # MAPsum<- cbind(MAPsum,data.frame(PoissonObs_PoissonMu=MAP))
  # MEANsum<- cbind(MEANsum,data.frame(PoissonObs_PoissonMu=means))
  # UPsum<- cbind(UPsum,data.frame(PoissonObs_PoissonMu=upper))
  # LOWsum<- cbind(LOWsum,data.frame(PoissonObs_PoissonMu=lower))
  
  # MAPsum<- cbind(MAPsum,data.frame(DetectionNumObs_MultinomMu=MAP))
  # MEANsum<- cbind(MEANsum,data.frame(DetectionNumObs_MultinomMu=means))
  # UPsum<- cbind(UPsum,data.frame(DetectionNumObs_MultinomMu=upper))
  # LOWsum<- cbind(LOWsum,data.frame(DetectionNumObs_MultinomMu=lower))
  
  # xtable(MAPsum)
  # xtable(MEANsum)
  # xtable(UPsum)
  # xtable(LOWsum)
}

for (i in 1:3){
  tmp<-get(paste0("Sheep",i))
  
  print(paste0("Max LL: ",max(tmp[,1]),", iteration: ",i))
  
  for (i in 2:14) tmp[,i] <- links[[i-1]](tmp[,i])
  
  # get the MAP parameters:
  ind <- tmp[,1] %>% which.max
  MAP <-  tmp[ind, -1]
  
  summaries <- data.frame(simulated = simulated, SVjitter=SVjitter, MLE=MLE, MAP = MAP, Bdiff=MAP-simulated, Mdiff=MLE-simulated)
  print(summaries)
}

combinedChosen<-rbind(Sheep1,Sheep2,Sheep3)

for (i in 2:14) combinedChosen[,i] <- links[[i-1]](combinedChosen[,i])

# get the MAP parameters:
ind <- combinedChosen[,1] %>% which.max
MAP <-  combinedChosen[ind, -1]

# get the median, mean and 95% CIs:
means <- apply(combinedChosen[,-1], 2, mean)
medis <- apply(combinedChosen[,-1], 2, median)
lower <- apply(combinedChosen[,-1], 2, quantile, probs = 0.025)
upper <- apply(combinedChosen[,-1], 2, quantile, probs = 0.975)

summaries <- data.frame(simulated = simulated, SVjitter=SVjitter, MLE=MLE, MAP = MAP, mean = means,
                        median = medis, lower = lower, upper = upper)

xtable(t(summaries)[,c(1:5, 13)])
xtable(t(summaries)[,-c(1:5, 13)])

library(GGally)
tmp%<>%as.data.frame.array()
ggpairs(tmp)


tmp<-rbind(Sheep1,Sheep2,Sheep3,Sheep4)

# library(modi)
for (i in 1:3){
  tmp<-get(paste0("Sheep",i))
  
  print(paste0("Max LL: ",max(tmp[,1]),", chain number: ",i))
  
  for (j in 2:ncol(tmp)) tmp[,j] <- links[[j-1]](tmp[,j])
  
  means <- apply(tmp[,-1], 2, mean,na.rm=T)
  medis <- apply(tmp[,-1], 2, median,na.rm=T)
  lower <- apply(tmp[,-1], 2, quantile, probs = 0.025,na.rm=T)
  upper <- apply(tmp[,-1], 2, quantile, probs = 0.975,na.rm=T)
  wmean <- apply(tmp[,-1],2,weighted.mean,w=exp(tmp[,1]),na.rm = T)
  # get the MAP parameters:
  ind <- tmp[,1] %>% which.max
  MAP <-  tmp[ind, -1]
  
  summaries <- data.frame(GLM=GLM, MAP = MAP, mean = means,
                          median = medis, lower = lower, upper = upper, weighted_mean=wmean) #, MLE=MLE)
  # summaries <- data.frame(simulated = simulated, MAP = MAP, mean = means,MLE=MLE)  
  
  print(summaries)
  
}

InsideCI<-summaries$GLM>=summaries$lower & summaries$GLM<=summaries$upper
summaries%<>%cbind(InsideCI)

print(summaries)

xtable(t(summaries)[,c(1:7)])
xtable(t(summaries)[,-c(1:7)])



Params<-list(Pstar=0.234, gamzy0=1, epsilon=2, minVar=1e-3)
tmp<-Sheep2
LL<-tmp[,1]
lenny<-length(LL)
cores<-4
propMu<- x0
propCOV<-diag(length(GLM))/60
gamzy <- Params$gamzy0/(1:(lenny+cores))^(seq(from=1/(1+Params$epsilon),to=1,length.out=(lenny+cores)))
Lgsf<-0
Pstar<-Params$Pstar
talpha<-c(1,exp(LL[1:lenny-1]-LL[2:lenny]))
alpha<-apply(cbind(rep(1,lenny),talpha),1,min)
Lgsf<-rep(0,lenny)
sCOV<-rep(0,lenny)
dmu<-rep(0,lenny)
for (i in 1:lenny){
  Lgsf[i] <- Lgsf[i] + gamzy[i]*(alpha[i]-Pstar)
  propMu <- propMu + gamzy[i]*(tmp[i,-1] - propMu)
  propCOV <- propCOV + gamzy[i]*((tmp[i,-1] - propMu)%*%t(tmp[i,-1] - propMu) - propCOV)
  sCOV[i]<-sum(propCOV*t(propCOV))
  dmu[i]<-sum((tmp[i,-1] - propMu)%*%t(tmp[i,-1] - propMu))
}

plot(log(sCOV))
plot(dmu)
plot(exp(Lgsf),ylim=c(0.8,1.2))
plot(log(alpha),ylim=c(-5,0))
# propMu
# max(propCOV)

plotNames <- c("Survival Intercept", "Survival Gradient", "Growth Intercept", "Growth Gradient", "log(Growth Sigma)", "Reproduction Intercept", "Reproduction Gradient", "Offspring Number",
               "Offspring Size Gradient", "Offspring Size Intercept", "log(Offspring Size Sigma)", "qlogis(Survival Offspring)", "Observed Probability")
plotNamesEx <-c(plotNames[1:(length(plotNames)-1)],"Obs Prob Shape 1","Obs Prob Shape 2")
  
tmp<-c()
for (i in 1:4) tmp%<>%rbind(get(paste0("Sheep",i)))   #[10000:30000,])
for (i in 2:ncol(tmp)) tmp[,i] <- links[[i-1]](tmp[,i])
tmp<-tmp[tmp[,1]> -275,]

colnames(tmp)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

# library(GGally)
# tmp%<>%as.data.frame.array()
# ggpairs(tmp)
# 
# x0<-colMeans(tmp[,-1])
# propCOV<-cov(tmp[,-1])
# 
# saveRDS(x0,"./Results/x0_GSF_fixed_multMu_multObs_GLMx0")
# saveRDS(propCOV,"./Results/propCOV_GSF_fixed_multMu_multObs_GLMx0")
# 
# cory<-Hmisc::rcorr(tmp[,-1])
# col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot(cory$r, method = "color", col = col(200),  
#          type = "upper", order = "hclust", 
#          addCoef.col = "black", # Add coefficient of correlation
#          tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
#          # Combine with significance level
#          p.mat = cory$P, sig.level = 0.01,  
#          # hide correlation coefficient on the principal diagonal
#          diag = FALSE 
# )

Variable<-Value<-c()
for(i in 1:ncol(tmp)) {
  Variable%<>%c(rep(colnames(tmp)[i],nrow(tmp)))
  Value%<>%c(tmp[,i])
}
shdens<-data.frame(Variable=Variable, Value=Value)

q<-ggplot(shdens,aes(Value,group=Variable))+geom_density(aes(colour=Variable,fill=Variable),alpha=0.5)
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  ggtitle("Poisson Obs with Multinomial Mu - Fixed Observation")+xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) 
ggsave("VariableDensities_MultMuPoisObsFixedAutoshift.png", plot=q,path = paste0(directory,'Plots/Hamish/'),width = 12,height = 8)






MultObs<-tmp[sample(1:nrow(tmp),size = 10000,prob = exp(tmp[,1]-min(tmp[,1])),replace = F),]
PoisObs<-tmp2[sample(1:nrow(tmp2),size = 10000,prob = exp(tmp2[,1]-min(tmp2[,1])),replace = F),]
MultObs[MultObs[,13]<0.99,]<-NA
MultObs[MultObs[,9]>1.15,]<-NA
PoisObs[PoisObs[,13]<0.99,]<-NA
PoisObs[PoisObs[,9]>1.15,]<-NA

PoisObs[,-1]%>%corrr::correlate() %>% corrr::network_plot(min_cor = .2)
MultObs[,-1]%>%corrr::correlate() %>% corrr::network_plot(min_cor = .2)

Variable<-Value<-c()
for(i in 1:ncol(MultObs)) {
  Variable%<>%c(rep(colnames(MultObs)[i],nrow(MultObs)))
  Value%<>%c(MultObs[,i])
}
Model<-rep("Multinomial Obs",nrow(MultObs))
  
for(i in 1:ncol(PoisObs)) {
  Variable%<>%c(rep(colnames(PoisObs)[i],nrow(PoisObs)))
  Value%<>%c(PoisObs[,i])
}
Model%<>%c(rep("Poisson Obs",nrow(PoisObs)))
  
shdens<-data.frame(Variable=Variable, Value=Value, Model=Model)
abliny<-data.frame(Variable=colnames(PoisObs)[-1],Z=unname(unlist(GLM)))

q<-shdens%>%filter(Variable!="Log-Target")%>%ggplot(aes(Value,group=Model))+geom_density(aes(colour=Model,fill=Model),alpha=0.5)+
  geom_vline(data = abliny, aes(xintercept = Z))
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q
ggsave("VariableDensities_MultObs_vs_PoisObs_FixedAutoshift.png", plot=q,path = paste0(directory,'Plots/Hamish/'),width = 12,height = 8)

wilcy<-data.frame(PValue=sapply(2:ncol(MultObs),function(i) {
    wilcy<-shdens%>%filter(Variable==colnames(tmp)[i])
    return(wilcox.test(Value ~ Model,data=wilcy)$p.value)}), Variable=colnames(tmp)[2:ncol(MultObs)])

ttestSameV<-data.frame(PValue=sapply(2:ncol(MultObs),function(i) {
  wilcy<-shdens%>%filter(Variable==colnames(tmp)[i])
  return(t.test(Value ~ Model,data=wilcy, var.equal = TRUE)$p.value)}), Variable=colnames(tmp)[2:ncol(MultObs)])

ttestDiffV<-data.frame(PValue=sapply(2:ncol(MultObs),function(i) {
  wilcy<-shdens%>%filter(Variable==colnames(tmp)[i])
  return(t.test(Value ~ Model,data=wilcy, var.equal = F)$p.value)}), Variable=colnames(tmp)[2:ncol(MultObs)])

# CONCLUSION FROM HYP TESTS - All are s.s. different except Offspring Survival parameter

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Sheep1<-readRDS("Results/REAL_GSF_beta_multMu_poisObs_GLMx0_final_autoshift_uniform_its40000_2021-12-15_231846")
Sheep2<-readRDS("Results/REAL_GSF_beta_multMu_poisObs_GLMx0_final_autoshift_uniform_its40000_2021-12-15_231950")
Sheep3<-readRDS("Results/REAL_GSF_beta_multMu_poisObs_GLMx0_final_autoshift_uniform_its40000_2021-12-15_231955")
Sheep4<-readRDS("Results/REAL_GSF_beta_multMu_poisObs_GLMx0_final_autoshift_uniform_its40000_2021-12-15_232154")

PoisObs<-c()
for (i in 2:4) PoisObs%<>%rbind(get(paste0("Sheep",i))[10000:40000,])
for (i in 2:ncol(PoisObs)) PoisObs[,i] <- links[[i-1]](PoisObs[,i])
colnames(PoisObs)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

cory<-Hmisc::rcorr(PoisObs,type = "spearman")
corrplot::corrplot(cory$r, method = "color", col = col(200),  
                   type = "upper",
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
                   # Combine with significance level
                   p.mat = cory$P, sig.level = 0.01,  
                   # hide correlation coefficient on the principal diagonal
                   diag = FALSE 
)
# ggsave("SpearmanCorrMat_MultMuPoisObsFixedAutoshift.png", plot=q,path = paste0(directory,'Plots/Hamish/'),width = 10,height = 8)

Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_092853")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_100003")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_100738")
Sheep4<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_102617")

MultObs<-c()
for (i in 1:4) MultObs%<>%rbind(get(paste0("Sheep",i))[7000:20000,])
for (i in 2:ncol(MultObs)) MultObs[,i] <- links[[i-1]](MultObs[,i])
colnames(MultObs)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])


cory<-Hmisc::rcorr(MultObs,type = "spearman")
corrplot::corrplot(cory$r, method = "color", col = col(200),  
                   type = "upper",
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
                   # Combine with significance level
                   p.mat = cory$P, sig.level = 0.01,  
                   # hide correlation coefficient on the principal diagonal
                   diag = FALSE 
)


diff<-abs(cory$r-cory2$r)
diff<-2*((diff-min(diff))/(max(diff)-min(diff))-0.5)
corrplot::corrplot(diff, method = "color", col = RColorBrewer::brewer.pal(5, "BrBG"),  
                   type = "upper",
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
                   # Combine with significance level
                   p.mat = cory$P, sig.level = 0.01,  
                   # hide correlation coefficient on the principal diagonal
                   diag = FALSE 
)



















Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_uniform_its30000_2021-12-12_004716")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_uniform_its30000_2021-12-12_004649")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_uniform_its30000_2021-12-12_004525")
Sheep4<-readRDS("./Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_uniform_its30000_2021-12-12_003838")

MultObs15<-c()
for (i in 1:4) MultObs15%<>%rbind(get(paste0("Sheep",i))[7000:30000,])
for (i in 2:ncol(MultObs15)) MultObs10[,i] <- links[[i-1]](MultObs15[,i])
colnames(MultObs15)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

Variable<-Value<-c()
for(i in 1:ncol(MultObs)) {MultObs10
  Variable%<>%c(rep(colnames(MultObs)[i],nrow(MultObs)))
  Value%<>%c(MultObs[,i])
}
Model<-rep("Multinomial Obs - 15 breaks 0.5 shift",length(MultObs))

for(i in 1:ncol(MultObs10)) {
  Variable%<>%c(rep(colnames(MultObs10)[i],nrow(MultObs10)))
  Value%<>%c(MultObs10[,i])
}
Model%<>%c(rep("Multinomial Obs - 10 breaks auto-shift",length(MultObs10)))

shdens<-data.frame(Variable=Variable, Value=Value, Model=Model)
abliny<-data.frame(Variable=colnames(MultObs10)[-1],Z=unname(unlist(GLM)))

q<-shdens%>%filter(Variable!="Log-Target")%>%ggplot(aes(Value,group=Model))+geom_density(aes(colour=Model,fill=Model),alpha=0.5)+
  geom_vline(data = abliny, aes(xintercept = Z))
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q


Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_053728")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_053835")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_054908")
Sheep4<-readRDS("./Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_063434")

MultObs10<-c()
for (i in 1:4) MultObs10%<>%rbind(get(paste0("Sheep",i))[7000:30000,])
for (i in 2:ncol(MultObs10)) MultObs10[,i] <- links[[i-1]](MultObs10[,i])
colnames(MultObs10)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])


Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_20brks_autoshift_uniform_its40000_2021-12-20_133729")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_20brks_autoshift_uniform_its40000_2021-12-20_134011")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_20brks_autoshift_uniform_its40000_2021-12-20_134058")

new20<-c()
for (i in 1:3) new20%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
for (i in 2:ncol(new20)) new20[,i] <- links[[i-1]](new20[,i])
colnames(new20)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])


Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_poisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-15_232154")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_poisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-15_231950")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_poisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-15_231955")

old15<-c()
for (i in 1:3) old15%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
for (i in 2:ncol(old15)) old15[,i] <- links[[i-1]](old15[,i])
colnames(old15)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])


Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_095702")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_095823")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_100914")
Sheep4<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_101019")

old15<-c()
for (i in 1:4) old15%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
for (i in 2:ncol(old15)) old15[,i] <- links[[i-1]](old15[,i])
colnames(old15)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_132955")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_133142")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_133614")
Sheep4<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_133928")
Sheep5<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_095702")
Sheep6<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_095823")
Sheep7<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_100914")
Sheep8<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_101019")

new15<-c()
for (i in 1:8) new15%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
for (i in 2:ncol(new15)) new15[,i] <- links[[i-1]](new15[,i])
colnames(new15)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

Sheep1<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_095626")
Sheep2<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_095906")
Sheep3<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_095948")
Sheep4<-readRDS("./Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_101148")

new10<-c()
for (i in 1:4) new10%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
for (i in 2:ncol(new10)) new10[,i] <- links[[i-1]](new10[,i])
colnames(new10)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])




apply(old15,2,weighted.mean,w=old15[,1])
apply(new15,2,weighted.mean,w=new15[,1])

Variable<-Value<-Model<-c()
for(i in 1:ncol(new10)) {
  Variable%<>%c(rep(colnames(new10)[i],nrow(new10)))
  Value%<>%c(new10[,i])
}
Model%<>%c(rep("Num Classes=10",length(new10)))

for(i in 1:ncol(new15)) {
  Variable%<>%c(rep(colnames(new15)[i],nrow(new15)))
  Value%<>%c(new15[,i])
}
Model%<>%c(rep("Num Classes=15",length(new15)))

for(i in 1:ncol(new20)) {
  Variable%<>%c(rep(colnames(new20)[i],nrow(new20)))
  Value%<>%c(new20[,i])
}
Model%<>%c(rep("Num Classes=20",length(new20)))

shdens<-data.frame(Variable=Variable, Value=Value, Model=Model)
# abliny<-data.frame(Variable=colnames(old15)[-1],Z=unname(unlist(GLM)))

q<-shdens%>%filter(Variable!="Log-Target")%>%ggplot(aes(Value,group=Model))+geom_density(aes(colour=Model,fill=Model),alpha=0.5)
  # geom_vline(data = abliny, aes(xintercept = Z))
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q

apply(new10,2,weighted.mean,w=new10[,1])
apply(new15,2,weighted.mean,w=new15[,1])
apply(new20,2,weighted.mean,w=new20[,1])


Sheep1<-readRDS("Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_092853")
Sheep2<-readRDS("Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_100003")
Sheep3<-readRDS("Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_100738")
Sheep4<-readRDS("Results/REAL_GSF_fixed_multMu_multObs_GLMx0_final_autoshift_uniform_its20000_2021-12-15_102617")
Sheep5<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_5brks_autoshift_uniform_its20000_2021-12-18_052003")
Sheep6<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_053728")
Sheep7<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_053835")
Sheep8<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_054908")
Sheep9<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-18_063434")
Sheep10<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_5brks_autoshift_uniform_its20000_2021-12-19_155657")
Sheep11<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_20brks_autoshift_uniform_its20000_2021-12-19_170551")
Sheep12<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_20brks_autoshift_uniform_its20000_2021-12-19_175303")
Sheep13<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_20brks_autoshift_uniform_its20000_2021-12-19_182403")
Sheep14<-readRDS("Results/REAL_GSF_fixed_multMu_MultObs_GLMx0_final_20brks_autoshift_uniform_its20000_2021-12-19_182905")

multnom<-c()
for (i in 1:14) multnom%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
# for (i in 2:ncol(multnom)) multnom[,i] <- links[[i-1]](multnom[,i])
colnames(multnom)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

multnom<-multnom[!(multnom[,6]< -5),]
multnom<-multnom[sample(1:nrow(multnom),size=10000,replace = F,prob = exp(multnom[,1])),]

Sheep1<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_095626")
Sheep2<-readRDS("Results/REAL_GSF_fixed_multMu_poisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-15_231950")
Sheep3<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_095702")
Sheep4<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_095823")
Sheep5<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_095906")
Sheep6<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_095948")
Sheep7<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_100914")
Sheep8<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_101019")
Sheep9<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_10brks_autoshift_uniform_its40000_2021-12-20_101148")
Sheep10<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_132955")
Sheep11<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_133142")
Sheep12<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_133614")
Sheep13<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_20brks_autoshift_uniform_its40000_2021-12-20_133729")
Sheep14<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-20_133928")
Sheep15<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_20brks_autoshift_uniform_its40000_2021-12-20_134011")
Sheep16<-readRDS("Results/REAL_GSF_fixed_multMu_PoisObs_GLMx0_final_20brks_autoshift_uniform_its40000_2021-12-20_134058")
Sheep17<-readRDS("Results/REAL_GSF_fixed_multMu_poisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-15_231955")
Sheep18<-readRDS("Results/REAL_GSF_fixed_multMu_poisObs_GLMx0_final_15brks_autoshift_uniform_its40000_2021-12-15_232154")

poiss<-c()
for (i in 1:14) poiss%<>%rbind(get(paste0("Sheep",i))) #[7000:30000,])
# for (i in 2:ncol(poiss)) poiss[,i] <- links[[i-1]](poiss[,i])
colnames(poiss)<- c("Log-Target",plotNames[1:(length(plotNames)-1)])

poiss<-poiss[!(poiss[,6]< -5),]
poiss<-poiss[sample(1:nrow(poiss),size=10000,replace = F,prob = exp(poiss[,1])),]



Variable<-Value<-Model<-c()
for(i in 1:ncol(poiss)) {
  Variable%<>%c(rep(colnames(poiss)[i],nrow(poiss)))
  Value%<>%c(poiss[,i])
}
Model%<>%c(rep("Poisson Obs Model",length(poiss)))

for(i in 1:ncol(multnom)) {
  Variable%<>%c(rep(colnames(multnom)[i],nrow(multnom)))
  Value%<>%c(multnom[,i])
}
Model%<>%c(rep("Multinomial Obs Model",length(multnom)))

shdens<-data.frame(Variable=Variable, Value=Value, Model=Model)

lSHEEP<-GetSoaySheep(directory,oneSex=T)
GLM<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb)))
GLMCI<-do.call(getInitialValues_R,c(lSHEEP[c("solveDF","detectedNum")],list(fixedObsProb=fixedObsProb,CI = T)))
loGLM<-c(GLMCI$survPars[,1],GLMCI$growthPars[1:2,1],log(GLMCI$growthPars[3,1]),GLMCI$reprPars[,1],GLMCI$offNumPars[1],GLMCI$offSizePars[1:2,1],log(GLMCI$offSizePars[3,1]),5.2)
hiGLM<-c(GLMCI$survPars[,2],GLMCI$growthPars[1:2,2],log(GLMCI$growthPars[3,2]),GLMCI$reprPars[,2],GLMCI$offNumPars[2],GLMCI$offSizePars[1:2,2],log(GLMCI$offSizePars[3,2]),6.1+(6.1-5.2))
abliny<-data.frame(Variable=colnames(old15)[-1],Z=unname(unlist(GLM)),Y=loGLM,X=hiGLM)

q<-shdens%>%filter(Variable!="Log-Target")%>%ggplot(aes(Value,group=Model))+geom_density(aes(colour=Model,fill=Model),alpha=0.5)+
geom_vline(data = abliny, aes(xintercept = Z))+geom_vline(data = abliny, aes(xintercept = X),linetype = "dashed")+geom_vline(data = abliny, aes(xintercept = Y),linetype = "dashed")
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q


Sheep1<-readRDS("Results/REAL_GSF_beta_poissonMu_poissonObs_GLMx0_30000_15brks_autoshift_uniform_its30000_2021-12-20_220239")
Sheep2<-readRDS("Results/REAL_GSF_beta_poissonMu_poissonObs_GLMx0_30000_15brks_autoshift_uniform_its30000_2021-12-20_220315")
Sheep3<-readRDS("Results/REAL_GSF_beta_poissonMu_poissonObs_GLMx0_30000_15brks_autoshift_uniform_its30000_2021-12-20_220340")
Sheep4<-readRDS("Results/REAL_GSF_beta_poissonMu_poissonObs_GLMx0_30000_15brks_autoshift_uniform_its30000_2021-12-20_220719")
Sheep5<-readRDS("Results/REAL_GSF_beta_poissonMu_poissonObs_GLMx0_30000_15brks_autoshift_uniform_its30000_2021-12-20_223013")

SIMpois<-c()
for (i in 1:5) SIMpois%<>%rbind(get(paste0("Sheep",i)))
# SIMpois<-SIMpois[sample(1:nrow(SIMpois),30000,replace = F,prob = exp(SIMpois[,1])),]

Sheep1<-readRDS("Results/REAL_GSF_beta_poissonMu_multinomialObs_GLMx0_15000_15brks_autoshift_uniform_its15000_2021-12-21_011039")
Sheep2<-readRDS("Results/REAL_GSF_beta_poissonMu_multinomialObs_GLMx0_15000_15brks_autoshift_uniform_its15000_2021-12-21_012958")
Sheep2<-Sheep2[1000:15000,]
Sheep1<-Sheep1[7000:15000,]

SIMmult<-c()
for (i in 1:2) SIMmult%<>%rbind(get(paste0("Sheep",i)))
# SIMmult<-SIMmult[sample(1:nrow(SIMmult),10000,replace = F,prob = exp(SIMmult[,1])),]

colnames(SIMpois)<-colnames(SIMmult)<- c("Log-Target",plotNamesEx)
Variable<-Value<-Model<-c()
for(i in 1:ncol(SIMmult)) {
  Variable%<>%c(rep(colnames(SIMmult)[i],nrow(SIMmult)))
  Value%<>%c(SIMmult[,i])
}
Model%<>%c(rep("Multinomial Obs Model",length(SIMpois)))

for(i in 1:ncol(SIMpois)) {
  Variable%<>%c(rep(colnames(SIMpois)[i],nrow(SIMpois)))
  Value%<>%c(SIMpois[,i])
}
Model%<>%c(rep("Poisson Obs Model",length(SIMmult)))

shdens<-data.frame(Variable=Variable, Value=Value, Model=Model)

vals<-c(-8.25, 3.77,
        1.41, 0.56, log(0.08),
        -7.23, 2.6,
        log(0.06),
        0.36, 0.71, log(0.16),
        qlogis(0.9),
        log(50),
        log(10))


abliny<-data.frame(Variable=colnames(SIMpois)[-1],Z=unname(unlist(vals)))

q<-shdens%>%filter(Variable!="Log-Target")%>%ggplot(aes(Value,group=Model))+geom_density(aes(colour=Model,fill=Model),alpha=0.5)+
  geom_vline(data = abliny, aes(xintercept = Z))
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q


Sheep1<-readRDS("Results/SIM_GSF_beta_poissonMu_multinomialObs_GLMx0_20000_15brks_autoshift_uniform_its20000_2021-12-21_175043")
Sheep2<-readRDS("Results/SIM_GSF_beta_poissonMu_multinomialObs_GLMx0_20000_15brks_autoshift_uniform_its20000_2021-12-21_180115")
Sheep3<-readRDS("Results/SIM_GSF_beta_poissonMu_multinomialObs_GLMx0_20000_15brks_autoshift_uniform_its20000_2021-12-21_180656")
Sheep4<-readRDS("Results/SIM_GSF_beta_poissonMu_multinomialObs_GLMx0_20000_15brks_autoshift_uniform_its20000_2021-12-21_185041")

SIMmult15<-c()
for (i in 1:4) SIMmult15%<>%rbind(get(paste0("Sheep",i))[500:20000,])
colnames(SIMmult15)<-colnames(SIMmult15)<- c("Log-Target",plotNamesEx)

Sheep1<-readRDS("Results/SIM_pop100_yr10_GSF_beta_poissonMu_multinomialObs_GLMx0_20000_20brks_autoshift_uniform_its20000_2021-12-21_192344")
Sheep2<-readRDS("Results/SIM_pop100_yr10_GSF_beta_poissonMu_multinomialObs_GLMx0_20000_20brks_autoshift_uniform_its20000_2021-12-21_194028")
SIMmult20<-c()
for (i in 1:2) SIMmult20%<>%rbind(get(paste0("Sheep",i))[500:20000,])
colnames(SIMmult20)<-colnames(SIMmult20)<- c("Log-Target",plotNamesEx)


Variable<-Value<-Model<-c()
for(i in 1:ncol(SIMmult15)) {
  Variable%<>%c(rep(colnames(SIMmult15)[i],nrow(SIMmult15)))
  Value%<>%c(SIMmult15[,i])
}
Model%<>%c(rep("Multinomial 15",length(SIMmult15)))

for(i in 1:ncol(SIMmult20)) {
  Variable%<>%c(rep(colnames(SIMmult20)[i],nrow(SIMmult20)))
  Value%<>%c(SIMmult20[,i])
}
Model%<>%c(rep("Multinomial 20",length(SIMmult20)))

shdens<-data.frame(Variable=Variable, Value=Value, Model=Model)

abliny<-data.frame(Variable=colnames(SIMmult20)[-1],Z=unname(unlist(vals)))

q<-shdens%>%filter(Variable!="Log-Target")%>%ggplot(aes(Value,group=Model))+geom_density(aes(colour=Model,fill=Model),alpha=0.5)+
  geom_vline(data = abliny, aes(xintercept = Z))
q<-q+facet_wrap(. ~ Variable,scales = "free") + theme(strip.text.x = element_text(size = 12))+
  xlab("MCMC Sample Value")+ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5)) ;q

FormPostSheep<-function(namer,chains=F,fixed=F){
  
  filerz<-list.files("./Results/",pattern = namer)
  if(fixed) plotNam<-plotNames[-13] else plotNam<-plotNamesEx
  
  output<-c()
  if(!chains){
    for (i in 1:length(filerz)) output%<>%rbind(readRDS(paste0("./Results/",filerz[i])))
    colnames(output)<-colnames(output)<- c("Log-Target",plotNam)
  } else {
    for (i in 1:length(filerz)) {
      tmp<-readRDS(paste0("./Results/",filerz[i]))
      colnames(tmp)<- c("Log-Target",plotNam)
      tmp%<>%mcmc()
      output%<>%c(list(tmp))
    }
    
    return(output)  
  }
}

vals<-c(-8.25, 3.77,
        1.41, 0.56, log(0.08),
        -7.23, 2.6,
        log(0.06),
        0.36, 0.71, log(0.16),
        qlogis(0.9),
        log(50),
        log(10)
)

# namer<-"SIM_pop100_yr10_GSF_beta_poissonMu_multinomialObs_GLMx0_10000_15brks_autoshift_uniform_its10000_2021-12"
namer<-"SIM_pop300_yr10_GSF_beta_poissonMu_multinomialObs_GLMx0_10000_10brks_autoshift_uniform_its10000_2021-12-"
out<-FormPostSheep(namer)
plot(out[,1])

vals
unname(apply(out,2,weighted.mean,w=out[,1]))[2:15]
unname(colMeans(out))[2:15]

namer<-"REAL_GSF_fixed_multMu_MultObs_GLMx0_final_10brks_autoshift_uniform_its20000_2021-12-1"

ExploreSheep<-function(namer){
  out<-FormPostSheep(namer,chains = T)
  outy<-grDiagnostic(out)
  # Number of effective samples
  Neff<-effectiveSize(out)
  # Rhat - Gelman-Rubin convergence parameter
  Rhat<-gelman.diag(out)$psrf[,1]
  # p-value & MSE determination
  
}






















yr<-c("5","10","20")
pop<-c("50","100","300","500")
fibe<-c("fixed","beta")
mu<-c("multinomialMu","poissonMu")
obs<-c("binomialObs","poissonObs","multinomialObs")
brks<-c("5","10","15","20","30")
samp<-c("sampleDTN","sampleNorm")
shift<-c("autoshift","manshift")
its<-c("10000","20000")

# ppp<-"100"
# yyy<-"10"
# fff<-"beta"
# mmm<-"poissonMu"
# ooo<-"multinomialObs"
# bbb<-"10"
# sss<-"sampleDTN"
# shsh<-"autoshift"
# iii<-"10000"

vals<-c(-8.25, 3.77,
        1.41, 0.56, log(0.08),
        -7.23, 2.6,
        log(0.06),
        0.36, 0.71, log(0.16),
        qlogis(0.9),
        log(50),
        log(10)
)

vals<-c(-7.25, 3.77,
            1.41, 0.56, log(0.08),
            -7.23, 2.6,
            log(0.06),
            0.36, 0.71, log(0.16),
            qlogis(0.9),
            log(50),
            log(10)
)

vals%<>%unlist()

allresults<-list()

for (ppp in pop){
  
  for (yyy in yr){
    
    for (fff in fibe){
      
      for (mmm in mu){
        
        for (ooo in obs){
          
          for (bbb in brks){
            
            for (sss in samp){
              
              for (shsh in shift){
                
                for (iii in its){
                  
                  pathy<-paste0("SIM_pop",ppp,"_yr",yyy,"_GSF_",fff,"_",mmm,"_",ooo,"_GLMx0_",iii,"_",bbb,"brks_",sss,"_",shsh)
                  filerz<-list.files(pattern = pathy,path = "./Results")
                  if(is.null(filerz)) sammy<-"sampleDTN" else sammy<-sss
                  pathy<-paste0("SIM_pop",ppp,"_yr",yyy,"_GSF_",fff,"_",mmm,"_",ooo,"_GLMx0_",iii,"_",bbb,"brks_",shsh)
                  filerz%<>%c(list.files(pattern = pathy,path = "./Results"))
                  filerz<-filerz[!grepl("SimulatedData_",filerz)]
                  
                  if(length(filerz)==0) next
                  print(" ")
                  print(pathy)
                  
                  lennie<-length(filerz)
                  
                  Totale<-sapply(1:lennie,function(i) readRDS(paste0("./Results/",filerz[i])),simplify = F)
                  
                  cBtot<-list()
                  for (i in 2:ncol(Totale[[1]])) cBtot%<>%c(list(sapply(1:lennie,function(j) Totale[[j]][,i])))
                  rBtot<-matrix(unlist(Totale),nrow = lennie*nrow(Totale[[1]]),ncol = ncol(Totale[[1]]))
                  
                  rBtot<-Totale[[1]]
                  if (lennie>1){
                    for (i in 2:lennie) rBtot%<>%rbind(Totale[[i]])
                  }
                  
                  output<-data.frame(True=vals[1:ifelse(fff=="fixed",length(vals)-2,length(vals))])
                  # output<-data.frame(predMean=colMeans(rBtot[,-1]))
                  # output$True<-vals[1:12]
                  
                  # Mean & ConfInt
                  output$predMean<-colMeans(rBtot[,-1])
                  output$predMAP<-rBtot[which.max(rBtot[,1]),-1]
                  output$predMedian<-apply(rBtot[,-1],2,median)
                  output$predLow<-apply(rBtot[,-1],2,quantile,probs=0.05)
                  output$predHigh<-apply(rBtot[,-1],2,quantile,probs=0.95)
                  
                  output$pvalue<-sapply(2:ncol(rBtot),function(i) 2*min(sum(rBtot[,i]<=vals[i],na.rm = T),
                                                                        sum(rBtot[,i]>=vals[i],na.rm = T))/sum(!is.na(rBtot[,i])))
                  
                  output$RMSEmean<-sqrt(colMeans((rBtot[,-1]-output$True)^2))
                  output$RMSEmeanMLE<-sqrt((rBtot[1,-1]-output$True)^2) # first value was generated by MLE on the simulated data
                  
                  # Need multiple chains per simulation for Rhat
                  output$Rhat<-sapply(1:length(cBtot),function(i) rstan::Rhat(cBtot[[i]]))
                  output$essBulk<-sapply(1:length(cBtot),function(i) rstan::ess_bulk(cBtot[[i]]))
                  output$essTail<-sapply(1:length(cBtot),function(i) rstan::ess_tail(cBtot[[i]]))

                  metadata<-list(Nchains=lennie,pop=ppp,yrs=yyy,fibe=fff,mu=mmm,obs=ooo,brks=bbb,samp=sammy,shift=shsh)
                  
                  allresults%<>%c(list(results=output,metadata=metadata))
                  
                  print(dplyr::select(output,True,predMean,pvalue,RMSEmean,Rhat,essBulk,essTail))
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
}



Sheep1<-readRDS("Results/SIM_pop100_yr10_GSF_fixed_poissonMu_multinomialObs_GLMx0_30000_5brks_regbinspaceFALSE_sampleDTN_autoshift_uniform_its30000_2022-07-15_030311")
Sheep2<-readRDS("Results/SIM_pop100_yr10_GSF_fixed_poissonMu_multinomialObs_GLMx0_30000_10brks_regbinspaceFALSE_sampleDTN_autoshift_uniform_its30000_2022-07-15_040235")
Sheep3<-readRDS("Results/SIM_pop100_yr10_GSF_fixed_poissonMu_multinomialObs_GLMx0_30000_15brks_regbinspaceFALSE_sampleDTN_autoshift_uniform_its30000_2022-07-15_043946")


