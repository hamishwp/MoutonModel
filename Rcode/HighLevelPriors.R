
# Do we make the HLP based off the GLM, such as by taking the values +/-5 sds?
# 

# Need to work on the following:
#   Surv & repr params (all 4)
#   Offspring size params


# Put in positive link functions for the following:
# 1) gradient of growth
# 2) gradient offspring size
# 3) gradient of survival
# 4) gradient of reproduction

# Find an upper bound on any of the above gradients?

######### make values by assuming NA values are all positive or all negative for survival and reproduction
lSHEEP$solveDF%>%mutate(s_bins = cut(size, breaks = IPMLTP$breaks))%>%
  group_by(s_bins)%>%summarise(repr=mean(reproduced,na.rm = T))

lSHEEP$solveDF%>%mutate(s_bins = cut(prev.size, breaks = IPMLTP$breaks))%>%
  group_by(s_bins)%>%summarise(survy=mean(survived,na.rm = T))


### SCHILD ###
# Offspring survival probability
SchildEst<-boot::boot(lSHEEP$solveDF$survived[lSHEEP$solveDF$prev.size<IPMLTP$sizes[2]],
                      function(u,i) mean(u[i],na.rm=T),R=1000)
SchildEst<-boot::boot.ci(SchildEst,type=c("norm","basic","perc"),conf = 0.95)
minSchild<-SchildEst$percent[4]
maxSchild<-SchildEst$percent[5]
schildHLP<-function(schild) ifelse(schild>maxSchild | schild<minSchild, F,T)

### SURVIVAL ###
# Survival Gradient
minSurv2<-0
maxSurv2<-Inf
# Survival Intercept, given that min gradient >=0
survy<-lSHEEP$solveDF%>%group_by(census.number)%>%
  summarise(survy=sum(survived,na.rm = T)/length(survived[!is.na(survived)]))%>%
  pull(survy)
if(mean(survy)>(1-1e-5) | median(survy)>(1-1e-5)) {
  maxSurv1 <- Inf
} else {
  sursize <- min(c(max(survy[survy!=1]),(1-1e-5)))
  maxSurv1 <- -log(1/sursize-1)
}
minSurv1<--Inf
# Higher level prior as a combination of the two parameters 
# (note this is based on survival of lowest size class - Schild - as we already have upper-surv1 and lower-surv2)
survHLP<-function(survPars) {
  ifelse(survPars[1]<minSurv1 | survPars[2]<minSurv2 | survPars[1]>maxSurv1 | survPars[2]>maxSurv2,F,T)
  # ifelse(1/(1+exp(-IPMLTP$sizes[1]*survPars[2]-survPars[1])) > SchildEst$percent[5],F,T)
}

### REPRODUCTION ###
# Reproduction Gradient
minRepr2<-0
maxRepr2<-Inf
# Reproduction Intercept, given that min gradient >=0
repry<-lSHEEP$solveDF%>%group_by(census.number)%>%
  summarise(repr=mean(reproduced,na.rm = T))%>%pull(repr)
if(mean(repry)>(1-1e-5) | median(repry)>(1-1e-5)) {
  maxRepr1 <- Inf
} else {
  repsize <- min(c(max(repry[repry!=1]),(1-1e-5)))
  maxRepr1 <- -log(1/repsize-1)
}
minRepr1<--Inf
# Higher level prior as a combination of the two parameters 
# (note this is based on reproduction of lowest size class - Schild - as we already have upper-surv1 and lower-surv2)
RchildEst<-boot::boot(lSHEEP$solveDF$reproduced[lSHEEP$solveDF$size<IPMLTP$sizes[2]],
                      function(u,i) mean(u[i],na.rm=T),R=1000)
RchildEst<-boot::boot.ci(RchildEst,type=c("norm","basic","perc"),conf = 0.95)
# Create the function
reprHLP<-function(reprPars) {
  ifelse(reprPars[1]<minRepr1 | reprPars[2]<minRepr2 | reprPars[1]>maxRepr1 | reprPars[2]>maxRepr2,F,T)
  # ifelse(1/(1+exp(-IPMLTP$sizes[1]*reprPars[2]-reprPars[1])) > RchildEst$percent[5],F,T)
}






# sizes<-IPMLTP$sizes
# minsize<-c(IPMLTP$DTN$L,IPMLTP$sizes[1:length(IPMLTP$sizes)-1])
# maxsize<-c(IPMLTP$sizes[2:length(IPMLTP$sizes)],IPMLTP$DTN$U)
# lSHEEP$solveDF%>%mutate(s_bins = cut(size, breaks = IPMLTP$breaks))%>%group_by(s_bins)%>%
#   summarise(mingrow=min(size-prev.size,na.rm = T),maxgrow=max(size-prev.size,na.rm = T))
# growers<-lSHEEP$solveDF%>%mutate(s_bins = cut(size, breaks = IPMLTP$breaks))%>%group_by(s_bins)%>%
#   summarise(mingrow=quantile(size-prev.size,0.05,na.rm = T),maxgrow=quantile(size-prev.size,0.95,na.rm = T))
# growers<-growers[-nrow(growers),]
# growers$realmin<-maxsize*vals$growthPars[2]+vals$growthPars[1]-5*vals$growthPars[3]-maxsize
# growers$realmax<-maxsize*vals$growthPars[2]+vals$growthPars[1]+5*vals$growthPars[3]-maxsize

# Growth: 
minGrow1<--Inf
minGrow2<-0
minGrow3<-0
maxGrow1<-Inf
maxGrow2<-Inf
maxGrow3<-Inf

growHLP<-function(growPars) {
  ifelse(growPars[1]<minGrow1 | growPars[2]<minGrow2 |
         growPars[3]<minGrow3 | growPars[3]>maxGrow3 |
         growPars[1]>maxGrow1 | growPars[2]>maxGrow2,F,T)
}

# Offspring size
minOsize1<--Inf
minOsize2<-0
minOsize3<-0
maxOsize1<-Inf
maxOsize2<-Inf
maxOsize3<-Inf

osizHLP<-function(offSizePars) {
  ifelse(offSizePars[1]<minOsize1 | offSizePars[2]<minOsize2 |
           offSizePars[3]<minOsize3 | offSizePars[3]>maxOsize3 |
           offSizePars[1]>maxOsize1 | offSizePars[2]>maxOsize2,F,T)
}

# Offspring number
OnumEst<-boot::boot(lSHEEP$solveDF$off.born,
                      function(u,i) mean(u[i],na.rm=T),R=1000)
OnumEst<-boot::boot.ci(OnumEst,type=c("norm","basic","perc"),conf = 0.95)

onumHLP<-function(offNumPars) ifelse(offNumPars<OnumEst$percent[4] | 
                                     offNumPars>OnumEst$percent[5], F, T)

# 1) Calculate the maximum and minimum growth/loss
# 2) Make sure that the expected value +/- 1 sigma is less/more than max/min resp.
# 3) Do step 2 for all size bins

# Offspring size:
# 1) Calculate the maximum and minimum parent-size/offspring-size
# 2) Make sure that the expected value +/- 1 sigma is less/more than max/min resp.
# 3) Do step 2 for all size bins which have reproducing parents

if(HLPon){
  IPMLTP$HLP<-function(theta,lTargPars){
    apply(theta,1, function(x0){
    vals<-Sample2Physical(x0,lTargPars)
    ifelse(!all(c(survHLP(vals$survPars),
                  growHLP(vals$growthPars),
                  reprHLP(vals$reprPars),
                  onumHLP(vals$offNumPars),
                  osizHLP(vals$offSizePars),
                  schildHLP(vals$Schild))),F,T)
  })}
} else IPMLTP$HLP<-function(theta,lTargPars) rep(T,nrow(theta))




