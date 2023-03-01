################################################################################
################### OBSERVATION CONDITIONAL LIKELIHOODS ########################
################################################################################

binomialObs <- function(Y, X, p){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability, so that the number of observed
  #           animals in each size class is binomially distributed with number
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # For a single time step:
  helper <- function(yx, p){
    n <- length(yx)
    y <- yx[1:n]
    x <- yx[-(1:n)]
    dbinom(y, x, p) %>% log %>% sum %>% `-` %>% return
  }
  
  # Loop through all time steps and return:
  rbind(Y, X) %>% apply(2, helper, p=p) %>% return
}

detectionNumObs1D <- function(Y, X, p, logy=F){
  # purpose : A simple model which ignores the probability of the distribution
  #           of sizes, and says that the probability of observing what we
  #           observe is simply the probability of detecting the number of 
  #           animals we did, given the true distribution.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  trueSizes <- sum(X)
  obsSizes <- sum(Y)
  
  return(dbinom(obsSizes, trueSizes, p, log=logy))
}

detectionNumObs <- function(Y, X, p, logy=T){
  # purpose : A simple model which ignores the probability of the distribution
  #           of sizes, and says that the probability of observing what we
  #           observe is simply the probability of detecting the number of 
  #           animals we did, given the true distribution.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  return(dbinom(colSums(Y), colSums(X), p, log=logy))
  
}


poissonObs <- function(Y, X, p, logy=T){
  # purpose : A model which assumes the number of individuals in each size class
  #           is poisson distributed with mean equal to the number of true
  #           individuals in the classx the detection probability.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  results <- dpois(Y+1, (X*p)+1, log=logy)
  # results %<>% replace(is.infinite(results), log(.Machine$double.xmin))
  if (logy) return(colSums(results,na.rm = T))
  else return(colprods(results,na.rm = T))
}

prob1Obs <- function(Y, X, p, logy=T){
  # purpose : returns a uniform probability so that we can test the particle 
  #           filter more easily
  # inputs  : Y   - The observation matrix (not used)
  #           X   - The True state matrix (not used)
  #           p   - The probability of detection (not used)
  #           log - if TRUE, returns the log probability instead
  Y %<>% as.matrix %>% ncol
  if (logy) return(runif(Y) %>% log)
  else return(runif(Y))
}

SplitMNom<-function(vec,logy) dmultinom(x = vec[,1],prob=vec[,2],log=logy)

dmnom<-function(Y,multiProbs,logy) apply(abind(Y,multiProbs,along=3),2,SplitMNom,log=logy)

multinomialObs <- function(Y, X, p, logy=T){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability. Then, given the total number of 
  #           observed animals, we assumed that they are multinomially
  #           distributed according to the same size distribution as the true
  #           population
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # Sometimes the particle filter requires all arguments to be a list, since
  # p is only a single parameter, this can cause it to be passed as a list of 
  # one element, this tries to elegantly cope with the edge case:
  p %<>% as.numeric
  
  # Some required quantities:
  truePopSizes <- colSums(X)
  obsPopSizes <- colSums(Y)
  multiProbs <- apply(X+1, 2, countsToProbs)

  out1<-dmnom(Y,multiProbs,logy)
  out2<-dbinom(obsPopSizes+1, truePopSizes+1, p, log = logy)
  
  if (logy) {
    outout<-out1+out2
    outout[outout<log(.Machine$double.xmin)]<-log(.Machine$double.xmin)
    return(outout)
  }else {
    return(out1*out2)
  }
  
}

multinomPoisObs <- function(Y, X, p, logy=T){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability. Then, given the total number of 
  #           observed animals, we assumed that they are multinomially
  #           distributed according to the same size distribution as the true
  #           population
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # Sometimes the particle filter requires all arguments to be a list, since
  # p is only a single parameter, this can cause it to be passed as a list of 
  # one element, this tries to elegantly cope with the edge case:
  p %<>% as.numeric
  
  # Some required quantities:
  truePopSizes <- colSums(X)
  obsPopSizes <- colSums(Y)
  multiProbs <- apply(X+1, 2, countsToProbs)
  # multiProbs[is.na(multiProbs)]<-1
  
  # out1<- dmultinomial(t(Y), prob = t(multiProbs), log = logy)
  out1<-dmnom(Y,multiProbs,logy)
  out2<-dpois(obsPopSizes+1, (truePopSizes*p)+1, log=logy)
  
  if (logy) {
    return(out1+out2)
  }else {
    return(out1*out2)
  }
}

beta_binObs <- function(Y, X, p, logy=T, time, NoParts){
  # purpose : A model which assumes the number of individuals in each size class
  #           is poisson distributed with mean equal to the number of true
  #           individuals in the classx the detection probability.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  pobs<-rbeta(ncol(Y),p[1],p[2])
  
  if (logy) return(dbinom(colSums(X), colSums(Y), pobs, log=logy))
  else return(dbinom(colSums(X), colSums(Y), pobs, log=logy))
}

beta_poisObs <- function(Y, X, p, logy=T, time, NoParts){
  # purpose : A model which assumes the number of individuals in each size class
  #           is poisson distributed with mean equal to the number of true
  #           individuals in the classx the detection probability.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  pobs<-rbeta(ncol(Y),p[1],p[2])
  results <- dpois(Y, (X*pobs), log=logy)
  # results %<>% replace(is.infinite(results), log(.Machine$double.xmin))
  if (logy) return(colSums(results))
  else return(colprods(results))
}

beta_mnomObs<-function(Y, X, p, logy=T, time, NoParts){
  # purpose : An observation model for size distributions. We assume we miss an 
  #           animal with fixed probability. Then, given the total number of 
  #           observed animals, we assumed that they are multinomially
  #           distributed according to the same size distribution as the true
  #           population, using a beta distribution to generate the observed probabilities
  #           of trials equal to the true number of animals in that size class.
  # inputs  : Y - A vector or matrix containing the observed sizes. If it is a
  #               matrix, it is assumed that each time step is a column.
  #           X - A vector or matrix containing the true sizes. The columns are
  #               also assumed to be for different time steps if this is a
  #               matrix.
  #           p - The probability of observing an animal given it is present in
  #               the true population.
  # output  : Either a single probability, or a vector with the probability for
  #           each time step. Log probabilities are used.
  
  # Sometimes the particle filter requires all arguments to be a list, since
  # p is only a single parameter, this can cause it to be passed as a list of 
  # one element, this tries to elegantly cope with the edge case:
  
  # Some required quantities:
  truePopSizes <- colSums(X)
  obsPopSizes <- colSums(Y)
  multiProbs <- apply(X, 2, countsToProbs)
  indy<-truePopSizes!=0
  
  pobs<-rbeta(length(truePopSizes),p[1],p[2])
  
  # out1<- dmultinomial(t(Y), prob = t(multiProbs), log = logy)
  # Make sure that if any X distributions have zero individuals, the likelihood equals zero
  out1<-out2<-rep(ifelse(logy,-Inf,0),length(truePopSizes))
  out1[indy]<-dmnom(Y[,indy],multiProbs[,indy],logy)
  out2[indy]<-dbinom(obsPopSizes[indy], truePopSizes[indy], pobs[indy], log = logy)
  
  if (logy) return(out1+out2)
  else return(out1*out2)
  
}

fixedMuObs<-function(Y, X, p, logy=T, time, NoParts){
  return(multinomialObs(Y[,rep(time,NoParts)], X, p[time], logy=logy))
}

fixedPoisObs<-function(Y, X, p, logy=T, time, NoParts){
  return(poissonObs(Y[,rep(time,NoParts)], X, p[time], logy=logy))
}

fixedBinObs<-function(Y, X, p, logy=T, time, NoParts){
  return(detectionNumObs(Y[,rep(time,NoParts)], X, p[time], logy=logy))
}

# Minkowski distance function
# Distances were taken from https://www.biorxiv.org/content/10.1101/2021.07.29.454327v1.full.pdf
# Note that if p=1 we are using the L1 distances - our default
MADadaptdist<-function(sest,sobs,p=1){
  # In case nothing works at all
  if(all(sest==0)) return(list(shat=rep(0,length(sobs)), sw=rep(-1e300,ncol(sest)))) 
  # Median of columns
  meds<-apply(sest,1,median,na.rm=T)
  # Median Absolute Deviation (MAD) to the sample
  MAD<-abs(sest-meds)
  # Median Absolute Deviation to the Observation (MADO)
  MADO<-abs(sest-sobs)
  # Don't punish the summary statistics that deviate more from obs data initially than others:
  if(sum(MADO>2*MAD,na.rm = T)/length(sobs)<1/3) PCMAD<-MAD+MADO else PCMAD<-MAD
  # Calculate the distance per summary statistic
  d_i<--MADO/(PCMAD+1e-5); d_i[,apply(PCMAD,2,sd)==0]<--1e300
  # Particle weight calculation
  sw<--abs(apply(d_i,2,function(dd) pracma::nthroot(median(dd^p,na.rm = T),p)))
  # output total distance
  return(list(shat=meds, sw=sw))
}

MAEdistVar<-function(sest,sobs,p=1){
  # Now calculate the variance for the distance normalisation
  sds<-apply(sest,2,sd,na.rm=T)
  # Median of columns
  meds<-apply(sest,1,median,na.rm=T)
  # Calculate the distance per summary statistic
  d_i<--abs(sest-sobs)/(sds+1e-5); d_i[,sds==0]<--1e300
  # Particle weight calculation
  sw<--abs(apply(d_i,2,function(dd) pracma::nthroot(mean(dd^p,na.rm = T),p)))
  # output total distance
  return(list(shat=meds, sw=sw)) 
}

MAEdist<-function(sest,sobs,p=1){
  # Median of columns
  meds<-apply(sest,1,median,na.rm=T)
  # Calculate the distance per summary statistic
  d_i<--abs(sest-sobs)
  # Particle weight calculation
  sw<--abs(apply(d_i,2,function(dd) pracma::nthroot(median(dd^p,na.rm = T),p)))
  # output total distance
  return(list(shat=meds, sw=sw)) 
}

################################################################################
######################### DEFINE DISTANCE METRIC ###############################
################################################################################

# Combining the distance function with the observation model - 
if(obsModel=="MultinomObs"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function - multinomial
    MultiMod<-function(i) {
      multiProbs <- countsToProbs(wArgs$Sd[,i,wArgs$time]+1)
      
      dbinom(sum(wArgs$Sd[,i,wArgs$time])+1, colSums(wArgs$Sstar[,i,])+1, pobs, log = T) +
        apply(wArgs$Sstar[,i,]+1,2,function(xx) dmultinom(x = xx,prob=multiProbs,log=T))
    }
    # For each type of summary statistic
    sw<-rowSums(sapply(1:3,MultiMod))
    # Calc median shat value
    shat<-apply(wArgs$Sstar,1:2,median,na.rm=T)*pobs
    
    return(list(sw=maxESS(sw,mag=30),shat=shat))
  }
  
} else if(obsModel=="multinomPoisObs"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function - multinomial-Poisson
    MultiMod<-function(i) {
      multiProbs <- countsToProbs(wArgs$Sd[,i,wArgs$time]+1)
      
      dpois(sum(wArgs$Sd[,i,wArgs$time])+1, (colSums(wArgs$Sstar[,i,])*pobs)+1, log=T) +
      apply(wArgs$Sstar[,i,]+1,2,function(xx) dmultinom(x = xx,prob=multiProbs,log=T))
    }
    # For each type of summary statistic
    sw<-rowSums(sapply(1:3,MultiMod))
    # Calc median shat value
    shat<-apply(wArgs$Sstar,1:2,median,na.rm=T)*pobs
    
    return(list(sw=maxESS(sw,mag=30),shat=shat))
  }
} else if(obsModel=="multinomMAE"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function - multinomial-Poisson
    MultiMod<-function(i) {
      multiProbs <- countsToProbs(wArgs$Sd[,i,wArgs$time]+1)
      
      -abs(sum(wArgs$Sd[,i,wArgs$time])-colSums(wArgs$Sstar[,i,])*pobs) +
        apply(wArgs$Sstar[,i,]+1,2,function(xx) dmultinom(x = xx,prob=multiProbs,log=T))
    }
    # For each type of summary statistic
    sw<-rowSums(sapply(1:3,MultiMod))
    # Calc median shat value
    shat<-apply(wArgs$Sstar,1:2,median,na.rm=T)*pobs
    
    return(list(sw=maxESS(sw,mag=30),shat=shat))
  }
}else if (obsModel=="PoisObs"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function - poisson
    MultiMod<-function(i) sum(poissonObs(wArgs$Sd[,,wArgs$time], 
                                         wArgs$Sstar[,,i], 
                                         pobs, 
                                         logy=T))
    # For each type of summary statistic
    sw<-sapply(1:dim(wArgs$Sstar)[3],MultiMod)
    # Calc median shat value
    shat<-apply(wArgs$Sstar,1:2,median,na.rm=T)*pobs
    
    return(list(sw=maxESS(sw,mag=30),shat=shat))
  }
  
} else if (obsModel=="BinomObs"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function - binomial
    MultiMod<-function(i) sum(detectionNumObs(wArgs$Sd[,,wArgs$time], 
                                              wArgs$Sstar[,,i],
                                              pobs, 
                                              logy=T))
    # For each type of summary statistic
    sw<-sapply(1:dim(wArgs$Sstar)[3],MultiMod)
    # Calc median shat value
    shat<-apply(wArgs$Sstar,1:2,median,na.rm=T)*pobs
    
    return(list(sw=maxESS(sw,mag=30),shat=shat))
  }
  
} else if (obsModel=="MADadaptdist"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function
    Sstar<-apply(wArgs$Sstar*pobs,3,rbind)
    # Calculate the distances - 
    disties<-MADadaptdist(sest = Sstar,
                      sobs = c(wArgs$Sd[,,wArgs$time]))
    # Make sure that anytime the summary stats are zero across all bins are zero
    mostlyempty<-all(apply(apply(wArgs$Sstar*pobs,2:3,function(x) sum(x)==0),2,function(x) any(x)))
    if(mostlyempty) disties$sw[]<--1e300

    return(list(sw=disties$sw,
                shat=array(disties$shat,dim(wArgs$Sstar)[1:2])*pobs))
  }
  
} else if (obsModel=="MAEdistVar"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function
    Sstar<-apply(wArgs$Sstar*pobs,3,rbind)
    # Calculate the distances - 
    disties<-MAEdistVar(sest = Sstar,
                     sobs = c(wArgs$Sd[,,wArgs$time]))
    # Make sure that anytime the summary stats are zero across all bins are zero
    mostlyempty<-all(apply(apply(wArgs$Sstar*pobs,2:3,function(x) sum(x)==0),2,function(x) any(x)))
    if(mostlyempty) disties$sw[]<--1e300
    
    return(list(sw=disties$sw,
                shat=array(disties$shat,dim(wArgs$Sstar)[1:2])*pobs))
  }
  
} else if (obsModel=="MAEdist"){
  
  ObsDistance<-function(wArgs,pobs){
    # Distance function
    Sstar<-apply(wArgs$Sstar*pobs,3,rbind)
    # Calculate the distances - 
    disties<-MAEdist(sest = Sstar,
                     sobs = c(wArgs$Sd[,,wArgs$time]))
    # Make sure that anytime the summary stats are zero across all bins are zero
    mostlyempty<-all(apply(apply(wArgs$Sstar*pobs,2:3,function(x) sum(x)==0),2,function(x) any(x)))
    if(mostlyempty) disties$sw[]<--1e300
    
    return(list(sw=disties$sw,
                shat=array(disties$shat,dim(wArgs$Sstar)[1:2])*pobs))
  }
  
}

# To choose the distance function, please see 'ObsDistance.R'
if(fixedObsProb){
  obsfun<-function(output,wArgs){
    # Calculate the distances
    diz<-ObsDistance(wArgs,pobs=wArgs$pobs[wArgs$time])
    # Calculate the total distances
    output$d<-output$d+sum(diz$sw,na.rm = T)
    # Exponentiate and scale the particle weights to become from 0 to 1
    output$sw<-exp(diz$sw-max(diz$sw,na.rm = T))
    summz<-sum(output$sw,na.rm = T)
    # Normalise it!
    if(summz==0 | sum(!is.na(output$sw))<2) {
      output$sw<-rep(1/length(output$sw),length(output$sw))
    } else output$sw<-output$sw/summz
    # Add the median summary statistics data as well:
    output$shat[,,wArgs$time]<-diz$shat
    
    return(output)
  }
} else {
  obsfun<-function(output,wArgs){
    # Simulate the obsProbPar from the beta distribution here
    pobs<-rbeta(dim(wArgs$Sstar)[3],wArgs$pobs[1],wArgs$pobs[2])
    # Calculate the total distances
    output$d<-output$d+sum(diz$sw,na.rm = T)
    # Exponentiate and scale the particle weights to become from 0 to 1
    output$sw<-exp(diz$sw-max(diz$sw,na.rm = T))
    summz<-sum(output$sw,na.rm = T)
    # Normalise it!
    if(summz==0 | sum(!is.na(output$sw))<2) {
      output$sw<-rep(1/length(output$sw),length(output$sw))
    } else output$sw<-output$sw/summz
    # Add the median summary statistics data as well:
    output$shat[,,wArgs$time]<-diz$shat
    
    return(output)
  }
}




















