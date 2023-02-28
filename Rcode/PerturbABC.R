################# Perturbation functions for the ABCSMC algorithm?
# MVN based on resampled theta with a global covariance matrix
pert_GlobCov<-function(outin,lTargPars){
  # Make sure only successful samples contribute to the perturbation
  inds<-!is.na(outin$distance) & outin$distance>=outin$delta[outin$iteration]
  # Reduce ewights and previous parameter samples for computation
  wewei<-outin$weightings[inds]; theta<-outin$theta[inds,]
  # Calculate the global weighted mean and covariance matrix
  COVhat<-cov.wt(theta,wewei)
  # Then sample from the MVN distribution
  ResampleSIR<-function(NN) {
    # First sample particles based on the distance function
    thth<-theta[sample(1:sum(inds),NN,T,prob = wewei),]
    # Perturb the sample
    thth<-t(apply(thth,1,function(tt) Rfast::rmvnorm(1, tt, COVhat$cov)))
    # Transition kernel K(theta_i | theta_j) for j = 1,..,N of old particles
    # NOTE WE RELY ON SYMMETRY OF MVN(theta1, theta2, COV) = MVN(theta2, theta1, COV)
    ptrans<-apply(thth,1,function(tt) Rfast::dmvnorm(theta, tt, COVhat$cov)%*%wewei)
    # Priors
    prpr<-exp(apply(thth,1,function(tt) lTargPars$priorF(tt)))
    # Associated weights NOTE THAT NORMALISATION IS DONE LATER
    newW<-prpr/ptrans
    
    return(list(theta=thth,weightings=newW))
  }
  # Return the function as output
  return(ResampleSIR)
}

# MV skew-norm based on resampled theta with a global skew-covariance matrix
pert_GlobSkewCov<-function(outin,lTargPars){
  # Make sure only successful samples contribute to the perturbation
  inds<-!is.na(outin$distance) & outin$distance>=outin$delta[outin$iteration]
  # Reduce ewights and previous parameter samples for computation
  wewei<-outin$weightings[inds]; theta<-outin$theta[inds,]
  # Calculate the global mean, skew and covariance matrix
  disty<-sn::msn.mle(y = theta)$dp
  # Adjust such that normally distributed values remain un-skewed
  for (i in 1:ncol(theta)) if(!lTargPars$skew[i]) disty$alpha[i]<-0
  # Then sample from the MV-SN distribution
  ResampleSIR<-function(NN) {
    # First sample particles based on the distance function
    thth<-theta[sample(1:nrow(theta),NN,T,prob = wewei),]
    # Perturb the sample
    thth<-t(apply(thth,1,function(tt) sn::rmsn(1, xi = tt, Omega = disty$Omega, alpha = disty$alpha)))
    acc<-lTargPars$HLP(thth,lTargPars)
    while(sum(!acc)>0){
      thth[!acc,]<-theta[sample(1:nrow(theta),sum(!acc),T,prob = wewei),]
      if(sum(!acc)>1){
        thth[!acc,]<-t(apply(thth[!acc,],1,function(tt) sn::rmsn(1, xi = tt, Omega = disty$Omega, alpha = disty$alpha)))
        acc[!acc]<-lTargPars$HLP(thth[!acc,],lTargPars)
      } else {
        thth[!acc,]<-as.numeric(t(sn::rmsn(1, xi = thth[!acc,], Omega = disty$Omega, alpha = disty$alpha)))
        acc<-lTargPars$HLP(thth,lTargPars)
      }
    }
    # Transition kernel K(theta_i | theta_j) for j = 1,..,N of old particles
    # NOTE WE RELY ON SYMMETRY OF MVN(theta1, theta2, COV) = MVN(theta2, theta1, COV)
    ptrans<-apply(thth,1,function(tt) sn::dmsn(theta, xi = tt, Omega = disty$Omega, alpha=disty$alpha)%*%wewei)
    # Priors
    prpr<-exp(apply(thth,1,function(tt) lTargPars$priorF(tt)))
    # Associated weights NOTE THAT NORMALISATION IS DONE LATER
    newW<-prpr/ptrans
    newW[is.infinite(newW)]<-max(newW[!is.infinite(newW)],na.rm = T)
    newW[newW==0]<-min(newW[newW!=0])
    
    return(list(theta=thth,weightings=newW))
  }
  # Return the function as output
  return(ResampleSIR)
}





















# MVN based on resampled theta with a partly guided centering and global covariance matrix
pert_GuidedGlobCov<-function(output,SumStats,weightings=NULL){
  # Make sure only successful samples contribute to the perturbation
  inds<-output$distance<output$delta
  # Calculate the global weighted mean and covariance matrix
  COVhat<-cov.wt(t(output$shat[,inds]),weightings); lennie<-ncol(output$theta[inds,])
  # Calculate the guiding mean:
  guidMu<-weighted.mean(output$theta[inds],output$distance[inds])
  # and its weight:
  guidW<-mean(output$distance[inds],na.rm=T)
  # Then sample from the rnorm
  ResampleSIR<-function(NN) {
    # First sample particles based on the distance function
    iis<-sample((1:lennie)[inds],NN,T,prob = output$distance[inds])
    theta<-output$theta[iiis,]
    # Guided mean + resampled theta location:
    thth<-(theta*output$distance[iiis]+guidMu*guidW)/(output$distance[iiis]+guidW)
    # Return samples from parameter space
    thetastar<-array(rnorm(NN*lennie, thth, COVhat$cov),dim=c(NN,lennie))
    # Prepare to modify the weights by the jump from theta to thetastar
    normdist<-pnorm(thetastar,thth, COVhat$cov)
    # Make sure no values are more than 0.5 - using lower tail or higher tail
    normdist[normdist>0.5]<-normdist[normdist>0.5]-0.5
    
    return(list(theta=thetastar,modWeights=normdist))
  }
  # Return the function as output
  return(ResampleSIR)
}

# MVN based on resampled theta with a local covariance matrix based on N nearest neighbours
pert_LocalNNCov<-function(output,SumStats,weightings=NULL){
  
  
  
  pNNs
  
  
  # Make sure only successful samples contribute to the perturbation
  inds<-output$distance<output$delta
  # Calculate the global weighted mean and covariance matrix
  COVhat<-cov.wt(t(output$shat[,inds]),weightings); lennie<-ncol(output$theta[inds,])
  # Then sample from the rnorm
  ResampleSIR<-function(NN) {
    # First sample particles based on the distance function
    theta<-output$theta[sample((1:lennie)[inds],NN,T,prob = output$distance[inds]),]
    # Return samples from parameter space
    thetastar<-array(rnorm(NN*lennie, theta, COVhat$cov),dim=c(NN,lennie))
    # Prepare to modify the weights by the jump from theta to thetastar
    normdist<-pnorm(thetastar,theta, COVhat$cov)
    # Make sure no values are more than 0.5 - using lower tail or higher tail
    normdist[normdist>0.5]<-normdist[normdist>0.5]-0.5
    
    return(list(theta=thetastar,modWeights=normdist))
  }
  # Return the function as output
  return(ResampleSIR)
}













# https://arxiv.org/abs/2206.12235
# We will try using the algorithm 'fullcond', equations 9 & 10. 

# NOOOOOOOOOOOOOOOO! Use this one instead due to covariance matrix non-positive definiteness:
# Generator Parameter Calibration by Adaptive Approximate Bayesian Computation with Sequential Monte Carlo Sampler
# Seyyed Rashid Khazeiynasab Student Member, IEEE and Junjian Qi, Senior Member, IEEE, 2021
# DOI:     10.1109/TSG.2021.3077734

pert_modDistNN<-function(output,SumStats,weightings=NULL){
  
  
  
  
  
  
  # Guided N-neareat neighbours approach:
  # For sampled parameter space vector theta*
  # Take N nearest neighbours of previous iteration particles
  # then recalculate the distance function (Minkowski) for each of them
  # Then make the mean & sd of these particles, weighted by the distance,
  # Use this as the perturbation kernel.
  
  
  
  
  
  
  
  
  
  # If no weights are provided for the samples:
  if(is.null(weightings)) weightings<-rep(1,ncol(output$shat))
  # Covariance matrix
  COVhat<-cov.wt(t(output$shat),weightings)$cov
  # Check for elements that cause the covariance matrix to be non-positive-definite
  eigy<-eigen(COVhat)
  # Remove them later using inds
  inds<-eigy$values<=1e-5
  # Average summary statistic per objective function element
  sy<-rowMeans(output$shat)
  # Checks to ensure covariance matrix is positive-definite:
  # Need to ensure that any permanently zero summary stats are removed for collinearity
  inds<-inds & apply(output$shat,1,sd)!=0 & sy!=0
  
  
  
  
  
  
  # Make sure all variables with perfect collinearity are removed
  cory<-cov.wt(t(output$shat),weightings,cor = T)$cor;
  # Remove for multicollinearity
  inds<-inds & !(1:nrow(cory)%in%caret::findCorrelation(cory, cutoff=0.99))
  
  
  Stheta<-COVhat[1:12,1:12] - COVhat[1:12,13:nrow(COVhat)]%*%COVhat[13:nrow(COVhat),13:nrow(COVhat)]%*%COVhat[13:nrow(COVhat),1:12]
  
  
  # Prevent zeros from impacting calculation of inverse of covariance matrix
  output$shat<-output$shat[inds,]; SumStats<-SumStats[inds]
  # Combine (theta,s) into a vector
  xxx<-cbind(output$theta,t(output$shat))
  # Calculate full (theta,s) covariance matrices
  COVhat<-cov.wt(xxx,weightings)
  # Take the mean from the weighted covariance calculation to save time
  Mhat<-COVhat$center; COVhat<-COVhat$cov
  # Store length of theta vector
  lennie<-ncol(output$theta); lsy<-nrow(output$shat)
  # This should be vectorised, but it actually won't take long compared to the particle filter
  Maddy<-rep(NA,lennie)
  Mmulty<-SigStar<-array(NA,c(lennie,lsy+lennie-1))
  for(k in 1:lennie){
    Mmulty[k,]<-c(COVhat[k,-k]%*%chol2inv(COVhat[-k,-k]))
    Maddy[k]<-Mhat[k] - c(Mmulty[k,]%*%Mhat[-k])
    SigStar[k]<-sqrt(COVhat[k,k] - c(Mmulty[k,]%*%COVhat[-k,k]))
  }
  # Then sample from the rnorm
  ResampleSIR<-function(NN) {
    # First sample particles based on the distance function
    theta<-output$theta[sample(seq_along(output$theta),NN,T,prob = output$distance),]
    # Then perturb the particles based on the location in parameter-summary statistic space
    Mstar<-Maddy+Multy%*%rbind(t(theta[,-k]), array(rep(SumStats,NN),c(lsy,NN)))
    # Return samples from parameter space
    array(rnorm(NN*lennie, Mstar, sigStar),dim=c(NN,lennie))
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  DF<-output$shat%>%t()%>%as.data.frame.array()
  namers<-paste0("V",1:ncol(DF))
  colnames(DF)<-namers
  tmp<-mclapply(seq_along(DF),function(i) {
    formy<-paste0(namers[i],"~",paste(namers[-i],collapse = "+"))
    return(mean(abs(lm(formy,DF)$residuals)))
  },mc.cores = ncores)%>%unlist()%>%c()
  
  
  
  
  
  
  OLS<-apply(output$shat,2,function(x) abs(x-c(SumStats)))%>%t()
  # Find the variance of each parameter on the summary statistics to see what most influences posterior
  standies<-sapply(seq_along(output$theta[1,]),function(i) {
    sapply(seq_along(OLS[1,]),function(j) abs(cor(output$theta[,i],OLS[,j])))
  })
  ParamVar<-apply(standies,2,median,na.rm=T)
  # Then also calculate the weight of the different summary statistics in terms of contribution to the posterior
  SSweight<-apply(standies,1,median,na.rm=T)
  
  
  
  
  # Recalculate the distances:
  distie<-Minkowski(output$shat,SumStats,rep(ncol(output$shat),3),p=1)
  # Reduce the objective function vector to remove values that didn't change
  sy<-Rfast::rowmeans(output$shat)
  # Checks to ensure covariance matrix is positive-definite:
  # Need to ensure that any permanently zero summary stats are removed for collinearity
  inds<-apply(output$shat,1,sd)!=0 & sy!=0
  # Now create a non-parametric kernel density estimate with the parameter space
  cory<-apply(OLS[,inds],2,function(x) cor(x,distie$sw))
  # Set the weights for the summary statistics 
  SSweight<-rep(0,length(inds))
  SSweight[inds]<-minmaxScale(cory)
  
  
  # Influence of the parameter on reducing all the summary statistics towards zero
  # Maybe calculate the distance and variance of the summary statistics, per position in parameter space
  # Create a kernel density estimation of the variance, but weighted by the distance function
  # Then create an optimisation technique to reduce the two?
  
  
  # Create a function, using GPR, of variance per position in parameter space, theta.
  # 
  # Create another function of mean? Use the weighted mean between the weighted mean of all particles
  # With that of the local point in parameter space: nearest neighbour?
  
  
  
  
  
  stop("insert sy not s_sample into mean equation")
  
  return(ResampleSIR)
}



