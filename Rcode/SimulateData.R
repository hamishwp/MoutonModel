library(magrittr) # for the ever useful pipe operator

linLogit <- function(z, par){
  # purpose : Calculates the probability an individual of size x survives to 
  #           the next time step. Uses a logistic form.
  # inputs  : z         - The size of the individual (continuous)
  #           par       - A vector with entries 'intercept' and 'gradient'
  # output  : The density of survival
  (par[1] + par[2]*z) %>% plogis() %>% return
}

sampleRec <- function(n, rate, L, U){
  # purpose : sample sizes from the recruitment function using a rejection 
  #           algorithm.
  # inputs  : n    - The desired number of samples
  #           rate - The rate parameter of the function
  #           L    - The lower bound of the size variable 
  #           U    - The upper bound of the size variable
  #
  # NOTE : log link on rate?
  samples <- rep(NA, n)
  for (i in 1:n){
    sampled.size <- rexp(1, rate)
    while(sampled.size<L | sampled.size>U) sampled.size <- rexp(1, rate)
    samples[i] <- sampled.size
  }
  
  return(samples)
}

sampleDTN <- function(x, pars){
  
  # range is a vector of two values
  mu<- x*pars[2] + pars[1]
  # F.a <- pnorm(pars[4], mean = mu, sd = pars[3])
  # F.b <- pnorm(pars[5], mean = mu, sd = pars[3])
  
  return(qnorm(runif(length(x), 
             min = pnorm(pars[4], mean = mu, sd = pars[3]), 
             max = pnorm(pars[5], mean = mu, sd = pars[3])),
        mean = mu, sd = pars[3]))
  
}

dgammaM<-function(x,mu,sig_percent){
  # rgamma(n shape = alpha, scale = theta)
  # Note that the sig_percent is the constant coefficient of variation
  # Therefore, it is like a percentage of the mean
  ssq<-sig_percent*sig_percent
  dgamma(x,shape=1./ssq,scale=mu*ssq)
}

# sampleDTN <- function(x, pars){
#   # purpose : produces samples from a doubly truncated normal distribution 
#   #           using a rejection sampling algorithm
#   # extract parameter values:
#   int <- pars[1] ; gra <- pars[2] ; sigma <- pars[3]
#   L <- pars[4] ; U <- pars[5] ;
#   
#   # obtain the samples:
#   mean <- x*gra + int
#   sampled.size <- rnorm(length(x), mean, sigma)
#   
#   
#   
#   # perform rejection sampling to ensure the samples are within (L, U):
#   while (any(sampled.size<L | sampled.size>U)){
#     badSamples <- sampled.size<L | sampled.size>U
#     replacements <- rnorm(sum(badSamples), mean[badSamples], sigma)
#     sampled.size[badSamples] <- replacements
#   }
#   
#   return(sampled.size)
# }



sampleNorm <- function(x, pars){
  # purpose : produces samples as DTN without the double truncation.
  # if(anyNA(rnorm(n = pars[4], mean = (x*pars[2] + pars[1]) , sd = pars[3]))) print(x)
  return(rnorm(n = length(x), mean = (x*pars[2] + pars[1]) , sd = pars[3]))
}

sampleOffNumMod <- function(z, rate) (rpois(length(z),rate-1L)+1L)

sampleOffNum <- function(z, rate){
  # purpose : Samples the number of offspring an individual has when they 
  #           reproduce, from a truncated poisson distribution
  # input   : The mean number of children each parent has when they reproduce
  # output  : The integer number of children had by the parent
  
  # if (class(rate)!="numeric") stop("invalid rate input")
  # rate <- exp(rate)
  
  # rejection algorithm that samples from the truncated Poisson with at least
  # 1 offspring produced:
  num <- rpois(length(z), rate)
  while(any(num==0)){
    badSamples <- num==0
    replacements <- rpois(sum(badSamples), rate)
    num[badSamples] <- replacements
  }
  
  return(num)
}

sampleExp <- function(z, pars){
  # purpose : Samples the size of an offspring from a truncated exponential
  #           distribution.
  # input   : The rate of the exponential distribution
  # output  : Samples of the children sizes
  
  rate <- exp(pars[1]) ; L <- pars[2] ; U <- pars[3]
  if (U<=L) stop('invalid size range')
  if (class(rate)!="numeric") stop("invalid rate input")
  if (rate<=0) stop("rate must be positive")
  
  # rejection algorithm that samples from the truncated exponential offspring
  # size distribution:
  samples <- rexp(length(z), rate)
  while(any(samples<L | samples>U)){
    badSamples <- samples < L | samples > U
    replacements <- rexp(sum(badSamples), rate)
    samples[badSamples] <- replacements
  }
  
  return(samples)
}

calcBreaks<-function(SHEEP,nbks,regbinspace=T){
  # Get breakpoints (histogram bin intervals), count data, and growth rate per census:
  if(regbinspace){
    breaks<-ggplot_build(ggplot()+geom_histogram(data=SHEEP$solveDF,mapping=aes(size),bins = nbks))$data[[1]]$x
    # Merge all bins that are in the lowest 5th percentile
    minnie<-quantile(SHEEP$solveDF$size,1/nbks,na.rm = T)%>%unname()
    breaks<-seq(from=minnie,to=max(breaks),length.out=nbks)
    # others<-seq(from=minnie,to=max(breaks),length.out=(nbks-1))
    # breaks<-c(min(breaks),others)
    if(breaks[nbks]>SHEEP$U) breaks[nbks]<-max(c(SHEEP$solveDF$size,SHEEP$solveDF$prev.size),na.rm = T)-1e-5
    if(breaks[1]>SHEEP$L) breaks[1]<-min(c(SHEEP$solveDF$size,SHEEP$solveDF$prev.size),na.rm = T)+1e-5
  } else {
    breaks<-quantile(c(SHEEP$solveDF$prev.size,SHEEP$solveDF$size),probs = 0:(nbks-1)/(nbks-1),na.rm = T)%>%unname
  }
  return(breaks)
}

GetSoaySheep_binned <-function(SHEEP,shift=0.5,oneSex=T,nbks=10,regbinspace=F){
  
  SHEEP$breaks<-calcBreaks(SHEEP,nbks,regbinspace)
  
  SHEEP$sizes <- SHEEP$breaks[-(nbks)] + shift*diff(SHEEP$breaks)
  breaks<-SHEEP$breaks; breaks[c(1,nbks)]<-c(-Inf,Inf)
  SHEEP$COUNTS <- getSizeDistns(SHEEP$solveDF, SHEEP$breaks)
  
  SHEEP$priorProbs<-rowSums(SHEEP$COUNTS)/sum(SHEEP$COUNTS)
  SHEEP$obsProbTime <- apply(SHEEP$COUNTS, 2, sum)/SHEEP$detectedNum
  
  # What form? array? We have size bins, variable and year, so 3D matrix?
  
  formData<-function(SHEEP, censy){
    # Filter by census first
    solveDF<-SHEEP$solveDF[SHEEP$solveDF$census.number==censy,]
    D<-length(SHEEP$breaks)
    # Initialise output array - dimension are c(size bins, validation variables)
    output<-array(NA,dim = c(D-1,3))
    # Sheep that survived from last census, based on current size bin
    output[,1]<-vectorToCounts(c(solveDF$size[solveDF$survived==1]), 
                               SHEEP$breaks) #- output[,1]
    # Sheep that were also parents in this census
    output[,2]<-vectorToCounts(c(solveDF$size[solveDF$reproduced==1]), SHEEP$breaks)
    # Number of offspring per size bin
    output[,3]<-vectorToCounts(c(solveDF$rec1.wt,solveDF$rec2.wt), SHEEP$breaks)
    
    return(output)
  }
  
  cen<-unique(SHEEP$solveDF$census.number)
  
  SHEEP$SumStats<-sapply(cen,
                         FUN = function(i) formData(SHEEP,i),
                         simplify = F) %>% unlist() %>%
    array(dim=c(nbks-1,3,length(cen)))
  
  return(SHEEP)
  
}

simulateIBM <- function(n, t, survFunc, survPars, growthSamp, growthPars,
                        reprFunc, reprPars, offNum=NULL, offNumSamp=NULL,
                        offNumPars=NULL, offSizeSamp, offSizePars, Start,
                        Schild=1, obsProb=NULL, OneGend=FALSE, thresh=Inf,
                        CalcOffSurv=TRUE, popPrint=TRUE, verbose=FALSE){
  # purpose : Uses an Individual Based Model (IBM) to simulate census data used
  #           in fitting Integral Projection Models.
  # inputs  : n           - The starting size of the population
  #           t           - The maximum number of censuses to simulate
  #           survFunc    - The function which determines the probability of 
  #                         survival. It takes two arguments, the size of the 
  #                         individual, and a vector of parameters.
  #           survPar     - The vector of parameters for the survival function
  #           growthSamp  - The function which samples new sizes for
  #                         individuals. It takes two arguments, the size of the 
  #                         individual, and a vector of parameters
  #           growthPars  - The parameters for the function which samples the
  #                         new size of individuals.
  #           reprFunc    - The function which determines the probability of 
  #                         reproduction. It takes two arguments, the size of
  #                         the individual, and a vector of parameters.
  #           reprPars    - The parmeters for the reproduction function.
  #           offNum      - If not NULL, indicates the fixed number of 
  #                         offspring that a parent has when they reproduce.
  #           offNumSamp  - If offNum is NULL, this is the function which
  #                         samples the number of offspring that a parent has
  #                         when they reproduce. It takes two arguments, the
  #                         size of the individual, and a vector of parameters
  #           offNumPars  - If offNum is NULL, this is the vector of parameters
  #                         for the number of offspring each parent has when 
  #                         they reproduce.
  #           offSizeSamp - The function which samples the size of the offspring
  #                         when they are first censused. It takes two
  #                         arguments, the size of the parent, and a vector of
  #                         parameters.
  #           offSizePars - The vector of parameters for the function which
  #                         samples the sizes of offspring.
  #           Start       - The size of the individual which fathers all the 
  #                         initial members of the population. Can be a vector,
  #                         if so, it must be of length n
  #           Schild      - The probability of survival of children to their 
  #                         first census.
  #           obsProb     - The probability of the individual to be observed at 
  #                         in any one census.
  #           OneGend     - If TRUE, the census only tracks one gender, and so 
  #                         only half of the children born are assimilated into
  #                         the census data.
  #           thresh      - If the size of the population ever exceeds this 
  #                         quantity, we stop simulating.
  #           CalcOffSurv - If TRUE will calculate the number of offspring
  #                         that recruited for each parent in the population.
  #                         if FALSE, the function will fill this column with 
  #                         all NAs to save on computation time.
  #           popPrint    - If TRUE, the size of the population is printed at
  #                         the end of each census.
  #           verbose     - if TRUE, will print the survival rate, reproduction
  #                         rate, avergae change in size, rate of births and
  #                         child cenus rate at every survey:
  
  # make sure the number of children per litter is specified correctly:
  if (is.null(offNum) & (is.null(offNumSamp) | is.null(offNumPars))){
    stop('offNum or both offNumFunc and offNumPars must be specified')
  }
  
  if (is.null(thresh)) thresh <- .Machine$integer.max
  
  # to allow names ot be passed:
  survFunc %<>% match.fun ; growthSamp %<>% match.fun ; reprFunc %<>% match.fun
  offSizeSamp %<>% match.fun
  
  # growthPars[3]<-exp(growthPars[3])
  # offSizePars[3]<-exp(offSizePars[3])
  
  # If a parent always has the same number of children, we make sure it is 
  # an integer:
  # if (!is.null(offNum)) offNum <- ceiling(offNum)
  
  # generate first generation from parent of fixed size:
  if (length(Start)!=n) {
    
    Start <- rep(Start[1], n)
    initial.sizes <- offSizeSamp(Start, offSizePars)
    
    DF <- data.frame(id=1:n, size=initial.sizes, survived=1, 
                     census.number=1, parent.size=Start, reproduced=NA, 
                     rec1.wt=NA, rec2.wt=NA,
                     prev.size=NA, off.survived=NA, off.born=NA)  
    
  } else DF <- data.frame(id=1:n, size=Start, survived=rep(1,length(Start)), 
                          census.number=rep(1,length(Start)), 
                          parent.size=sample(Start,length(Start),replace = T), 
                          reproduced=NA, rec1.wt=NA, rec2.wt=NA,
                          prev.size=NA, off.survived=NA, off.born=NA)
  
  # create a data.frame which contains only the individuals from the previous
  # census (and will be constantly updated):
  previous.census <- DF
  
  time <- 2 # set the current census number (will be continually updated)
  while((time < t + 2) & nrow(previous.census)<thresh){
    
    # Select only the survivors of the previous census:
    survivorsDF <- subset(previous.census, previous.census$survived==1)
    current.pop.size <- nrow(survivorsDF)
    
    # Determine who survives this time:
    survivors<-newSizes<-rep(NA,current.pop.size)
    survivors[!is.na(survivorsDF$size)] <- rbinom(length(survivorsDF$size[!is.na(survivorsDF$size)]), 1, 
                survFunc(survivorsDF$size[!is.na(survivorsDF$size)], survPars))
    
    if (sum(survivors,na.rm = T)==0) break
    
    iids<-survivors==1 & !is.na(survivors)
    
    # Determine their new sizes:
    newSizes[iids] <- growthSamp(survivorsDF$size[iids],growthPars)

    # Determine who reproduces (out of the survivors):
    reproducers <- rep(NA, current.pop.size)
    # reproducers[iids] <- rbinom(sum(survivors,na.rm = T), 1, reprFunc(na.omit(newSizes), reprPars))
    reproducers[iids] <- rbinom(sum(survivors,na.rm = T), 1, reprFunc(survivorsDF$size[iids], reprPars))
    
    # Determine the sex of the offspring:
    probGend <- ifelse(isTRUE(OneGend), 0.5, 1)
    if(sum(reproducers==1 & !is.na(reproducers))>0) reproducers[reproducers==1 & !is.na(reproducers)]<-
      rbinom(sum(reproducers==1 & !is.na(reproducers)), 1, probGend)
    
    iids<-reproducers==1 & !is.na(reproducers)
    
    # nRepr <- reproducers %>% na.omit %>% `==`(1) %>% sum
    parentsDF <- survivorsDF[iids,]
    
    # Determine the number of offspring for each parent:
    if (!is.null(offNumPars)) {censusOffNum <-offNumSamp(reproducers[iids], offNumPars)
    }else if(!is.null(offNum)) {censusOffNum <- rep(reproducers[iids],offNum)
    }else stop("No offspring parameters provided, check offNum or offNumPars")
    # parentSizes <- offNumSamp(reproducers, offNumPars)

    # Determine the sizes of the offspring:
    reprSizes <- survivorsDF$size[iids]
    parentSizes <- rep(reprSizes, censusOffNum)
    parentIDs <- rep(parentsDF$id, censusOffNum)
    offspringSizes <- offSizeSamp(parentSizes, offSizePars)
    
    # Determine the survival of the offspring:
    nOffspring <- length(offspringSizes)
    survivingChildren <- rbinom(nOffspring, 1, Schild)
    
    # Find out which children make it to census:
    newIDStart <- (max(DF$id)+1)
    censusedChildren <- which(survivingChildren==1)

    # Calculate how many surviving children per parent with a helper:
    helper <- function(x){
      ifelse(!x %in% parentIDs,
             NA,
             sum(x==parentIDs & survivingChildren==1))
    }#helper

    helper2 <- function(x){
      ifelse(!x %in% parentIDs,
             NA,
             sum(x==parentIDs))
    }#helper2

    helper3<-function(x){
      
      if(x %in% parentIDs){
        # Which parents had offspring?
        tmp<-offspringSizes[parentIDs==x]
        if(length(tmp)>=2) {
          return(c(rec1.wt=tmp[1],rec2.wt=tmp[2]))
        } else if(length(tmp)==1){
          return(c(rec1.wt=tmp,rec2.wt=NA))
        } else stop("Something wrong with number of offspring in simulated data")
          
      } else return(c(rec1.wt=NA,rec2.wt=NA))
      
    }
    # Only use the helper if the user wants us to. Otherwise make a vector of
    # NAs for speed of calculation:
    if (isTRUE(CalcOffSurv)){
      # off.born <- rep(0,length(survivorsDF$id))
      # tably<-c(table(parentIDs))
      # off.born[survivorsDF$id%in%names(tably)]<-tably
      # # 
      # off.survived <- rep(0,length(survivorsDF$id))
      # off.survived[off.born==0]<-NA
      # tably<-c(table(parentIDs[survivingChildren==1]))
      # off.survived[survivorsDF$id%in%names(tably)]<-unname(tably)
      #       
      off.born<-sapply(survivorsDF$id,helper2)
      off.survived<-sapply(survivorsDF$id,helper)
      rec.wt<- t(sapply(survivorsDF$id,helper3))
    } else{
      off.survived <- rep(NA, current.pop.size)
      off.born <- rep(NA, current.pop.size)
      rec.wt <- matrix(NA, nrow = current.pop.size, ncol = 2)
    }
    
    # Print out summaries of the simulation if desired:
    if (verbose){
      cat("survival rate is", sum(survivors)/current.pop.size,"\n")
      oldSizes <- survivorsDF$size[which(survivors==1)]
      cat("average growth is", mean(na.omit(newSizes) - oldSizes), "\n")
      cat("rate of births is", nOffspring/current.pop.size, "\n")
      cat("no. offspring is", mean(off.born[off.born>0],na.rm=T), "\n")
      cat("no. surviving offspring is", mean(off.survived/off.born,na.rm=T), "\n")
      cat("child censusing rate is", length(censusedChildren)/nOffspring,"\n")
      gr <- (length(censusedChildren)+sum(survivors))/current.pop.size
      cat("growth rate is", gr, "\n")
      cat("\n")
    }
    
    # Create the DF for the parents:
    currentDF <- data.frame(id=survivorsDF$id,
                            size=newSizes,
                            survived=survivors,
                            census.number=time,
                            parent.size=NA,
                            reproduced=reproducers,
                            rec1.wt=rec.wt[,1],
                            rec2.wt=rec.wt[,2],
                            prev.size=survivorsDF$size,
                            off.survived=off.survived,
                            off.born=off.born)
    
    # Update the previous.census data.frame:
    if (length(censusedChildren)==0) {
      previous.census <- currentDF
    } else {
      newIDs <- newIDStart:(length(censusedChildren)+newIDStart-1)
      offspringDF <- data.frame(id = newIDs,
                                size = offspringSizes[censusedChildren],
                                survived = 1, census.number = time,
                                parent.size = parentSizes[censusedChildren],
                                reproduced = NA, rec1.wt=NA, rec2.wt=NA, prev.size = NA,
                                off.survived = NA, off.born = NA)
      
      # Update the previous.census data.frame:
      previous.census <- rbind(currentDF, offspringDF)
    }
    
    # Update the overall data:
    DF <- rbind(DF, previous.census)
    
    # Iterate the census.number:
    time <- time + 1
    
    # Let the user know the size of the population:
    if(popPrint) print(nrow(previous.census))
    
  }
  
  # If we need to remove some of the observations
  if(!is.null(obsProb)){
    if (length(obsProb)>1){
      for(yryr in unique(DF$census.number)){
        yrID<-which(DF$census.number==yryr)
        DF[sample(yrID,size = round((1-obsProb[yryr])*length(yrID)),replace = F),-c(1,3,4)]<-NA      
      }
    } else {
      DF[sample(1:nrow(DF),size = round((1-obsProb)*nrow(DF)),replace = F),-c(1,3,4)]<-NA
    }
  } 
  
  DF%<>%filter(census.number!=1); DF$census.number<-DF$census.number-1
  
  detectedLiveFemales <- DF%>%filter(survived==1)%>%group_by(census.number)%>%
    summarise(detectedNum=length(size),.groups = 'drop_last')%>%pull(detectedNum)%>%unname()
  
  return(list(solveDF=DF,
              detectedNum=detectedLiveFemales,
              L=min(c(DF$size,DF$prev.size),na.rm = T),
              U=max(c(DF$size,DF$prev.size),na.rm = T)))

}
















simulateIBM_old <- function(n, t, survFunc, survPars, growthSamp, growthPars,
                        reprFunc, reprPars, offNum=NULL, offNumSamp=NULL,
                        offNumPars=NULL, offSizeSamp, offSizePars, Start,
                        Schild=1, obsProb=0.85, OneGend=FALSE, thresh=Inf,
                        CalcOffSurv=TRUE, popPrint=TRUE, verbose=FALSE){
  # purpose : Uses an Individual Based Model (IBM) to simulate census data used
  #           in fitting Integral Projection Models.
  # inputs  : n           - The starting size of the population
  #           t           - The maximum number of censuses to simulate
  #           survFunc    - The function which determines the probability of 
  #                         survival. It takes two arguments, the size of the 
  #                         individual, and a vector of parameters.
  #           survPar     - The vector of parameters for the survival function
  #           growthSamp  - The function which samples new sizes for
  #                         individuals. It takes two arguments, the size of the 
  #                         individual, and a vector of parameters
  #           growthPars  - The parameters for the function which samples the
  #                         new size of individuals.
  #           reprFunc    - The function which determines the probability of 
  #                         reproduction. It takes two arguments, the size of
  #                         the individual, and a vector of parameters.
  #           reprPars    - The parmeters for the reproduction function.
  #           offNum      - If not NULL, indicates the fixed number of 
  #                         offspring that a parent has when they reproduce.
  #           offNumSamp  - If offNum is NULL, this is the function which
  #                         samples the number of offspring that a parent has
  #                         when they reproduce. It takes two arguments, the
  #                         size of the individual, and a vector of parameters
  #           offNumPars  - If offNum is NULL, this is the vector of parameters
  #                         for the number of offspring each parent has when 
  #                         they reproduce.
  #           offSizeSamp - The function which samples the size of the offspring
  #                         when they are first censused. It takes two
  #                         arguments, the size of the parent, and a vector of
  #                         parameters.
  #           offSizePars - The vector of parameters for the function which
  #                         samples the sizes of offspring.
  #           Start       - The size of the individual which fathers all the 
  #                         initial members of the population. Can be a vector,
  #                         if so, it must be of length n
  #           Schild      - The probability of survival of children to their 
  #                         first census.
  #           obsProb     - The probability of the individual to be observed at 
  #                         in any one census.
  #           OneGend     - If TRUE, the census only tracks one gender, and so 
  #                         only half of the children born are assimilated into
  #                         the census data.
  #           thresh      - If the size of the population ever exceeds this 
  #                         quantity, we stop simulating.
  #           CalcOffSurv - If TRUE will calculate the number of offspring
  #                         that recruited for each parent in the population.
  #                         if FALSE, the function will fill this column with 
  #                         all NAs to save on computation time.
  #           popPrint    - If TRUE, the size of the population is printed at
  #                         the end of each census.
  #           verbose     - if TRUE, will print the survival rate, reproduction
  #                         rate, avergae change in size, rate of births and
  #                         child cenus rate at every survey:
  
  # make sure the number of children per litter is specified correctly:
  if (is.null(offNum) & (is.null(offNumSamp) | is.null(offNumPars))){
    stop('offNum or both offNumFunc and offNumPars must be specified')
  }
  
  if (is.null(thresh)) thresh <- .Machine$integer.max
  
  # to allow names ot be passed:
  survFunc %<>% match.fun ; growthSamp %<>% match.fun ; reprFunc %<>% match.fun
  offSizeSamp %<>% match.fun
  
  # growthPars[3]<-exp(growthPars[3])
  # offSizePars[3]<-exp(offSizePars[3])
  
  # If a parent always has the same number of children, we make sure it is 
  # an integer:
  # if (!is.null(offNum)) offNum <- ceiling(offNum)
  
  # generate first generation from parent of fixed size:
  if (length(Start)!=n) Start <- rep(Start[1], n)
  initial.sizes <- offSizeSamp(Start, offSizePars)
  
  # generate the data.frame to store simulated data:
  DF <- data.frame(individual=1:n, size=initial.sizes, survived=1, 
                   census.number=1, parent.size=Start, reproduced=NA, 
                   prev.size=NA, off.survived=NA, off.born=NA)
  
  # create a data.frame which contains only the individuals from the previous
  # census (and will be constantly updated):
  previous.census <- DF
  
  time <- 2 # set the current census number (will be continually updated)
  while((time < t + 1) & nrow(previous.census)<thresh){
    
    # Select only the survivors of the previous census:
    survivorsDF <- subset(previous.census, previous.census$survived==1)
    current.pop.size <- nrow(survivorsDF)
    
    # Determine who survives this time:
    survProbs <- survFunc(survivorsDF$size, survPars)
    survivors <- rbinom(current.pop.size, 1, survProbs)
    if (sum(survivors)==0) break
    
    # Determine their new sizes:
    newSizes <- growthSamp(survivorsDF$size, growthPars)
    newSizes[-which(survivors==1)] <- NA # overwrite for dead individuals
    
    # Determine who reproduces (out of the survivors):
    reproducers <- rep(NA, current.pop.size)
    reprProbs <- reprFunc(na.omit(newSizes), reprPars)
    reproducers[survivors==1] <- rbinom(sum(survivors==1), 1, reprProbs)
    nRepr <- reproducers %>% na.omit %>% `==`(1) %>% sum
    parentsDF <- subset(survivorsDF, reproducers==1 & !is.na(reproducers))
    
    # Determine the number of offspring for each parent:
    if (!is.null(offNum)) censusOffNum <- rep(offNum, nRepr)
    else censusOffNum <-(rpois(nRepr,offNumPars-1L)+1L)
    # else censusOffNum <- offNumSamp(parentsDF$size, offNumPars)
    
    # Determine the sizes of the offspring:
    reprSizes <- survivorsDF$size[reproducers==1 & !is.na(reproducers)]
    parentSizes <- rep(reprSizes, censusOffNum)
    parentIDs <- rep(parentsDF$id, censusOffNum)
    offspringSizes <- offSizeSamp(parentSizes, offSizePars)
    
    # Determine the survival of the offspring:
    nOffspring <- length(offspringSizes)
    survivingChildren <- rbinom(nOffspring, 1, Schild)
    
    # Determine the sex of the offspring:
    probGend <- ifelse(isTRUE(OneGend), 0.5, 1)
    gender <- rbinom(nOffspring, 1, probGend)
    
    # Find out which children make it to census:
    newIDStart <- (max(DF$id)+1)
    censusedChildren <- which(gender==1 & survivingChildren==1)
    
    # Print out summaries of the simulation if desired:
    if (verbose){
      cat("survival rate is", sum(survivors)/current.pop.size,"\n")
      oldSizes <- survivorsDF$size[which(survivors==1)]
      cat("average growth is", mean(na.omit(newSizes) - oldSizes), "\n")
      cat("reproduction rate is", nRepr/current.pop.size, "\n")
      cat("rate of births is", nOffspring/current.pop.size, "\n")
      cat("child censusing rate is", length(censusedChildren)/nOffspring,"\n")
      gr <- (length(censusedChildren)+sum(survivors))/current.pop.size
      cat("growth rate is", gr, "\n")
      cat("\n")
    }
    
    # Calculate how many surviving children per parent with a helper:
    helper <- function(x){
      ifelse(!x %in% parentIDs,
             NA,
             sum(x==parentIDs & survivingChildren==1))
    }#helper
    
    helper2 <- function(x){
      ifelse(!x %in% parentIDs,
             NA,
             sum(x==parentIDs))
    }#helper2
    
    # Only use the helper if the user wants us to. Otherwise make a vector of
    # NAs for speed of calculation:
    if (isTRUE(CalcOffSurv)){
      off.survived <- sapply(survivorsDF$id, helper)
      off.born <- sapply(survivorsDF$id, helper2)
    } else{
      off.survived <- rep(NA, current.pop.size)
      off.born <- rep(NA, current.pop.size)
    }
    
    # Create the DF for the parents:
    currentDF <- data.frame(id=survivorsDF$id,
                            size=newSizes,
                            survived=survivors,
                            census.number=time,
                            parent.size=NA,
                            reproduced=reproducers,
                            prev.size=survivorsDF$size,
                            off.survived=off.survived,
                            off.born=off.born)
    
    # Update the previous.census data.frame:
    if (length(censusedChildren)==0) previous.census <- currentDF
    
    else{
      newIDs <- newIDStart:(length(censusedChildren)+newIDStart-1)
      offspringDF <- data.frame(id = newIDs,
                                size = offspringSizes[censusedChildren],
                                survived = 1, census.number = time,
                                parent.size = parentSizes[censusedChildren],
                                reproduced = NA, prev.size = NA,
                                off.survived = NA, off.born = NA)
      
      # Update the previous.census data.frame:
      previous.census <- rbind(currentDF, offspringDF)
    }
    
    # Update the overall data:
    DF <- rbind(DF, previous.census)
    
    # Iterate the census.number:
    time <- time + 1
    
    # Let the user know the size of the population:
    if(popPrint) print(nrow(previous.census))
  }
  
  # NotObs<-c()
  # for (c in unique(DF$census.number)){
  #   indy<-DF$census.number==c
  #   indy<-(1:length(indy))[indy]
  #   lenS<-length(indy)
  #   subprob<-max(c(min(c(rnorm(n = 1,mean = (1-obsProb),sd = 0.05),0.99)),0.01))
  #   subs<-rbinom(n = 1,size = lenS,prob = subprob)
  #   NotObs%<>%c(sample(x = indy,size = subs,replace = F))
  # }
  
  # DF[NotObs,-c(1,4)]<-NA
  
  return(DF)
}


# UNCOMMENT the below sections for example datasets to be generated:

# Try to simulate with similar values as the MSc Thesis from Duke:
# set.seed(102938)
# testSim <- simulateIBM(n=50, t=100,
#                        # set survival details:
#                        survFunc = linLogit, survPars = c(2.26, 0.23),
#                        # set growth details:
#                        growthSamp = sampleDTN,
#                        growthPars = c(0, 1, log(0.5), 0, 10),
#                        # set reproduction details:
#                        reprFunc = linLogit, reprPars = c(-3.5, 0.25),
#                        # set offspring number and size distribution details:
#                        offNumSamp = sampleOffNum, offNumPars = log(4),
#                        offSizeSamp = sampleExp,
#                        offSizePars = c(log(10), 0, 10),
#                        Schild=0.873,
#                        # set other miscelaneous parameters:
#                        Start=3, thresh=1000, OneGend = TRUE, popPrint = FALSE)
# 
# max.cens <- testSim %>% `$`(census.number) %>% max
# subset(testSim,  testSim$census.number==max.cens) %>% nrow %>% print

# # Look at the last census to give an idea on how the population grows:
# statement <- testSim$census.number==tsteps & testSim$survived==1
# statement %>% subset(x=testSim) %>% nrow %>% print

# Try to simulate Sheep data from 2.6 of Ellner, Childs & Rees:
# simmedData <- simulateIBM(n=500, t=1000,
#                           # set survival details:
#                           survFunc = linLogit, survPars = c(-9.65, 3.77),
#                           # set growth details:
#                           growthSamp = sampleDTN,
#                           growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
#                           # set reproduction details:
#                           reprFunc = linLogit, reprPars = c(-7.23, 2.6),
#                           # set offspring number and size distribution details:
#                           offNum=1, offSizeSamp = sampleDTN,
#                           offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
#                           # Child survival probability:
#                           Schild=0.873,
#                           # set other miscelaneous parameters:
#                           Start=3, thresh=5000, OneGend = TRUE, popPrint = T)
