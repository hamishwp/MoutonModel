###############    CHECKS   ##################
# particleFilter :
if (checks){
  if (!class(mu) %in% c('character','function')) stop('invalid mu')
  if (!class(mu) %in% c('character','function')) stop('invalid sampleState')
  if (class(muPar)!='list') stop('muPar must be passed as a list')
  if (class(sampleStatePar)!='list') stop('sampleStatePar must be a list')
  if ((class(b)!='integer' & b%%1!=0) | b<=0) stop('Invalid choice of b')
  if ((class(t)!='integer' & t%%1!=0) | t<=0) stop('Invalid choice of t')
}

# sampleStateIPM :
if (checks){
  if (class(previousState)!='integer') stop('Invalid previous state')
  if (any(previousState%%1!=0)) stop('State space can contain only integers')
  if (!class(survFunc) %in% c('function','character')) stop('Invalid input')
  if (!class(growthFunc) %in% c('function','character')) stop('Invalid input')
  if (!class(reprFunc) %in% c('function','character')) stop('Invalid input')
  if (!class(offNumSamp) %in% c('function','character')) stop('Invalid input')
  if (!class(offSizeSamp) %in% c('function','character'))stop('Invalid input')
}

# sampleNorm and sampleDTN :
if (length(pars)!=6) stop('Incorrect number of parameters')
if (U<L) stop("Upper limit smaller than lower limit")
