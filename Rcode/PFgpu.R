lenny<-1000

###################################################
ptm <- proc.time()
### non-parallelised, classical PF-style method ###
x<-matrix(exp(-(grunif(100000)-0.5)^2),10,10000)
r<-matrix(0, 10, lenny)
xsam<-sample(x,prob=x)
proc.time() - ptm
###################################################


##################################################
ptm <- proc.time()
### GPU method for PF-style required functions ###
library("gmatrix")
# Uniform RNG, default values \in (0,1)
#tmp<-exp(-(grunif(100000)-0.5)^2)
# Order the vector and project as a gmatrix (order is row,column)
#x<-gmatrix(tmp[gorder(tmp)],10,10000)
g
# Sample values from each row
r<-matrix(0, 10, lenny)
for (j in 1:lenny){
  r[,j]<-as.numeric(rsample(x,log=FALSE))
}
xy<-as.matrix(x)
for (i in 1:10){xsam[i]<-xy[i,r[i,]]}
ggc()
proc.time() - ptm
##################################################