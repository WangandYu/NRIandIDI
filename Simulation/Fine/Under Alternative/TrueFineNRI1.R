
library(survival)
library(cmprsk)

# sample size
ndim = 1000
nsim = 1000

### true parameters ###
beta11=0.2
beta12=-0.5
beta13=1
p=0.65

T0 = c(20,21,22)

for(m in 1:nsim)
{
  ### generate the covariates ###
  
  z1=rnorm(ndim)
  z1=pmax(pmin(z1,3.5),-3.5)
  z2=ifelse(runif(ndim)<0.7,1,0)
  z3=rnorm(ndim)
  z3=pmax(pmin(z3,3.5),-3.5)
  
  F1=1-(1-p)^(exp(z1*beta11+z2*beta12+z3*beta13))
  
  W=runif(ndim)
  epsilon=ifelse(W<F1,1,2)
  TT=rep(0,ndim)
  TT[W<F1]=-log(1-(1-(1-W[W<F1])^(1/exp(z1[W<F1]*beta11+z2[W<F1]*beta12+z3[W<F1]*beta13)))/p)
  zz=exp(z1[W>=F1]*beta11/10+z2[W>=F1]*beta12/10+z3[W>=F1]*beta13/10)
  TT[W>=F1]=sapply(zz,function(o) return(rexp(1,o)))
  TT=TT^0.2
  TT=pmin(TT*20,100)
  hist(TT)
  
  for(t0 in T0){
    n1 = sum(ifelse(TT<=t0 & epsilon==1,1,0))
    n2 = sum(ifelse(TT<=t0 & epsilon==2,1,0))
    n3 = ndim - n1 - n2
    
    write.table(t(c(n1,n2,n3)),paste("NRI_n_",t0,".txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE,append = TRUE)
  }
  #############################################################################
  ### Fit the cox proportional hazard regression model without z3 (model 1) ###
  #############################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  cov = cbind(z1,z2)
  crr1 = crr(TT,epsilon,cov)
  pred1 = predict(crr1,cov)
  # summary(crr1)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred1[,1]<T0[1])),max(which(pred1[,1]<T0[2])),max(which(pred1[,1]<T0[3])))
  
  # estimate the overall survival function #
  p1.hat.m1 = (pred1[t0.index,-1])
  
  crr2 = crr(TT,epsilon,cov,failcode=2)
  pred2 = predict(crr2,cov)
  # summary(crr2)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred2[,1]<T0[1])),max(which(pred2[,1]<T0[2])),max(which(pred2[,1]<T0[3])))
  
  # estimate the overall survival function #
  p2.hat.m1 = (pred2[t0.index,-1])
  
  p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
  pind.m1=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m1[i,],p2.hat.m1[i,], p3.hat.m1[i,])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m1 = cbind(pind.m1,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  
  #######################################################################################
  ### Fit the cox proportional hazard regression model using all covariates (model 2) ###
  #######################################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  cov = cbind(z1,z2,z3)
  crr1 = crr(TT,epsilon,cov,failcode=1)
  pred1 = predict(crr1,cov)
  # summary(crr1)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred1[,1]<T0[1])),max(which(pred1[,1]<T0[2])),max(which(pred1[,1]<T0[3])))
  
  # estimate the overall survival function #
  p1.hat.m2 = (pred1[t0.index,-1])
  
  crr2 = crr(TT,epsilon,cov,failcode=2)
  pred2 = predict(crr2,cov)
  # summary(crr2)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred2[,1]<T0[1])),max(which(pred2[,1]<T0[2])),max(which(pred2[,1]<T0[3])))
  
  # estimate the overall survival function #
  p2.hat.m2 = (pred2[t0.index,-1])
  
  p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
  pind.m2=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m2[i,],p2.hat.m2[i,], p3.hat.m2[i,])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m2 = cbind(pind.m2,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  
  for(i in 1:3){
    t0=T0[i]
    sum1 = sum(ifelse(pind.m2[,i]==1 & pind.m1[,i]!=1 & TT<=t0 & epsilon ==1, 1, 0)) - sum(ifelse(pind.m2[,i]!=1 & pind.m1[,i]==1 & TT<=t0 & epsilon ==1, 1, 0))
    sum2 = sum(ifelse(pind.m2[,i]==2 & pind.m1[,i]!=2 & TT<=t0 & epsilon ==2, 1, 0)) - sum(ifelse(pind.m2[,i]!=2 & pind.m1[,i]==2 & TT<=t0 & epsilon ==2, 1, 0))
    sum3 = sum(ifelse(pind.m2[,i]==3 & pind.m1[,i]!=3 & TT>t0, 1, 0)) - sum(ifelse(pind.m2[,i]!=3 & pind.m1[,i]==3 & TT>t0, 1, 0))
    
    write.table(t(c(sum1,sum2,sum3)),paste("NRI_sum_",t0,".txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE,append = TRUE)
  }
}

dir2=''
ns = read.table(paste(dir2,"NRI_n_20.txt",sep=""),header = FALSE)
dim(ns)
ns = apply(ns,2,sum)

sums = read.table(paste(dir2,"NRI_sum_20.txt",sep=""),header = FALSE)
dim(sums)
sums = apply(sums,2,sum)

sum(sums/ns)/3 
ns = read.table(paste(dir2,"NRI_n_21.txt",sep=""),header = FALSE)
dim(ns)
ns = apply(ns,2,sum)

sums = read.table(paste(dir2,"NRI_sum_21.txt",sep=""),header = FALSE)
dim(sums)
sums = apply(sums,2,sum)

sum(sums/ns)/3 
ns = read.table(paste(dir2,"NRI_n_22.txt",sep=""),header = FALSE)
dim(ns)
ns = apply(ns,2,sum)

sums = read.table(paste(dir2,"NRI_sum_22.txt",sep=""),header = FALSE)
dim(sums)
sums = apply(sums,2,sum)

sum(sums/ns)/3 