library(survival)
library(cmprsk)

# the following function computes the extended IDI based on the data on the observed event times and cause indicators

estIDI_Fine = function(X,delta, Z, t0, withsd)
{ndim = length(X)
npts = length(t0)

#############################################################################
### Fit the cox proportional hazard regression model without z3 (model 1) ###
#############################################################################

# calcualte the cumulative baseline hazard function for cause 1 #
cov = cbind(Z[,1],Z[,2],Z[,4])
crr1 = crr(X,delta,cov)
pred1 = predict(crr1,cov)
# summary(crr1)

# indecies for the predified t0 #
t0.index= c(max(which(pred1[,1]<t0[1])),max(which(pred1[,1]<t0[2])),max(which(pred1[,1]<t0[3])))

# estimate the overall survival function #
p1.hat.m1 = t(pred1[t0.index,-1])

crr2 = crr(X,delta,cov,failcode=2)
pred2 = predict(crr2,cov)
# summary(crr1)

# indecies for the predified t0 #
t0.index= c(max(which(pred2[,1]<t0[1])),max(which(pred2[,1]<t0[2])),max(which(pred2[,1]<t0[3])))

# estimate the overall survival function #
p2.hat.m1 = t(pred2[t0.index,-1])
YY=(p1.hat.m1+p2.hat.m1>1)
KK=(p1.hat.m1/(p1.hat.m1+p2.hat.m1))[p1.hat.m1+p2.hat.m1>1]
p1.hat.m1[YY]=KK
p2.hat.m1[YY]=1-KK

p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1



#######################################################################################
### Fit the cox proportional hazard regression model using all covariates (model 2) ###
#######################################################################################
# calcualte the cumulative baseline hazard function for cause 1 #
cov = cbind(Z[,1],Z[,2],Z[,3],Z[,4],Z[,6])
crr1 = crr(X,delta,cov)
pred1 = predict(crr1,cov)
# summary(crr1)

# indecies for the predified t0 #
t0.index= c(max(which(pred1[,1]<t0[1])),max(which(pred1[,1]<t0[2])),max(which(pred1[,1]<t0[3])))

# estimate the overall survival function #
p1.hat.m2 = t(pred1[t0.index,-1])

cov = cbind(Z[,1],Z[,2],Z[,3],Z[,4])
crr2 = crr(X,delta,cov,failcode=2)
pred2 = predict(crr2,cov)
# summary(crr1)

# indecies for the predified t0 #
t0.index= c(max(which(pred2[,1]<t0[1])),max(which(pred2[,1]<t0[2])),max(which(pred2[,1]<t0[3])))

# estimate the overall survival function #
p2.hat.m2 = t(pred2[t0.index,-1])
YY=(p1.hat.m2+p2.hat.m2>1)
KK=(p1.hat.m2/(p1.hat.m2+p2.hat.m2))[p1.hat.m2+p2.hat.m2>1]
p1.hat.m2[YY]=KK
p2.hat.m2[YY]=1-KK

p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2





#############################################################################
#### Kaplan-Meier estimator for the censoring distribution, G(t)=pr(C>t) ####
#############################################################################


censor.S1= NULL
censor.lambdastar = NULL
censored.index = which(delta == 0)

if (length(censored.index)>0) ## handle the case there is no censoring
{
  sort.time = sort(X)
  
  censor.tmp1 = matrix(rep(X,each = ndim),ndim,ndim)
  censor.tmp2 = matrix(rep(delta,each = ndim),ndim,ndim)
  
  censor.Y  = apply(ifelse(censor.tmp1 >= sort.time,1,0),1,sum)
  censor.d  = apply(ifelse(censor.tmp1 ==  sort.time & censor.tmp2 == 0,1,0),1,sum)
  censor.lambdastar = censor.d/censor.Y
  censor.S1 =  cumprod(1-censor.lambdastar)
  
  ### add 1 into censor.S1 and 0 to censor.uniquetime ###
  censor.S1 = c(1,censor.S1[1:(ndim-1)])
}

# position in the estimated survival and hazard for censoring for ith obs
pos = floor(rank(X))
est.S.censoring = censor.S1[pos]
est.lambda.censoring=censor.lambdastar[pos]


# first time pt
den1.t1 = ifelse(X<=t0[1] & delta ==1,1,0)/est.S.censoring
den1.t1 = ifelse(den1.t1=='NaN',0,den1.t1)


den2.t1 = ifelse(X<=t0[1] & delta ==2,1,0)/est.S.censoring
den2.t1 = ifelse(den2.t1=='NaN',0,den2.t1)

censor.t1 = max(which(sort(X) <=t0[1]))
den3.t1 = ifelse(X>t0[1],1,0)/censor.S1[censor.t1+1]
den3.t1 = ifelse(den3.t1=='NaN',0,den3.t1)


den.t1 = c(sum(den1.t1), sum(den2.t1), sum(den3.t1))
n.t1 = sum(den.t1)
pi.t1 = den.t1/ndim


pp.m2.t1 = cbind(p1.hat.m2[,1],p2.hat.m2[,1], p3.hat.m2[,1])
pp.m1.t1 = cbind(p1.hat.m1[,1],p2.hat.m1[,1], p3.hat.m1[,1])


num.t1.d1 = (pp.m1.t1 - matrix(rep(apply(pp.m1.t1,2,mean),ndim),ndim,3,byrow=TRUE))^2
num.t1.d2 = (pp.m2.t1 - matrix(rep(apply(pp.m2.t1,2,mean),ndim),ndim,3,byrow=TRUE))^2

num.t1 = apply(num.t1.d2-num.t1.d1,2,mean)
estIDI.t1 = mean(num.t1/pi.t1/(1-pi.t1))


# second time point
den1.t2 = ifelse(X<=t0[2] & delta ==1,1,0)/est.S.censoring
den1.t2 = ifelse(den1.t2=='NaN',0,den1.t2)


den2.t2 = ifelse(X<=t0[2] & delta ==2,1,0)/est.S.censoring
den2.t2 = ifelse(den2.t2=='NaN',0,den2.t2)

censor.t2 = max(which(sort(X) <=t0[2]))
den3.t2 = ifelse(X>t0[2],1,0)/censor.S1[censor.t2+1]
den3.t2 = ifelse(den3.t2=='NaN',0,den3.t2)


den.t2 = c(sum(den1.t2), sum(den2.t2), sum(den3.t2))
n.t2 = sum(den.t2)
pi.t2 = den.t2/ndim


pp.m2.t2 = cbind(p1.hat.m2[,2],p2.hat.m2[,2], p3.hat.m2[,2])
pp.m1.t2 = cbind(p1.hat.m1[,2],p2.hat.m1[,2], p3.hat.m1[,2])


num.t2.d1 = (pp.m1.t2 - matrix(rep(apply(pp.m1.t2,2,mean),ndim),ndim,3,byrow=TRUE))^2
num.t2.d2 = (pp.m2.t2 - matrix(rep(apply(pp.m2.t2,2,mean),ndim),ndim,3,byrow=TRUE))^2

num.t2 = apply(num.t2.d2-num.t2.d1,2,mean)
estIDI.t2 = mean(num.t2/pi.t2/(1-pi.t2))

# third time point
den1.t3 = ifelse(X<=t0[3] & delta ==1,1,0)/est.S.censoring
den1.t3 = ifelse(den1.t3=='NaN',0,den1.t3)


den2.t3 = ifelse(X<=t0[3] & delta ==2,1,0)/est.S.censoring
den2.t3 = ifelse(den2.t3=='NaN',0,den2.t3)

censor.t3 = max(which(sort(X) <=t0[3]))
den3.t3 = ifelse(X>t0[3],1,0)/censor.S1[censor.t3+1]
den3.t3 = ifelse(den3.t3=='NaN',0,den3.t3)


den.t3 = c(sum(den1.t3), sum(den2.t3), sum(den3.t3))
n.t3 = sum(den.t3)
pi.t3 = den.t3/ndim


pp.m2.t3 = cbind(p1.hat.m2[,3],p2.hat.m2[,3], p3.hat.m2[,3])
pp.m1.t3 = cbind(p1.hat.m1[,3],p2.hat.m1[,3], p3.hat.m1[,3])

num.t3.d1 = (pp.m1.t3 - matrix(rep(apply(pp.m1.t3,2,mean),ndim),ndim,3,byrow=TRUE))^2
num.t3.d2 = (pp.m2.t3 - matrix(rep(apply(pp.m2.t3,2,mean),ndim),ndim,3,byrow=TRUE))^2

num.t3 = apply(num.t3.d2-num.t3.d1,2,mean)
estIDI.t3 = mean(num.t3/pi.t3/(1-pi.t3))


return(c(estIDI.t1, estIDI.t2, estIDI.t3))

}
library(parallel)
DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
ndim=nrow(DATA) 
t0 = c(8,10,12)
ests=estIDI_Fine(DATA$t,DATA$indicator,DATA[,4:9], t0, TRUE)
IDI =ests
nboot=1000


gene=function(s){
  DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
  ndim=nrow(DATA) 
  set.seed(2020)
  seeds=floor(runif(1000000,1,100000000))
  set.seed(seeds[s])
  index = sample(1:ndim, ndim, replace=TRUE)
  DATA.boot = DATA[index,]
  t0 = c(8,10,12)
  return(estIDI_Fine(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estIDI_Fine","crr"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
boot.IDI=parLapply(clust,1:nboot, gene)

stopCluster(clust)


boot.IDI=matrix(unlist(boot.IDI),nrow=3)

gene=function(s){
  DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
  DATA.boot = DATA[-s,]
  t0 = c(8,10,12)
  return(estIDI_Fine(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estIDI_Fine","crr"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
jack.IDI=parLapply(clust,1:nrow(DATA), gene)

stopCluster(clust)


jack.IDI=matrix(unlist(jack.IDI),nrow=3)

z0.hat=qnorm(colMeans(t(boot.IDI)<matrix(rep(t(IDI),each=1000),ncol=3)))
theta0.hat=apply(jack.IDI,1,mean)
temp=-t(jack.IDI)+matrix(rep(theta0.hat,each=ncol(jack.IDI)),ncol=3)
a.hat=colSums(temp^3)/6/colSums(temp^2)^(1.5)
alpha1=pnorm(z0.hat+(z0.hat+qnorm(0.025))/(1-a.hat*(z0.hat+qnorm(0.025))))
alpha2=pnorm(z0.hat+(z0.hat+qnorm(0.975))/(1-a.hat*(z0.hat+qnorm(0.975))))

c(quantile(boot.IDI[1,],alpha1[1]),quantile(boot.IDI[2,],alpha1[2]),quantile(boot.IDI[3,],alpha1[3]))
c(quantile(boot.IDI[1,],alpha2[1]),quantile(boot.IDI[2,],alpha2[2]),quantile(boot.IDI[3,],alpha2[3]))
c(quantile(boot.IDI[1,],0.025),quantile(boot.IDI[2,],0.025),quantile(boot.IDI[3,],0.025))
c(quantile(boot.IDI[1,],0.975),quantile(boot.IDI[2,],0.975),quantile(boot.IDI[3,],0.975))
IDI

# 
# c(8,10,12)
# 0.012706894	0.02065968	0.026973727
# 0.024634119	0.039031737	0.049132403
# 0.013362345	0.021634538	0.028279836
# 0.026191144	0.040664884	0.051044189
# 0.018374752	0.029783626	0.038109037
# 
# 
# 0.052146364	0.055531574	0.056876383
# 0.095408578	0.097908958	0.095474739
# 0.062156771	0.065380445	0.067033944
# 0.110231356	0.114885025	0.11919159
# 0.078077816	0.081100939	0.081721068
