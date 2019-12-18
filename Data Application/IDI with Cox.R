library(survival)

# the following function computes the extended IDI based on the data on the observed event times and cause indicators

estIDI_Weibull = function(X,delta, Z, t0, withsd)
{ndim = length(X)
npts = length(t0)
z1=Z[,1]
z2=Z[,2]
z3=Z[,3]
z4=Z[,4]
z5=Z[,6]
#############################################################################
### Fit the cox proportional hazard regression model without z3 (model 1) ###
#############################################################################

# calcualte the cumulative baseline hazard function for cause 1 #
data1 = data.frame(cbind(X, event=ifelse(delta==1,1,0), z1=z1, z2=z2, z3=z3, z4=z4,z5=z5))
cox1 = coxph(Surv(X, event) ~ z1  + z2+z4, data1)	
# summary(cox1)
cox.basehaz1 = basehaz(cox1, centered=FALSE)
cox.basehaz1 = cox.basehaz1[cox.basehaz1$time<=t0[npts],]	# only keep the baseline hazard until the last point of t0

# Then calculate the comulative baseline hazard fuction for cause 2 #
data2 = data.frame(cbind(X, event=ifelse(delta==2,1,0), z1=z1, z2=z2, z3=z3, z4=z4))
cox2 = coxph(Surv(X, event) ~ z1  + z2+z4 , data2)	
# summary(cox2)
cox.basehaz2 = basehaz(cox2, centered=FALSE)
cox.basehaz2 = cox.basehaz2[cox.basehaz2$time<=t0[npts],]	# only keep the baseline hazard until the last point of t0

# indecies for the predified t0 #
t0.index1 = max(which(cox.basehaz1$time<=t0[1]))
t0.index2 = max(which(cox.basehaz1$time<=t0[2]))
t0.index3 = max(which(cox.basehaz1$time<=t0[3]))
t0.index = c(t0.index1, t0.index2, t0.index3)



dim.uniquetime = length(which(cox.basehaz1$time<=t0[npts]))	


# estimate the cumulative cause specific hazard for cause 1 #
cox.haz1 = matrix(rep(cox.basehaz1$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z1,z2,z4)%*%cox1$coef),dim.uniquetime),ndim,dim.uniquetime)
# calcualte the cumulative cause specific hazard for cause 2 #
cox.haz2 = matrix(rep(cox.basehaz2$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z1,z2,z4)%*%cox2$coef),dim.uniquetime),ndim,dim.uniquetime)

# estimate the overall survival function #
S = exp(-cox.haz1-cox.haz2)
S = cbind(1,S[,1:(dim(S)[2]-1)])		# S(t-;z)
# dLemda1(t;z)
dLamda1 = cox.haz1 - cbind(0,cox.haz1[,1:(dim(cox.haz1)[2]-1)])
# estimtate cumulative incidence function for cause1 #
F1 = t(apply(S*dLamda1,1,cumsum))
# estimates of F1 at t0 #
p1.hat.m1 = F1[,t0.index]

dLamda2 = cox.haz2 - cbind(0,cox.haz2[,1:(dim(cox.haz2)[2]-1)])
# estimtate cumulative incidence function for cause1 #
F2 = t(apply(S*dLamda2,1,cumsum))
# estimates of F1 at t0 #
p2.hat.m1 = F2[,t0.index]
YY=(p1.hat.m1+p2.hat.m1>1)
KK=(p1.hat.m1/(p1.hat.m1+p2.hat.m1))[p1.hat.m1+p2.hat.m1>1]
p1.hat.m1[YY]=KK
p2.hat.m1[YY]=1-KK

p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1





#######################################################################################
### Fit the cox proportional hazard regression model using all covariates (model 2) ###
#######################################################################################

# calcualte the cumulative baseline hazard function for cause 1 #
# data1 = data.frame(cbind(X, event=ifelse(delta==1,1,0), z1=z1, z2=z2, z3=z3, z4=z4, z5=z5))
cox1 = coxph(Surv(X, event) ~ z1  + z2+ z3+z4+z5, data1)	
#summary(cox1)
cox.basehaz1 = basehaz(cox1, centered=FALSE)
cox.basehaz1 = cox.basehaz1[cox.basehaz1$time<=t0[npts],]	# only keep the baseline hazard until time t0

# Then calculate the comulative baseline hazard fuction for cause 2 #
# data2 = data.frame(cbind(X, event=ifelse(delta==2,1,0), z1=z1, z2=z2, z3=z3, z4=z4, z5=z5))
cox2 = coxph(Surv(X, event) ~ z1  + z2+ z3+z4, data2)	
#summary(cox2)
cox.basehaz2 = basehaz(cox2, centered=FALSE)
cox.basehaz2 = cox.basehaz2[cox.basehaz2$time<=t0[npts],]	# only keep the baseline hazard until time t0

# indecies for the predified t0 #
t0.index1 = max(which(cox.basehaz1$time<=t0[1]))
t0.index2 = max(which(cox.basehaz1$time<=t0[2]))
t0.index3 = max(which(cox.basehaz1$time<=t0[3]))
t0.index = c(t0.index1, t0.index2, t0.index3)

dim.uniquetime = length(which(cox.basehaz1$time<=t0[npts]))

# estimate the cumulative cause specific hazard for cause 1 #
cox.haz1 = matrix(rep(cox.basehaz1$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z1,z2,z3,z4,z5)%*%cox1$coef),dim.uniquetime),ndim,dim.uniquetime)

# calcualte the cumulative cause specific hazard for cause 2 #
cox.haz2 = matrix(rep(cox.basehaz2$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z1,z2,z3,z4)%*%cox2$coef),dim.uniquetime),ndim,dim.uniquetime)

# estimate the overall survial fuction #
S = exp(-cox.haz1-cox.haz2)
S = cbind(1,S[,1:(dim(S)[2]-1)])		# S(t-;z)
# dLamda1(t;z)
dLamda1 = cox.haz1 - cbind(0,cox.haz1[,1:(dim(cox.haz1)[2]-1)])
# esimtate cumulative incidence function for cause1 #
F1 = t(apply(S*dLamda1,1,cumsum))
# estimates of F1 at t0 #
p1.hat.m2 = F1[,t0.index]

# dLamda2(t;z)
dLamda2 = cox.haz2 - cbind(0,cox.haz2[,1:(dim(cox.haz2)[2]-1)])
# esimtate cumulative incidence function for cause1 #
F2 = t(apply(S*dLamda2,1,cumsum))

# estimates of F2 at t0 #
p2.hat.m2 = F2[,t0.index]
YY=(p1.hat.m2+p2.hat.m2>1)
KK=(p1.hat.m2/(p1.hat.m2+p2.hat.m2))[p1.hat.m2+p2.hat.m2>1]
p1.hat.m2[YY]=KK
p2.hat.m2[YY]=1-KK
p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2



#############################################################################
#### Kaplan-Meier estimator for the censoring distribution, G(t)=pr(C>t) ####
#############################################################################

#Censoring Survival
censor.S1= NULL

#Censoring hazard
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

DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
ndim=nrow(DATA) 
t0 = c(8,10,12)
ests=estIDI_Weibull(DATA$t,DATA$indicator,DATA[,4:9], t0, TRUE)
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
  return(estIDI_Weibull(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
library(parallel)
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estIDI_Weibull","coxph","Surv","basehaz"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
boot.IDI=parLapply(clust,1:nboot, gene)

stopCluster(clust)


boot.IDI=matrix(unlist(boot.IDI),nrow=3)

gene=function(s){
  DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
  DATA.boot = DATA[-s,]
  t0 = c(8,10,12)
  return(estIDI_Weibull(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estIDI_Weibull","coxph","Surv","basehaz"))
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
# 0.024534857	0.037215686	0.045469625
# 0.042471439	0.061435107	0.074188287
# 0.026712818	0.039089345	0.047723336
# 0.046158966	0.065260385	0.077365418
# 0.033606709	0.049474129	0.059935764
# 
# 
# 0.069607492	0.072623188	0.075716967
# 0.110223762	0.113003416	0.112373506
# 0.079607186	0.083730669	0.086246633
# 0.132231687	0.137168357	0.142922621
# 0.095264055	0.09881989	0.099560316






