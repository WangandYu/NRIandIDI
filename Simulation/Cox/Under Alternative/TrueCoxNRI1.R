
library(survival)
# sample size
ndim = 1000
nsim = 1000

### true parameters ###
beta0 = 2.5
beta1 = 0.05
beta2 = -0.05
beta3 = 0.15
sigma = 0.2
beta = c(beta1,beta2,beta3)			# regression coefficients


T0 = c(11,12,13)

for(m in 1:nsim)
{
  npts=3
  ### generate the covariates ###
  z1 = rnorm(ndim)					# continuous 
  z1 = pmax(pmin(z1,3.5),-3.5)
  z2 = ifelse(runif(ndim,0,1)<=0.7,1,0)		# 0,1 variable
  z3 = rnorm(ndim)					# continuous with large coefficient
  z3 = pmax(pmin(z3,3.5),-3.5)
  
  ## generate event time using Weibull model ##
  z = cbind(z1,z2,z3)					# covariates
  U = runif(ndim, 0, 1)				# uniform(0,1)
  w = log(-log(1-U))					# inverse cdf of extreme value dist
  TT = exp(beta0 + z%*%beta + sigma*w)		# Weibull model
 
  
  epsilon = ifelse(runif(ndim)<0.5,1,2)
  
  
  for(t0 in T0){
    n1 = sum(ifelse(TT<=t0 & epsilon==1,1,0))
    n2 = sum(ifelse(TT<=t0 & epsilon==2,1,0))
    n3 = ndim - n1 - n2
    
    write.table(t(c(n1,n2,n3)),paste("NRI_n_",t0,"_",name,".txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE,append = TRUE)
  }
  #############################################################################
  ### Fit the cox proportional hazard regression model without z3 (model 1) ###
  #############################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  data1 = data.frame(cbind(TT, event=ifelse(epsilon==1,1,0), z1=z1, z2=z2))
  cox1 = coxph(Surv(TT, event) ~ z1 + z2, data1)	
  # summary(cox1)
  cox.basehaz1 = basehaz(cox1, centered=FALSE)
  cox.basehaz1 = cox.basehaz1[cox.basehaz1$time<=t0,]	# only keep the baseline hazard until time t0
  
  # Then calculate the comulative baseline hazard fuction for cause 2 #
  data2 = data.frame(cbind(TT, event=ifelse(epsilon==2,1,0), z1=z1, z2=z2))
  cox2 = coxph(Surv(TT, event) ~ z1 + z2, data2)	
  # summary(cox2)
  cox.basehaz2 = basehaz(cox2, centered=FALSE)
  cox.basehaz2 = cox.basehaz2[cox.basehaz2$time<=t0,]	# only keep the baseline hazard until time t0[2]
  
  # indecies for the predified t0 #
  
  t0.index1 = max(which(cox.basehaz1$time<=T0[1]))
  t0.index2 = max(which(cox.basehaz1$time<=T0[2]))
  t0.index3 = max(which(cox.basehaz1$time<=T0[3]))
  t0.index = c(t0.index1, t0.index2, t0.index3)
  
  dim.uniquetime = length(which(cox.basehaz1$time<=T0[npts]))	
  
  
  # estimate the cumulative cause specific hazard for cause 1 #
  cox.haz1 = matrix(rep(cox.basehaz1$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
  cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z1,z2)%*%cox1$coef),dim.uniquetime),ndim,dim.uniquetime)
  # calcualte the cumulative cause specific hazard for cause 2 #
  cox.haz2 = matrix(rep(cox.basehaz2$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
  cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z1,z2)%*%cox2$coef),dim.uniquetime),ndim,dim.uniquetime)
  
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
  p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
  
  pind.m1=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m1[,i],p2.hat.m1[,i], p3.hat.m1[,i])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m1 = cbind(pind.m1,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  
  
  #######################################################################################
  ### Fit the cox proportional hazard regression model using all covariates (model 2) ###
  #######################################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  data1 = data.frame(cbind(TT, event=ifelse(epsilon==1,1,0), z1=z1, z2=z2, z3=z3))
  cox1 = coxph(Surv(TT, event) ~ z1 + z2 + z3, data1)	
  #summary(cox1)
  cox.basehaz1 = basehaz(cox1, centered=FALSE)
  cox.basehaz1 = cox.basehaz1[cox.basehaz1$time<=t0,]	# only keep the baseline hazard until time t0
  
  # Then calculate the comulative baseline hazard fuction for cause 2 #
  data2 = data.frame(cbind(TT, event=ifelse(epsilon==2,1,0), z1=z1, z2=z2, z3=z3))
  cox2 = coxph(Surv(TT, event) ~ z1 + z2 + z3, data2)	
  #summary(cox2)
  cox.basehaz2 = basehaz(cox2, centered=FALSE)
  cox.basehaz2 = cox.basehaz2[cox.basehaz2$time<=t0,]	# only keep the baseline hazard until time t0
  
  # index for cause 1 and 2 events at the predified t0 #
  t0.index1 = max(which(cox.basehaz1$time<=T0[1]))
  t0.index2 = max(which(cox.basehaz1$time<=T0[2]))
  t0.index3 = max(which(cox.basehaz1$time<=T0[3]))
  t0.index = c(t0.index1, t0.index2, t0.index3)
  
  dim.uniquetime = length(which(cox.basehaz1$time<=T0[npts]))	
  
  # estimate the cumulative cause specific hazard for cause 1 #
  cox.haz1 = matrix(rep(cox.basehaz1$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
  cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z1,z2,z3)%*%cox1$coef),dim.uniquetime),ndim,dim.uniquetime)
  
  # calcualte the cumulative cause specific hazard for cause 2 #
  cox.haz2 = matrix(rep(cox.basehaz2$hazard,each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
  cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z1,z2,z3)%*%cox2$coef),dim.uniquetime),ndim,dim.uniquetime)
  
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
  p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
  
  pind.m2=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m2[,i],p2.hat.m2[,i], p3.hat.m2[,i])
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
    
    write.table(t(c(sum1,sum2,sum3)),paste("NRI_sum_",t0,"_",name,".txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE,append = TRUE)
  }
}
dir2=''
ns = read.table(paste(dir2,"NRI_n_11.txt",sep=""),header = FALSE)
dim(ns)
ns = apply(ns,2,sum)

sums = read.table(paste(dir2,"NRI_sum_11.txt",sep=""),header = FALSE)
dim(sums)
sums = apply(sums,2,sum)

sum(sums/ns)/3 
ns = read.table(paste(dir2,"NRI_n_12.txt",sep=""),header = FALSE)
dim(ns)
ns = apply(ns,2,sum)

sums = read.table(paste(dir2,"NRI_sum_12.txt",sep=""),header = FALSE)
dim(sums)
sums = apply(sums,2,sum)

sum(sums/ns)/3 
ns = read.table(paste(dir2,"NRI_n_13.txt",sep=""),header = FALSE)
dim(ns)
ns = apply(ns,2,sum)

sums = read.table(paste(dir2,"NRI_sum_13.txt",sep=""),header = FALSE)
dim(sums)
sums = apply(sums,2,sum)

sum(sums/ns)/3 