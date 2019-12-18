library(survival)
estNRI_Weibull = function(X,delta, Z, t0, withsd)
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
data2 = data.frame(cbind(X, event=ifelse(delta==2,1,0), z1=z1, z2=z2, z3=z3, z4=z4,z5=z5))
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
# define a function to assign the subject to the class with the highest prob

loc = function(tt)
{which(tt==max(tt))}
if(withsd){
  
  #The following calculates sample influence functions
  #Estimate P(X>t), which means survival of observed time no matter what result
  
  T.tmp1 = matrix(rep(X,each = ndim),ndim,ndim)
  
  est.S.T  = apply(ifelse(T.tmp1 > X,1,0),1,mean)
  
  Martin.ST = est.lambda.censoring/est.S.T 
  Martin.ST[is.nan(Martin.ST)|is.na(Martin.ST)]=0
  M.S.integral=matrix(0,ndim,ndim)
  for(i in 1:ndim)
    for(k in 1:ndim)
      M.S.integral[i,k]=-sum(Martin.ST[which(X<=min(X[i],X[k]))])+ifelse(delta[k]==0 & X[k]<= X[i],1,0)/est.S.T[k]
  M.S.integral[is.nan(M.S.integral)|is.na(M.S.integral)]=0
}
# evaluate the NRI at each of time points
estNRI = NULL
estsd = NULL
for (j in 1:npts){
  
  # for each subject, assign to one of the 3 categories based on the estimated risks of failing from cause 1, 2 or none
  
  pp.m1 = cbind(p1.hat.m1[,j],p2.hat.m1[,j], p3.hat.m1[,j])
  pind.m1 = apply(pp.m1,1,loc)
  
  pp.m2 = cbind(p1.hat.m2[,j],p2.hat.m2[,j], p3.hat.m2[,j])
  pind.m2 = apply(pp.m2,1,loc)
  
  
  num1 = ifelse(pind.m2==1 & pind.m1!=1 & X<=t0[j] & delta ==1,1,0)/est.S.censoring-ifelse(pind.m2!=1 & pind.m1==1 & X<=t0[j] & delta ==1,1,0)/est.S.censoring
  num1 = ifelse(num1=='NaN',0,num1)
  
  den1 = ifelse(X<=t0[j] & delta ==1,1,0)/est.S.censoring
  den1 = ifelse(den1=='NaN',0,den1)
  ratio1 = sum(num1)/sum(den1)
  
  num2 = ifelse(pind.m2==2 & pind.m1!=2 & X<=t0[j] & delta ==2,1,0)/est.S.censoring-ifelse(pind.m2!=2 & pind.m1==2 & X<=t0[j] & delta ==2,1,0)/est.S.censoring
  num2 = ifelse(num2=='NaN',0,num2)
  
  den2 = ifelse(X<=t0[j] & delta ==2,1,0)/est.S.censoring
  den2 = ifelse(den2=='NaN',0,den2)
  ratio2 = sum(num2)/sum(den2)
  
  censor.ind = max(which(sort(X) <=t0[j]))
  num3 = ifelse(pind.m2==3 & pind.m1!=3 & X>t0[j],1,0)-ifelse(pind.m2!=3 & pind.m1==3 & X>t0[j],1,0)
  num3 = ifelse(num3=='NaN',0,num3)
  
  den3 = ifelse(X>t0[j],1,0)
  den3 = ifelse(den3=='NaN',0,den3)
  ratio3 = sum(num3)/sum(den3)
  
  
  estNRI = c(estNRI, (ratio1+ratio2+ratio3)/3)
  
  if(withsd){
    Psi=matrix(0,ndim,ndim)
    
    for(j in 1:ndim)
      Psi[,j]=(num1-mean(num1)/mean(den1)*den1)*(1+M.S.integral[,j])/mean(den1)+
        (num2-mean(num2)/mean(den2)*den2)*(1+M.S.integral[,j])/mean(den2)+
        (num3-mean(num3)/mean(den3)*den3)/mean(den3)
    Psi=Psi/3
    IF=rep(0,ndim)
    
    for(i in 1:ndim)
      IF[i]=mean(Psi[i,-i]+Psi[-i,i])
    estsd=c(estsd,sqrt(mean(IF^2)))
  }
}

return(cbind(estNRI,estsd))

}

DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")

ndim=nrow(DATA) 
t0 = c(8,10,12)
ests=estNRI_Weibull(DATA$t,DATA$indicator,DATA[,4:9], t0, TRUE)
NRI = t(ests[,1])
NRI_sd = t(ests[,2])
NRI
NRI_sd/sqrt(ndim)
NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)

# > NRI
# [,1]       [,2]       [,3]
# [1,] -0.0006583278 0.04171752 0.07892697
# > NRI_sd/sqrt(ndim)
# [,1]        [,2]        [,3]
# [1,] 0.0003795176 0.007467051 0.008485864
# > NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]       [,2]       [,3]
# [1,] -0.001402169 0.02708237 0.06229498
# > NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]       [,2]       [,3]
# [1,] 8.551297e-05 0.05635267 0.09555896


# 
# > NRI
# [,1]       [,2]      [,3]
# [1,] 0.09836959 0.09972487 0.1056949
# > NRI_sd/sqrt(ndim)
# [,1]        [,2]        [,3]
# [1,] 0.008571474 0.009081042 0.009143123
# > NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]       [,2]       [,3]
# [1,] 0.08156981 0.08192636 0.08777474
# > NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]      [,2]      [,3]
# [1,] 0.1151694 0.1175234 0.1236151







