library(survival)
library(cmprsk)

estNRI_Fine = function(X,delta, Z, t0, withsd)
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
ests=estNRI_Fine(DATA$t,DATA$indicator,DATA[,4:9], t0, TRUE)
NRI = t(ests[,1])
NRI_sd = t(ests[,2])
NRI
NRI_sd/sqrt(ndim)
NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)

# 
# > NRI
# [,1]       [,2]       [,3]
# [1,]    0 0.02005908 0.07266381
# > NRI_sd/sqrt(ndim)
# [,1]        [,2]        [,3]
# [1,]    0 0.004619255 0.007309066
# > NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]       [,2]       [,3]
# [1,]    0 0.01100551 0.05833831
# > NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]       [,2]       [,3]
# [1,]    0 0.02911265 0.08698932
# 
# 
# > NRI
# [,1]       [,2]       [,3]
# [1,] 0.09400682 0.09739952 0.09805999
# > NRI_sd/sqrt(ndim)
# [,1]        [,2]        [,3]
# [1,] 0.00816582 0.008032636 0.008542305
# > NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]       [,2]       [,3]
# [1,] 0.07800211 0.08165584 0.08131738
# > NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)
# [,1]      [,2]      [,3]
# [1,] 0.1100115 0.1131432 0.1148026



