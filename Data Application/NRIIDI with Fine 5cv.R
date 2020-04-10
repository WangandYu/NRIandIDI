library(survival)
library(cmprsk)

get_probs = function(X0,delta0, Z0, t0){
  ids=NULL
  ids[[1]]=1:394
  ids[[2]]=395:788
  ids[[3]]=789:1182
  ids[[4]]=1183:1577
  ids[[5]]=1578:nrow(Z0)
  p11=p21=p31=p12=p22=p32=NULL
  for(i in 1:5){
    ndim = length(X0)-length(ids[[i]])
    npts = length(t0)
    Z=Z0[-ids[[i]],]
    Y=Z0[ids[[i]],]
    X=X0[-ids[[i]]]
    delta=delta0[-ids[[i]]]
    
    # calcualte the cumulative baseline hazard function for cause 1 #
    cov = cbind(Z[,1],Z[,2],Z[,4])
    crr1 = crr(X,delta,cov)
    pred1 = predict(crr1,cbind(Y[,1],Y[,2],Y[,4]))
    # summary(crr1)
    
    # indecies for the predified t0 #
    t0.index= c(max(which(pred1[,1]<t0[1])),max(which(pred1[,1]<t0[2])),max(which(pred1[,1]<t0[3])))
    
    # estimate the overall survival function #
    p1.hat.m1 = t(pred1[t0.index,-1])
    
    crr2 = crr(X,delta,cov,failcode=2)
    pred2 = predict(crr2,cbind(Y[,1],Y[,2],Y[,4]))
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
    
    p11=rbind(p11,p1.hat.m1)
    p21=rbind(p21,p2.hat.m1)
    p31=rbind(p31,p3.hat.m1)
    
    
    # calcualte the cumulative baseline hazard function for cause 1 #
    cov = cbind(Z[,1],Z[,2],Z[,3],Z[,4],Z[,6])
    crr1 = crr(X,delta,cov)
    pred1 = predict(crr1,cbind(Y[,1],Y[,2],Y[,3],Y[,4],Y[,6]))
    # summary(crr1)
    
    # indecies for the predified t0 #
    t0.index= c(max(which(pred1[,1]<t0[1])),max(which(pred1[,1]<t0[2])),max(which(pred1[,1]<t0[3])))
    
    # estimate the overall survival function #
    p1.hat.m2 = t(pred1[t0.index,-1])
    
    cov = cbind(Z[,1],Z[,2],Z[,3],Z[,4],Z[,6])
    crr2 = crr(X,delta,cov,failcode=2)
    pred2 = predict(crr2,cbind(Y[,1],Y[,2],Y[,3],Y[,4],Y[,6]))
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
    
    
    p12=rbind(p12,p1.hat.m2)
    p22=rbind(p22,p2.hat.m2)
    p32=rbind(p32,p3.hat.m2)
  }
  return(list(p11,p21,p31,p12,p22,p32))
  
}
estNRI_Fine = function(X,delta, Z, t0, withsd)
{ndim = length(X)
npts = length(t0)

res=get_probs(X,delta,Z,t0)
save(res,file = 'E:Legacy_Fine_probs.rda')
p1.hat.m1=res[[1]]
p2.hat.m1=res[[2]]
p3.hat.m1=res[[3]]
p1.hat.m2=res[[4]]
p2.hat.m2=res[[5]]
p3.hat.m2=res[[6]]


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

return(cbind(estNRI,c(estIDI.t1, estIDI.t2, estIDI.t3),estsd))

}

DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
set.seed(2020)
DATA=DATA[sample(1:nrow(DATA)),]
ndim=nrow(DATA) 
t0 = c(8,10,12)
ests=estNRI_Fine(DATA$t,DATA$indicator,DATA[,4:9], t0, TRUE)
NRI = t(ests[,1])
IDI = t(ests[,2])
# NRI_sd = t(ests[,3])
# NRI
# NRI_sd/sqrt(ndim)
# NRI-qnorm(0.975)*NRI_sd/sqrt(ndim)
# NRI+qnorm(0.975)*NRI_sd/sqrt(ndim)


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
  return(estNRI_Fine(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estNRI_Fine","crr","get_probs"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
boot.est=parLapply(clust,1:nboot, gene)

stopCluster(clust)


boot.NRI=matrix(unlist(lapply(boot.est,function(v) v[,1])),nrow=3)

boot.IDI=matrix(unlist(lapply(boot.est,function(v) v[,2])),nrow=3)

gene=function(s){
  DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
  DATA.boot = DATA[-s,]
  t0 = c(8,10,12)
  set.seed(2020)
  seeds=floor(runif(1000000,1,100000000))
  set.seed(seeds[s])
  DATA.boot=DATA.boot[sample(1:nrow(DATA.boot)),]
  return(estNRI_Fine(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estNRI_Fine","crr","get_probs"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
jack.est=parLapply(clust,1:nrow(DATA), gene)

stopCluster(clust)


jack.NRI=matrix(unlist(lapply(jack.est,function(v) v[,1])),nrow=3)

jack.IDI=matrix(unlist(lapply(jack.est,function(v) v[,2])),nrow=3)

z0.hat=qnorm(colMeans(t(boot.NRI)<matrix(rep(t(NRI),each=1000),ncol=3)))
theta0.hat=apply(jack.NRI,1,mean)
temp=-t(jack.NRI)+matrix(rep(theta0.hat,each=ncol(jack.NRI)),ncol=3)
a.hat=colSums(temp^3)/6/colSums(temp^2)^(1.5)
alpha1=pnorm(z0.hat+(z0.hat+qnorm(0.025))/(1-a.hat*(z0.hat+qnorm(0.025))))
alpha2=pnorm(z0.hat+(z0.hat+qnorm(0.975))/(1-a.hat*(z0.hat+qnorm(0.975))))
res=rbind(c(quantile(boot.NRI[1,],alpha1[1]),quantile(boot.NRI[2,],alpha1[2]),quantile(boot.NRI[3,],alpha1[3])),
          c(quantile(boot.NRI[1,],alpha2[1]),quantile(boot.NRI[2,],alpha2[2]),quantile(boot.NRI[3,],alpha2[3])),
          c(quantile(boot.NRI[1,],0.025),quantile(boot.NRI[2,],0.025),quantile(boot.NRI[3,],0.025)),
          c(quantile(boot.NRI[1,],0.975),quantile(boot.NRI[2,],0.975),quantile(boot.NRI[3,],0.975)),
          NRI)


z0.hat=qnorm(colMeans(t(boot.IDI)<matrix(rep(t(IDI),each=1000),ncol=3)))
theta0.hat=apply(jack.IDI,1,mean)
temp=-t(jack.IDI)+matrix(rep(theta0.hat,each=ncol(jack.IDI)),ncol=3)
a.hat=colSums(temp^3)/6/colSums(temp^2)^(1.5)
alpha1=pnorm(z0.hat+(z0.hat+qnorm(0.025))/(1-a.hat*(z0.hat+qnorm(0.025))))
alpha2=pnorm(z0.hat+(z0.hat+qnorm(0.975))/(1-a.hat*(z0.hat+qnorm(0.975))))
res=rbind(res,
          c(quantile(boot.IDI[1,],alpha1[1]),quantile(boot.IDI[2,],alpha1[2]),quantile(boot.IDI[3,],alpha1[3])),
          c(quantile(boot.IDI[1,],alpha2[1]),quantile(boot.IDI[2,],alpha2[2]),quantile(boot.IDI[3,],alpha2[3])),
          c(quantile(boot.IDI[1,],0.025),quantile(boot.IDI[2,],0.025),quantile(boot.IDI[3,],0.025)),
          c(quantile(boot.IDI[1,],0.975),quantile(boot.IDI[2,],0.975),quantile(boot.IDI[3,],0.975)),
          IDI)

write.csv(res,"E:Fine_macs.csv")

