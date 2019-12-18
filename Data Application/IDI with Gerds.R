library(maxLik)
get_probs = function(X0,delta0, Z0, t0){
  
  omega=function(i,j,percs,t){
    if(percs[i]==percs[i+j]) return(0)
    return((t-percs[i])/(percs[i+j]-percs[i]))
  }
  omega.diff=function(i,j,percs){
    if(percs[i]==percs[i+j]) return(0)
    return(1/(percs[i+j]-percs[i]))
  }
  b.spline=function(i,j,percs,t){
    if(t<percs[i]|t>percs[i+j+1]) return(0)
    if(j==0)
      return(1)
    return(omega(i,j,percs,t)*b.spline(i,j-1,percs,t)+(1-omega(i+1,j,percs,t))*b.spline(i+1,j-1,percs,t))
  }
  b.spline.diff=function(i,j,percs,t){
    if(t<percs[i]|t>percs[i+j+1]) return(0)
    if(j==0)
      return(0)
    return(omega(i,j,percs,t)*b.spline.diff(i,j-1,percs,t)+(1-omega(i+1,j,percs,t))*b.spline.diff(i+1,j-1,percs,t)+
             omega.diff(i,j,percs)*b.spline(i,j-1,percs,t)-omega.diff(i+1,j,percs)*b.spline(i+1,j-1,percs,t))
  }
  dt=data.frame(time=X0,event=delta0,Z1=Z0[,1],Z2=Z0[,2],Z3=Z0[,3],Z4=Z0[,4],Z5=Z0[,5],Z6=Z0[,6])
  
  ids=1:length(X0)
  ndim = length(X0)-length(ids)
  npts = length(t0)
  survival=dt[-ids,]
  
  splines=function(data,ndim.splines=6){
    names(data)=c('t','cause','Z1','Z2','Z3','Z4','Z5','Z6')
    q1=q2=quantile(data[['t']],(0:3)/3)
    percs1=c(rep(min(q1),3),q1,rep(max(q1),3))
    percs2=c(rep(min(q2),3),q2,rep(max(q2),3))
    res=NULL
    for(i in 1:nrow(data)){
      A1=A2=A1.diff=A2.diff=NULL
      for(ii in 1:(2+length(q1))){
        A1=c(A1,b.spline(ii,3,percs1,data[['t']][i]))
        A2=c(A2,b.spline(ii,3,percs2,data[['t']][i]))
        A1.diff=c(A1.diff,b.spline.diff(ii,3,percs1,data[['t']][i]))
        A2.diff=c(A2.diff,b.spline.diff(ii,3,percs2,data[['t']][i]))
      }
      res=rbind(res,c(A1,A2,A1.diff,A2.diff))
    }
    return(res)
  }
  sp=as.data.frame(splines(survival))
  names(sp)=c("spline11", "spline12", "spline13", "spline14", "spline15", "spline16",  
              "spline21", "spline22", "spline23", "spline24", "spline25", "spline26",  
              "splinediff11", "splinediff12", "splinediff13", "splinediff14", "splinediff15", "splinediff16", 
              "splinediff21", "splinediff22", "splinediff23", "splinediff24", "splinediff25", "splinediff26")
  survival=cbind(survival,sp)
  
  
  f=function(theta){
    if("survival"%in%search()) detach(survival)
    attach(survival)
    beta11=theta[1] 
    beta12=theta[2] 
    beta13=theta[3] 
    beta14=theta[4] 
    beta15=theta[5] 
    beta16=theta[6]  
    beta21=theta[7] 
    beta22=theta[8] 
    beta23=theta[9] 
    beta24=theta[10] 
    beta25=theta[11] 
    beta26=theta[12] 
    
    b11=theta[13] 
    b12=theta[14] 
    b13=theta[15] 
    b21=theta[16] 
    b22=theta[17] 
    b23=theta[18] 
    
    A1=beta11*spline11+ beta12*spline12+ beta13*spline13+ beta14*spline14+ beta15*spline15+ beta16*spline16
    A2=beta21*spline21+ beta22*spline22+ beta23*spline23+ beta24*spline24+ beta25*spline25+ beta26*spline26
    A1diff=beta11*splinediff11+ beta12*splinediff12+ beta13*splinediff13+ beta14*splinediff14+ beta15*splinediff15+ beta16*splinediff16
    A2diff=beta21*splinediff21+ beta22*splinediff22+ beta23*splinediff23+ beta24*splinediff24+ beta25*splinediff25+ beta26*splinediff26
    E1=exp(A1+b11*Z1+b12*Z2+b13*Z4)
    E2=exp(A2+b21*Z1+b22*Z2+b23*Z4)
    ll=0
    for(k in 1:length(event)){
      if( event[k]==1 ) 
        lh = (A1diff[k]*E1[k]+(A1diff[k]-A2diff[k])*E1[k]*E2[k])/(1+E1[k]+E2[k])^2
      else if(event[k]==2 )
        lh = (A2diff[k]*E2[k]+(A2diff[k]-A1diff[k])*E1[k]*E2[k])/(1+E1[k]+E2[k])^2
      else if(event[k]==0 )
        lh = 1/(1+E1[k]+E2[k]) 
      
      if (lh > 0) 
        ll = ll+log(lh)
      else 
        ll = ll-1e20
    }
    detach(survival)
    return(ll)
  }
  
  
  f2=function(theta){
    if("survival"%in%search()) detach(survival)
    attach(survival)
    beta11=theta[1] 
    beta12=theta[2] 
    beta13=theta[3] 
    beta14=theta[4] 
    beta15=theta[5] 
    beta16=theta[6]  
    beta21=theta[7] 
    beta22=theta[8] 
    beta23=theta[9] 
    beta24=theta[10] 
    beta25=theta[11] 
    beta26=theta[12] 
    
    b11=theta[13] 
    b12=theta[14] 
    b13=theta[15] 
    b14=theta[16] 
    b15=theta[17] 
    b21=theta[18] 
    b22=theta[19] 
    b23=theta[20] 
    b24=theta[21] 
    
    A1=beta11*spline11+ beta12*spline12+ beta13*spline13+ beta14*spline14+ beta15*spline15+ beta16*spline16
    A2=beta21*spline21+ beta22*spline22+ beta23*spline23+ beta24*spline24+ beta25*spline25+ beta26*spline26
    A1diff=beta11*splinediff11+ beta12*splinediff12+ beta13*splinediff13+ beta14*splinediff14+ beta15*splinediff15+ beta16*splinediff16
    A2diff=beta21*splinediff21+ beta22*splinediff22+ beta23*splinediff23+ beta24*splinediff24+ beta25*splinediff25+ beta26*splinediff26
    E1=exp(A1+b11*Z1+b12*Z2+b13*Z3+b14*Z4+b15*Z6)
    E2=exp(A2+b21*Z1+b22*Z2+b23*Z3+b24*Z4)
    ll=0
    for(k in 1:length(event)){
      if( event[k]==1 ) 
        lh = (A1diff[k]*E1[k]+(A1diff[k]-A2diff[k])*E1[k]*E2[k])/(1+E1[k]+E2[k])^2
      else if(event[k]==2 )
        lh = (A2diff[k]*E2[k]+(A2diff[k]-A1diff[k])*E1[k]*E2[k])/(1+E1[k]+E2[k])^2
      else
        lh = 1/(1+E1[k]+E2[k]) 
      
      if (lh > 0) 
        ll = ll+log(lh)
      else 
        ll = ll-1e20
    }
    detach(survival)
    return(ll)
  }
  
  A <- matrix(c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0
  ), 10, 18,byrow = TRUE)
  B <- rep(0,10)
  res=maxBFGS(f, start=c(-2,-1,0,1,2,3,-2,-1,0,1,2,3,0,0,0,0,0,0),constraints=list(ineqA=A, ineqB=B),control=list(printLevel=1))
  
  
  
  A2 <- matrix(c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0
  ), 10, 21,byrow = TRUE)
  B2 <- rep(0,10)
  res2=maxBFGS(f2, start=c(-2,-1,0,1,2,3,-2,-1,0,1,2,3,0,0,0,0,0,0,0,0,0),constraints=list(ineqA=A2, ineqB=B2),control=list(printLevel=1))
  
  
  
  cif=function(t0){
    dt2=dt[ids,]
    pars2=summary(res)[["estimate"]][,1]
    pars3=summary(res2)[["estimate"]][,1]
    B11=B12=B21=B22=NULL
    q1=q2=quantile(dt2[['time']],(0:3)/3)
    percs1=c(rep(min(q1),3),q1,rep(max(q1),3))
    percs2=c(rep(min(q2),3),q2,rep(max(q2),3))
    for(i in t0){
      A1=A2=NULL
      for(ii in 1:6){
        A1=c(A1,b.spline(ii,3,percs1,i))
        A2=c(A2,b.spline(ii,3,percs2,i))
      }
      B11=c(B11,sum(A1*pars2[1:6]))
      B12=c(B12,sum(A2*pars2[7:12]))
      B21=c(B21,sum(A1*pars3[1:6]))
      B22=c(B22,sum(A2*pars3[7:12]))
    }
    cvr21=rep(as.matrix(dt[ids,c('Z1','Z2','Z4')])%*%pars2[13:15],length(t0))
    cvr22=rep(as.matrix(dt[ids,c('Z1','Z2','Z4')])%*%pars2[16:18],length(t0))
    cvr31=rep(as.matrix(dt[ids,c('Z1','Z2','Z3','Z4','Z6')])%*%pars3[13:17],length(t0))
    cvr32=rep(as.matrix(dt[ids,c('Z1','Z2','Z3','Z4')])%*%pars3[18:21],length(t0))
    
    F21=matrix(exp(rep(B11,each=length(ids))+cvr21)/(1+exp(rep(B11,each=length(ids))+cvr21)+exp(rep(B12,each=length(ids))+cvr22)),ncol=length(t0))
    F22=matrix(exp(rep(B12,each=length(ids))+cvr22)/(1+exp(rep(B11,each=length(ids))+cvr21)+exp(rep(B12,each=length(ids))+cvr22)),ncol=length(t0))
    
    F31=matrix(exp(rep(B21,each=length(ids))+cvr31)/(1+exp(rep(B21,each=length(ids))+cvr31)+exp(rep(B22,each=length(ids))+cvr32)),ncol=length(t0))
    F32=matrix(exp(rep(B22,each=length(ids))+cvr32)/(1+exp(rep(B21,each=length(ids))+cvr31)+exp(rep(B22,each=length(ids))+cvr32)),ncol=length(t0))
    return(list(F21,F22,1-F21-F22,F31,F32,1-F31-F32))
  }
  
  
  npts=length(t0)
  
  cifs=cif(t0)
  # p1.hat.m1=cifs[[1]]
  # p2.hat.m1=cifs[[2]]
  # p3.hat.m1=cifs[[3]]
  # p1.hat.m2=cifs[[4]]
  # p2.hat.m2=cifs[[5]]
  # p3.hat.m2=cifs[[6]]
  return(cifs)
  
}
estIDI_Gerds = function(X,delta, Z, t0, withsd)
{
  ndim=length(X)
  npts = length(t0)
  res=get_probs(X,delta,Z,t0)
  save(res,file = 'E:Legacy_Gerds_probs.rda')
  p1.hat.m1=res[[1]]
  p2.hat.m1=res[[2]]
  p3.hat.m1=res[[3]]
  p1.hat.m2=res[[4]]
  p2.hat.m2=res[[5]]
  p3.hat.m2=res[[6]]
  

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
ests=estIDI_Gerds(DATA$t,DATA$indicator,DATA[,4:9], t0, TRUE)
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
  return(estIDI_Gerds(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
library(parallel)
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estIDI_Gerds","get_probs","maxBFGS"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
boot.IDI=parLapply(clust,1:nboot, gene)

stopCluster(clust)


boot.IDI=matrix(unlist(boot.IDI),nrow=3)

gene=function(s){
  DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
  DATA.boot = DATA[-s,]
  t0 = c(8,10,12)
  return(estIDI_Gerds(DATA.boot$t,DATA.boot$indicator,DATA.boot[,4:9], t0, FALSE))
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estIDI_Gerds","get_probs","maxBFGS"))
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






