




est_Gerds = function(name){
  library(maxLik)
  ndim = 1000
  
  z1 = rnorm(ndim)					# continuous 
  z1 = pmax(pmin(z1,3.5),-3.5)
  z2 = ifelse(runif(ndim,0,1)<=0.7,1,0)		# 0,1 variable
  z3 = rnorm(ndim)					# continuous with large coefficient
  z3 = pmax(pmin(z3,3.5),-3.5)
  beta=c(1,1,3)
  t=(log((1/runif(ndim)-1)/2)-cbind(z1,z2,z3)%*%beta)
  t=pmax(-100,pmin(t,100))
  #cen=runif(ndim,-11,21)
  cause=ifelse(runif(ndim)<0.5,1,2)
  #cause=ifelse(cen<t,0,cause)
  #t=pmin(cen,t)
  dt=data.frame(time=t,event=cause,Z1=z1,Z2=z2,Z3=z3)
  t0=c(-2,-1,0)
  
  
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
  survival=dt
  
  splines=function(data,ndim.splines=6){
    names(data)=c('t','cause','z1','z2','z3')
    q1=q2=quantile(data$t,(0:3)/3)
    percs1=c(rep(min(q1),3),q1,rep(max(q1),3))
    percs2=c(rep(min(q2),3),q2,rep(max(q2),3))
    res=NULL
    for(i in 1:nrow(data)){
      A1=A2=A1.diff=A2.diff=NULL
      for(ii in 1:(2+length(q1))){
        A1=c(A1,b.spline(ii,3,percs1,data$t[i]))
        A2=c(A2,b.spline(ii,3,percs2,data$t[i]))
        A1.diff=c(A1.diff,b.spline.diff(ii,3,percs1,data$t[i]))
        A2.diff=c(A2.diff,b.spline.diff(ii,3,percs2,data$t[i]))
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
    E1=exp(A1+b11*Z1+b12*Z2+b13*Z3)
    E2=exp(A2+b21*Z1+b22*Z2+b23*Z3)
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
    b21=theta[15] 
    b22=theta[16] 
    
    A1=beta11*spline11+ beta12*spline12+ beta13*spline13+ beta14*spline14+ beta15*spline15+ beta16*spline16
    A2=beta21*spline21+ beta22*spline22+ beta23*spline23+ beta24*spline24+ beta25*spline25+ beta26*spline26
    A1diff=beta11*splinediff11+ beta12*splinediff12+ beta13*splinediff13+ beta14*splinediff14+ beta15*splinediff15+ beta16*splinediff16
    A2diff=beta21*splinediff21+ beta22*splinediff22+ beta23*splinediff23+ beta24*splinediff24+ beta25*splinediff25+ beta26*splinediff26
    E1=exp(A1+b11*Z1+b12*Z2)
    E2=exp(A2+b21*Z1+b22*Z2)
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
  res=try(maxBFGS(f, start=c(-13,-8,-6,1,3,4,-13,-8,-6,1,3,4,1,1,3,1,1,3),constraints=list(ineqA=A, ineqB=B)),silent = TRUE)
  if(class(res)=="try-error") return(NULL)
  #print(summary(res))
  
  
  A2 <- matrix(c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0
  ), 10, 16,byrow = TRUE)
  B2 <- rep(0,10)
  res2=try(maxBFGS(f2, start=c(-8,-5,-2,0,1,2,-8,-5,-2,0,1,2,0.5,0.5,0.5,0.5),constraints=list(ineqA=A2, ineqB=B2)),silent = TRUE)
  #print(summary(res2))
  if(class(res2)=="try-error") return(NULL)
  
  
  
  cif=function(t0){
    dt=survival
    pars2=summary(res2)[["estimate"]][,1]
    pars3=summary(res)[["estimate"]][,1]
    B11=B12=B21=B22=NULL
    q1=q2=quantile(dt$time,(0:3)/3)
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
    cvr21=rep(as.matrix(dt[,c('Z1','Z2')])%*%pars2[13:14],length(t0))
    cvr22=rep(as.matrix(dt[,c('Z1','Z2')])%*%pars2[15:16],length(t0))
    cvr31=rep(as.matrix(dt[,c('Z1','Z2','Z3')])%*%pars3[13:15],length(t0))
    cvr32=rep(as.matrix(dt[,c('Z1','Z2','Z3')])%*%pars3[16:18],length(t0))
    
    F21=matrix(exp(rep(B11,each=nrow(dt))+cvr21)/(1+exp(rep(B11,each=nrow(dt))+cvr21)+exp(rep(B12,each=nrow(dt))+cvr22)),ncol=length(t0))
    F22=matrix(exp(rep(B12,each=nrow(dt))+cvr22)/(1+exp(rep(B11,each=nrow(dt))+cvr21)+exp(rep(B12,each=nrow(dt))+cvr22)),ncol=length(t0))
    
    F31=matrix(exp(rep(B21,each=nrow(dt))+cvr31)/(1+exp(rep(B21,each=nrow(dt))+cvr31)+exp(rep(B22,each=nrow(dt))+cvr32)),ncol=length(t0))
    F32=matrix(exp(rep(B22,each=nrow(dt))+cvr32)/(1+exp(rep(B21,each=nrow(dt))+cvr31)+exp(rep(B22,each=nrow(dt))+cvr32)),ncol=length(t0))
    return(list(F21,F22,1-F21-F22,F31,F32,1-F31-F32))
  }
  
  
  npts=length(t0)
  
  cifs=cif(t0)
  p1.hat.m1=cifs[[1]]
  p2.hat.m1=cifs[[2]]
  p3.hat.m1=cifs[[3]]
  p1.hat.m2=cifs[[4]]
  p2.hat.m2=cifs[[5]]
  p3.hat.m2=cifs[[6]]
  
  TT=dt$time
  epsilon=dt$event
  ndim=nrow(dt)
  
  n1=n2=n3=NULL
  for(t00 in t0){
    n1 = c(n1,sum(ifelse(TT<=t00 & epsilon==1,1,0)))
    n2 = c(n2,sum(ifelse(TT<=t00 & epsilon==2,1,0)))
    n3 = c(n3,ndim - sum(ifelse(TT<=t00 & epsilon==1,1,0)) - sum(ifelse(TT<=t00 & epsilon==2,1,0)))
  }
  
  sum1=sum2=sum3=NULL
  for(i in 1:3){
    sum1 = c(sum1,sum((p1.hat.m2[,i]-mean(p1.hat.m2[,i]))^2 - (p1.hat.m1[,i]-mean(p1.hat.m1[,i]))^2))
    sum2 = c(sum2,sum((p2.hat.m2[,i]-mean(p2.hat.m2[,i]))^2 - (p2.hat.m1[,i]-mean(p2.hat.m1[,i]))^2))
    sum3 = c(sum3,sum((p3.hat.m2[,i]-mean(p3.hat.m2[,i]))^2 - (p3.hat.m1[,i]-mean(p3.hat.m1[,i]))^2))
  }
  
  pind.m1=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m1[,i],p2.hat.m1[,i], p3.hat.m1[,i])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m1 = cbind(pind.m1,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  
  pind.m2=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m2[,i],p2.hat.m2[,i], p3.hat.m2[,i])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m2 = cbind(pind.m2,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  
  sum10=sum20=sum30=NULL
  for(i in 1:3){
    t00=t0[i]
    sum10 = c(sum10,sum(ifelse(pind.m2[,i]==1 & pind.m1[,i]!=1 & TT<=t00 & epsilon ==1, 1, 0)) - sum(ifelse(pind.m2[,i]!=1 & pind.m1[,i]==1 & TT<=t00 & epsilon ==1, 1, 0)))
    sum20 = c(sum20,sum(ifelse(pind.m2[,i]==2 & pind.m1[,i]!=2 & TT<=t00 & epsilon ==2, 1, 0)) - sum(ifelse(pind.m2[,i]!=2 & pind.m1[,i]==2 & TT<=t00 & epsilon ==2, 1, 0)))
    sum30 = c(sum30,sum(ifelse(pind.m2[,i]==3 & pind.m1[,i]!=3 & TT>t00, 1, 0)) - sum(ifelse(pind.m2[,i]!=3 & pind.m1[,i]==3 & TT>t00, 1, 0)))
  }
  
  return(rbind(n1,n2,n3,sum1,sum2,sum3,sum10,sum20,sum30))
  
}
library(parallel)
# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
cl <- makeCluster(no_cores)
GerdsTrueNRIIDI=parLapply(cl,1:1000,est_Gerds)
save(GerdsTrueNRIIDI,file = 'E:\\Research\\NRIIDI revisions\\Gerds\\GerdsTrueNRIIDI.rda')

load("E:\\Research\\NRIIDI revisions\\Gerds\\GerdsTrueNRIIDI.rda")
NRI=IDI=NULL
for(j in 1:3){
  ns=sums=sums0=NULL
  for(i in 1:1000){
    if(is.null(GerdsTrueNRIIDI[[i]])) next
    ns=rbind(ns,GerdsTrueNRIIDI[[i]][1:3,j])
    sums=rbind(sums,GerdsTrueNRIIDI[[i]][4:6,j])
    sums0=rbind(sums0,GerdsTrueNRIIDI[[i]][7:9,j])
  }
  ns=apply(ns,2,sum)
  NRI=c(NRI,sum(apply(sums0,2,sum)/ns)/3)
  IDI=c(IDI,sum(apply(sums,2,sum)/ns/(1-ns/sum(ns)))/3)
}
NRI
IDI
