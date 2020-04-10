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
  
  ids=NULL
  ids[[1]]=1:394
  ids[[2]]=395:788
  ids[[3]]=789:1182
  ids[[4]]=1183:1577
  ids[[5]]=1578:nrow(Z0)
  p11=p21=p31=p12=p22=p32=NULL
  for(cv_num in 1:5){
    ndim = length(X0)-length(ids[[cv_num]])
    npts = length(t0)
    survival=dt[-ids[[cv_num]],]
    
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
      b25=theta[22] 
      
      A1=beta11*spline11+ beta12*spline12+ beta13*spline13+ beta14*spline14+ beta15*spline15+ beta16*spline16
      A2=beta21*spline21+ beta22*spline22+ beta23*spline23+ beta24*spline24+ beta25*spline25+ beta26*spline26
      A1diff=beta11*splinediff11+ beta12*splinediff12+ beta13*splinediff13+ beta14*splinediff14+ beta15*splinediff15+ beta16*splinediff16
      A2diff=beta21*splinediff21+ beta22*splinediff22+ beta23*splinediff23+ beta24*splinediff24+ beta25*splinediff25+ beta26*splinediff26
      E1=exp(A1+b11*Z1+b12*Z2+b13*Z3+b14*Z4+b15*Z6)
      E2=exp(A2+b21*Z1+b22*Z2+b23*Z3+b24*Z4+b25*Z6)
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
    res=maxBFGS(f, start=c(-2,-1,0,1,2,3,-2,-1,0,1,2,3,0,0,0,0,0,0),constraints=list(ineqA=A, ineqB=B))
    
    
    
    A2 <- matrix(c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0
    ), 10, 22,byrow = TRUE)
    B2 <- rep(0,10)
    res2=maxBFGS(f2, start=c(-2,-1,0,1,2,3,-2,-1,0,1,2,3,0,0,0,0,0,0,0,0,0,0),constraints=list(ineqA=A2, ineqB=B2))
    
    
    
    cif=function(t0){
      dt2=dt[ids[[cv_num]],]
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
      cvr21=rep(as.matrix(dt[ids[[cv_num]],c('Z1','Z2','Z4')])%*%pars2[13:15],length(t0))
      cvr22=rep(as.matrix(dt[ids[[cv_num]],c('Z1','Z2','Z4')])%*%pars2[16:18],length(t0))
      cvr31=rep(as.matrix(dt[ids[[cv_num]],c('Z1','Z2','Z3','Z4','Z6')])%*%pars3[13:17],length(t0))
      cvr32=rep(as.matrix(dt[ids[[cv_num]],c('Z1','Z2','Z3','Z4','Z6')])%*%pars3[18:22],length(t0))
      
      F21=matrix(exp(rep(B11,each=length(ids[[cv_num]]))+cvr21)/(1+exp(rep(B11,each=length(ids[[cv_num]]))+cvr21)+exp(rep(B12,each=length(ids[[cv_num]]))+cvr22)),ncol=length(t0))
      F22=matrix(exp(rep(B12,each=length(ids[[cv_num]]))+cvr22)/(1+exp(rep(B11,each=length(ids[[cv_num]]))+cvr21)+exp(rep(B12,each=length(ids[[cv_num]]))+cvr22)),ncol=length(t0))
      
      F31=matrix(exp(rep(B21,each=length(ids[[cv_num]]))+cvr31)/(1+exp(rep(B21,each=length(ids[[cv_num]]))+cvr31)+exp(rep(B22,each=length(ids[[cv_num]]))+cvr32)),ncol=length(t0))
      F32=matrix(exp(rep(B22,each=length(ids[[cv_num]]))+cvr32)/(1+exp(rep(B21,each=length(ids[[cv_num]]))+cvr31)+exp(rep(B22,each=length(ids[[cv_num]]))+cvr32)),ncol=length(t0))
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
    p11=rbind(p11,p1.hat.m1)
    p21=rbind(p21,p2.hat.m1)
    p31=rbind(p31,p3.hat.m1)
    p12=rbind(p12,p1.hat.m2)
    p22=rbind(p22,p2.hat.m2)
    p32=rbind(p32,p3.hat.m2)
    cat(cv_num,'\n')
  }
  return(list(p11,p21,p31,p12,p22,p32))
  
}
estNRI_Gerds = function(X,delta, Z, t0, withsd)
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
ests=estNRI_Gerds(DATA$t,DATA$indicator,scale(DATA[,4:9]), t0, TRUE)
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
  temp=try(estNRI_Gerds(DATA.boot$t,DATA.boot$indicator,scale(DATA.boot[,4:9]), t0, FALSE))
  if(class(temp)=="try-error") return(NULL)
  return(temp)
}
library(parallel)
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estNRI_Gerds","get_probs","maxBFGS"))
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
boot.est=parLapply(clust,1:nboot, gene)

stopCluster(clust)


boot.NRI=matrix(unlist(lapply(boot.est,function(v) v[,1])),nrow=3)

boot.IDI=matrix(unlist(lapply(boot.est,function(v) v[,2])),nrow=3)

gene=function(s){
  DATA=read.csv("C:\\Users\\Zander Wang\\Box\\ZhengWang\\MACS\\Legacy\\survivaldata.csv")
  DATA.boot = DATA[-s,]
  set.seed(2020)
  seeds=floor(runif(1000000,1,100000000))
  set.seed(seeds[s])
  DATA.boot=DATA.boot[sample(1:nrow(DATA.boot)),]
  t0 = c(8,10,12)
  temp=try(estNRI_Gerds(DATA.boot$t,DATA.boot$indicator,scale(DATA.boot[,4:9]), t0, FALSE))
  if(class(temp)=="try-error") return(NULL)
  return(temp)
}
no_cores <- detectCores()

# Setup cluster
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, c("estNRI_Gerds","get_probs","maxBFGS"))
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

write.csv(res,"E:Gerds_macs.csv")

