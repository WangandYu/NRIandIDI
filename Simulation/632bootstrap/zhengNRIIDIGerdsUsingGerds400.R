#!/bin/bash
#SBATCH --job-name=NRIGG400
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=25
#SBATCH --cluster=mpi
#SBATCH --partition=opa
#SBATCH --mail-user=zhw61@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --time=6-00:00:00

module purge
module load gcc/8.2.0
module load r/3.6.0

#This is a comment

export OMPI_MCA_rmaps_base_oversubscribe=1

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl="self,tcp"

#cp $SLURM_SUBMIT_DIR/* $SLURM_SCRATCH
#cd $SLURM_SCRATCH

cd $SLURM_SUBMIT_DIR

#echo $SLURM_JOB_NODELIST > slurm_nodefile.txt

scontrol show hostnames | awk "{for(i=0;i<$SLURM_NTASKS_PER_NODE;i++)print}" > node_list.txt

# Telling cluster that you are using R
R --vanilla > mysnow.out <<EOF

#################################################################################
# Loading snowfall and looking for available nodes, then initializing:
#################################################################################
Sys.setenv(OMPI_MCA_btl="tcp,self")

library(snowfall)
##pbsnodefile = Sys.getenv("PBS_NODEFILE")
##pbsnodefile = Sys.getenv("SLURM_JOB_NODELIST")
##machines <- scan(pbsnodefile, what="")

machines <- scan("$SLURM_SUBMIT_DIR/node_list.txt", what="")
machines
nmach = length(machines)

sfInit(parallel=TRUE,type='MPI',cpus=nmach,socketHosts=machines)

#################################################################################
# Loading Other necessary R Libraries:
#################################################################################
sfLibrary(maxLik)

#################################################################################
#  Load External Data
#################################################################################


#################################################################################
#  Main function
#################################################################################



#################################################################################
#  Exporting the data and functions needed by the worker nodes:
#################################################################################
sfExportAll()
#################################################################################
#  Creating a wrapper function
#################################################################################


wrapper <- function(name){
  
  estNRI_Gerds = function(dt, t0, to_boot=FALSE)
  {
    ndim = nrow(dt)
    npts = length(t0)
    
    Time=dt[['time']]
    Event=dt[['event']]
    
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
    
    if(to_boot){
      idx_train=sample(1:ndim, ndim, replace=TRUE)
    } else {
      idx_train=1:ndim
    }
    X=Time[idx_train]
    delta=Event[idx_train]
    survival=dt
    
    splines=function(data,ndim.splines=6){
      names(data)=c('t','cause','z1','z2','z3')
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
    surv_temp=cbind(survival,sp)
    survival=surv_temp[idx_train,]
    
    
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
    if(class(res2)=="try-error") return(NULL)
    
    
    
    cif=function(t0){
      dt=survival_test
      pars2=summary(res2)[["estimate"]][,1]
      pars3=summary(res)[["estimate"]][,1]
      B11=B12=B21=B22=NULL
      q1=q2=quantile(dt[['time']],(0:3)/3)
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
    output_final=NULL
    
    idx_test=idx_train
    npts=length(t0)
    
    survival_test=surv_temp[idx_test,]
    cifs=cif(t0)
    if(is.null(cifs)) return(NULL)
    p1.hat.m1=cifs[[1]]
    p2.hat.m1=cifs[[2]]
    p3.hat.m1=cifs[[3]]
    p1.hat.m2=cifs[[4]]
    p2.hat.m2=cifs[[5]]
    p3.hat.m2=cifs[[6]]
    #############################################################################
    #### Kaplan-Meier estimator for the censoring distribution, G(t)=pr(C>t) ####
    #############################################################################
    
    censor.S1= NULL
    censor.lambdastar = NULL
    censored.index = which(delta == 0)
    
    sort.time = sort(unique(X))
    
    if (length(censored.index)>0) ## handle the case there is no censoring
    {
      censor.tmp1 = matrix(rep(X,each=length(sort.time)),length(sort.time),length(idx_train))
      censor.tmp2 = matrix(rep(delta,each = length(sort.time)),length(sort.time),length(idx_train))
      
      censor.Y  = apply(ifelse(censor.tmp1 >= sort.time,1,0),1,sum)
      censor.d  = apply(ifelse(censor.tmp1 ==  sort.time & censor.tmp2 == 0,1,0),1,sum)
      censor.lambdastar = censor.d/censor.Y
      censor.S1 =  cumprod(1-censor.lambdastar)
    } else {
      censor.S1=rep(1,length(sort.time))
    }
    est.S.censoring=NULL
    for(idx in idx_test){
      if(Time[idx]<sort.time[1])
        est.S.censoring=c(est.S.censoring,1)
      else
        est.S.censoring=c(est.S.censoring,censor.S1[max(which(sort.time<=Time[idx]))])
    }
    # define a function to assign the subject to the class with the highest prob
    
    loc = function(tt)
    {which(tt==max(tt))[1]}
    # evaluate the NRI at each of time points
    estNRI = NULL
    estsd = NULL
    for (j in 1:npts){
      
      # for each subject, assign to one of the 3 categories based on the estimated risks of failing from cause 1, 2 or none
      
      pp.m1 = cbind(p1.hat.m1[,j],p2.hat.m1[,j], p3.hat.m1[,j])
      pind.m1 = apply(pp.m1,1,loc)
      
      pp.m2 = cbind(p1.hat.m2[,j],p2.hat.m2[,j], p3.hat.m2[,j])
      pind.m2 = apply(pp.m2,1,loc)
      
      
      num1 = ifelse(pind.m2==1 & pind.m1!=1 & Time[idx_test]<=t0[j] & Event[idx_test] ==1,1,0)/est.S.censoring-
        ifelse(pind.m2!=1 & pind.m1==1 & Time[idx_test]<=t0[j] & Event[idx_test] ==1,1,0)/est.S.censoring
      num1 = ifelse(num1=='NaN',0,num1)
      
      den1 = ifelse(Time[idx_test]<=t0[j] & Event[idx_test] ==1,1,0)/est.S.censoring
      den1 = ifelse(den1=='NaN',0,den1)
      ratio1 = sum(num1)/sum(den1)
      
      num2 = ifelse(pind.m2==2 & pind.m1!=2 & Time[idx_test]<=t0[j] & Event[idx_test] ==2,1,0)/est.S.censoring-
        ifelse(pind.m2!=2 & pind.m1==2 & Time[idx_test]<=t0[j] & Event[idx_test] ==2,1,0)/est.S.censoring
      num2 = ifelse(num2=='NaN',0,num2)
      
      den2 = ifelse(Time[idx_test]<=t0[j] & Event[idx_test] ==2,1,0)/est.S.censoring
      den2 = ifelse(den2=='NaN',0,den2)
      ratio2 = sum(num2)/sum(den2)
      
      num3 = ifelse(pind.m2==3 & pind.m1!=3 & Time[idx_test]>t0[j],1,0)-ifelse(pind.m2!=3 & pind.m1==3 & Time[idx_test]>t0[j],1,0)
      num3 = ifelse(num3=='NaN',0,num3)
      
      den3 = ifelse(Time[idx_test]>t0[j],1,0)
      den3 = ifelse(den3=='NaN',0,den3)
      ratio3 = sum(num3)/sum(den3)
      
      
      estNRI = c(estNRI, (ratio1+ratio2+ratio3)/3)
      
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
    pi.t1 = den.t1/n.t1
    
    
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
    pi.t2 = den.t2/n.t2
    
    
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
    pi.t3 = den.t3/n.t3
    
    
    pp.m2.t3 = cbind(p1.hat.m2[,3],p2.hat.m2[,3], p3.hat.m2[,3])
    pp.m1.t3 = cbind(p1.hat.m1[,3],p2.hat.m1[,3], p3.hat.m1[,3])
    
    num.t3.d1 = (pp.m1.t3 - matrix(rep(apply(pp.m1.t3,2,mean),ndim),ndim,3,byrow=TRUE))^2
    num.t3.d2 = (pp.m2.t3 - matrix(rep(apply(pp.m2.t3,2,mean),ndim),ndim,3,byrow=TRUE))^2
    
    num.t3 = apply(num.t3.d2-num.t3.d1,2,mean)
    estIDI.t3 = mean(num.t3/pi.t3/(1-pi.t3))
    output_final=rbind(estNRI,c(estIDI.t1, estIDI.t2, estIDI.t3))
    if(to_boot){
      
      idx_test=(1:ndim)[!(1:ndim)%in%idx_train]
      npts=length(t0)
      
      survival_test=surv_temp[idx_test,]
      cifs=cif(t0)
      if(is.null(cifs)) return(NULL)
      p1.hat.m1=cifs[[1]]
      p2.hat.m1=cifs[[2]]
      p3.hat.m1=cifs[[3]]
      p1.hat.m2=cifs[[4]]
      p2.hat.m2=cifs[[5]]
      p3.hat.m2=cifs[[6]]
      #############################################################################
      #### Kaplan-Meier estimator for the censoring distribution, G(t)=pr(C>t) ####
      #############################################################################
      
      censor.S1= NULL
      censor.lambdastar = NULL
      censored.index = which(delta == 0)
      
      sort.time = sort(unique(X))
      
      if (length(censored.index)>0) ## handle the case there is no censoring
      {
        censor.tmp1 = matrix(rep(X,each=length(sort.time)),length(sort.time),length(idx_train))
        censor.tmp2 = matrix(rep(delta,each = length(sort.time)),length(sort.time),length(idx_train))
        
        censor.Y  = apply(ifelse(censor.tmp1 >= sort.time,1,0),1,sum)
        censor.d  = apply(ifelse(censor.tmp1 ==  sort.time & censor.tmp2 == 0,1,0),1,sum)
        censor.lambdastar = censor.d/censor.Y
        censor.S1 =  cumprod(1-censor.lambdastar)
      } else {
        censor.S1=rep(1,length(sort.time))
      }
      est.S.censoring=NULL
      for(idx in idx_test){
        if(Time[idx]<sort.time[1])
          est.S.censoring=c(est.S.censoring,1)
        else
          est.S.censoring=c(est.S.censoring,censor.S1[max(which(sort.time<=Time[idx]))])
      }
      # define a function to assign the subject to the class with the highest prob
      
      loc = function(tt)
      {which(tt==max(tt))}
      # evaluate the NRI at each of time points
      estNRI = NULL
      estsd = NULL
      for (j in 1:npts){
        
        # for each subject, assign to one of the 3 categories based on the estimated risks of failing from cause 1, 2 or none
        
        pp.m1 = cbind(p1.hat.m1[,j],p2.hat.m1[,j], p3.hat.m1[,j])
        pind.m1 = apply(pp.m1,1,loc)
        
        pp.m2 = cbind(p1.hat.m2[,j],p2.hat.m2[,j], p3.hat.m2[,j])
        pind.m2 = apply(pp.m2,1,loc)
        
        
        num1 = ifelse(pind.m2==1 & pind.m1!=1 & Time[idx_test]<=t0[j] & Event[idx_test] ==1,1,0)/est.S.censoring-
          ifelse(pind.m2!=1 & pind.m1==1 & Time[idx_test]<=t0[j] & Event[idx_test] ==1,1,0)/est.S.censoring
        num1 = ifelse(num1=='NaN',0,num1)
        
        den1 = ifelse(Time[idx_test]<=t0[j] & Event[idx_test] ==1,1,0)/est.S.censoring
        den1 = ifelse(den1=='NaN',0,den1)
        ratio1 = sum(num1)/sum(den1)
        
        num2 = ifelse(pind.m2==2 & pind.m1!=2 & Time[idx_test]<=t0[j] & Event[idx_test] ==2,1,0)/est.S.censoring-
          ifelse(pind.m2!=2 & pind.m1==2 & Time[idx_test]<=t0[j] & Event[idx_test] ==2,1,0)/est.S.censoring
        num2 = ifelse(num2=='NaN',0,num2)
        
        den2 = ifelse(Time[idx_test]<=t0[j] & Event[idx_test] ==2,1,0)/est.S.censoring
        den2 = ifelse(den2=='NaN',0,den2)
        ratio2 = sum(num2)/sum(den2)
        
        num3 = ifelse(pind.m2==3 & pind.m1!=3 & Time[idx_test]>t0[j],1,0)-ifelse(pind.m2!=3 & pind.m1==3 & Time[idx_test]>t0[j],1,0)
        num3 = ifelse(num3=='NaN',0,num3)
        
        den3 = ifelse(Time[idx_test]>t0[j],1,0)
        den3 = ifelse(den3=='NaN',0,den3)
        ratio3 = sum(num3)/sum(den3)
        
        
        estNRI = c(estNRI, (ratio1+ratio2+ratio3)/3)
        
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
      pi.t1 = den.t1/n.t1
      
      
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
      pi.t2 = den.t2/n.t2
      
      
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
      pi.t3 = den.t3/n.t3
      
      
      pp.m2.t3 = cbind(p1.hat.m2[,3],p2.hat.m2[,3], p3.hat.m2[,3])
      pp.m1.t3 = cbind(p1.hat.m1[,3],p2.hat.m1[,3], p3.hat.m1[,3])
      
      num.t3.d1 = (pp.m1.t3 - matrix(rep(apply(pp.m1.t3,2,mean),ndim),ndim,3,byrow=TRUE))^2
      num.t3.d2 = (pp.m2.t3 - matrix(rep(apply(pp.m2.t3,2,mean),ndim),ndim,3,byrow=TRUE))^2
      
      num.t3 = apply(num.t3.d2-num.t3.d1,2,mean)
      estIDI.t3 = mean(num.t3/pi.t3/(1-pi.t3))
     
    output_final=rbind(output_final,estNRI,c(estIDI.t1, estIDI.t2, estIDI.t3))
    }
    return(output_final)
    
  }
  
  
  # number of simulation
  nsim = 1
  # sample size
  ndim = 400
  nboot = 1000
  
  
  ### true parameters ###
  beta=c(1,1,3)
  LC=c(9,5,-11,-11,-11)
  RC=c(40,40,22,27,8.2)
  per=c(10,20,30,40,50)
  t0 = c(-2,-1,0)	
  res=NULL
  for(ww in c(3,5)){
    
    for ( m in 1:nsim)	# number of simulations
    {
      ### generate the covariates ###
      
      z1 = rnorm(ndim)					# continuous 
      z1 = pmax(pmin(z1,3.5),-3.5)
      z2 = ifelse(runif(ndim,0,1)<=0.7,1,0)		# 0,1 variable
      z3 = rnorm(ndim)					# continuous with large coefficient
      z3 = pmax(pmin(z3,3.5),-3.5)
      beta=c(1,1,3)
      t=(log((1/runif(ndim)-1)/2)-cbind(z1,z2,z3)%*%beta)
      TT=pmax(-100,pmin(t,100))
      
      epsilon = ifelse(runif(ndim)<0.5,1,2)
      C = runif(ndim, LC[ww], RC[ww])		
      
      X = as.numeric(pmin(TT,C))
      delta = as.numeric(TT<=C)*epsilon
      
      ests=try(estNRI_Gerds(data.frame(time=X,event=delta,Z1=z1,Z2=z2,Z3=z3), t0),silent=TRUE)
      if(class(ests)=="try-error") return(NULL)
      NRI = ests[1,]
      IDI = ests[2,]
      
      
      boot.NRI = NULL
      zheng.NRI = NULL
      boot.IDI = NULL
      zheng.IDI = NULL
      for(i in 1:nboot)
      {
        zanshi=try(estNRI_Gerds(data.frame(time=X,event=delta,Z1=z1,Z2=z2,Z3=z3), t0, TRUE),silent=TRUE)
        if(class(zanshi)=="try-error") next
        boot.NRI = rbind(boot.NRI,t(zanshi[1,]))
        zheng.NRI = rbind(zheng.NRI,t(zanshi[3,]))
        boot.IDI = rbind(boot.NRI,t(zanshi[2,]))
        zheng.IDI = rbind(zheng.NRI,t(zanshi[4,]))
      }
      
      
      
      boot.NRI = na.omit(boot.NRI)
      zheng.NRI = na.omit(zheng.NRI)
      boot.IDI = na.omit(boot.IDI)
      zheng.IDI = na.omit(zheng.IDI)
      if(nrow(boot.NRI)<=nboot/2 | nrow(zheng.NRI)<=nboot/2) return(NULL)
      
      sd.NRI = t(apply(boot.NRI,2,sd))
      
      NRI_zheng=0.632*colMeans(zheng.NRI)+(1-0.632)*NRI
      bias=NRI-NRI_zheng
      
      covs_NRI=(NRI+qnorm(0.025)*sd.NRI-bias < c(0.2090376,0.1886032,0.1688827))&(NRI+qnorm(0.975)*sd.NRI-bias > c(0.2090376,0.1886032,0.1688827))
      
      
      sd.IDI = t(apply(boot.IDI,2,sd))
      
      IDI_zheng=0.632*colMeans(zheng.IDI)+(1-0.632)*IDI
      bias=IDI-IDI_zheng
      
      covs_IDI=(IDI+qnorm(0.025)*sd.IDI-bias < c(0.2810522,0.2672072,0.2404910))&(IDI+qnorm(0.975)*sd.IDI-bias > c(0.2810522,0.2672072,0.2404910))
      
      res=rbind(res,rbind(NRI,sd.NRI,NRI_zheng,covs_NRI),rbind(IDI,sd.IDI,IDI_zheng,covs_IDI))
    }
    
  }
  return(res)
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
NRIIDIGerdsUsingGerds400_Alternative_zheng <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(NRIIDIGerdsUsingGerds400_Alternative_zheng,file = 'NRIIDIGerdsUsingGerds400_Alternative_zheng.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF