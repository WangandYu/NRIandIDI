#!/bin/bash
#SBATCH --job-name=GerdsShi
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=25
#SBATCH --cluster=mpi
#SBATCH --partition=opa-high-mem
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
  
  
  cummatrix = function(ee)
  {n2 = dim(ee)[1]
  n4 = dim(ee)[2]
  
  for (i in 1:n4)
    ee[,i] = cumsum(ee[,i])
  
  for (i in 1:n2)
    ee[i,] = cumsum(ee[i,])
  return(ee)
  }
  
  cumprodmatrix = function(ee)
  {n2 = dim(ee)[1]
  n4 = dim(ee)[2]
  
  for (i in 1:n4)
    ee[,i] = cumprod(ee[,i])
  
  for (i in 1:n2)
    ee[i,] = cumprod(ee[i,])
  return(ee)
  }
  
  
  est_Dabrowska = function(DD, timept1, timept2)
  { numtime = length(timept1)
  out = rep(-1,numtime*2)
  ndim = dim(DD)[1]        
  uniquetime1 = unique(sort(DD[,1]))
  uniquetime2 = unique(sort(DD[,2]))
  
  npts1 = length(uniquetime1)
  npts2 = length(uniquetime2)
  
  tmp1 = matrix(rep(DD[,1],each = npts1),npts1,ndim)
  tmp3 = matrix(rep(DD[,3],each = npts1),npts1,ndim)
  tmp2 = matrix(rep(DD[,2],each = npts2),npts2,ndim)
  tmp4 = matrix(rep(DD[,4],each = npts2),npts2,ndim)
  
  
  
  numrisk = ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))
  lambda11 = ifelse(tmp1==uniquetime1 & tmp3==1,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4==1,1,0))/numrisk
  
  lambda10.left =ifelse(tmp1==uniquetime1 & tmp3 ==1,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))/numrisk   
  lambda01.up = ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4==1,1,0))/numrisk
  
  numrisk = NULL
  
  lambda10.left = ifelse(lambda10.left == "NaN", 0, lambda10.left)
  lambda01.up = ifelse(lambda01.up == "NaN", 0, lambda01.up)
  lambda11 = ifelse(lambda11 == "NaN", 0, lambda11)
  
  
  
  lambda1 = apply(ifelse(tmp1==uniquetime1 & tmp3 == 1,1,0),1,sum)/apply(ifelse(tmp1>=uniquetime1,1,0),1,sum)   
  S1 = cumprod(1-lambda1)
  lambda2 = apply(ifelse(tmp2==uniquetime2 & tmp4 == 1,1,0),1,sum)/apply(ifelse(tmp2>=uniquetime2,1,0),1,sum)   
  S2 = cumprod(1-lambda2)
  tmp1 = tmp2 = tmp3 = tmp4  = NULL
  
  factor = ifelse(lambda10.left!=1 & lambda01.up!=1, 1-(lambda10.left*lambda01.up-lambda11)/(1-lambda10.left)/(1-lambda01.up), 1)
  
  
  factor = cumprodmatrix(factor)   
  factor = ifelse(is.na(factor),0,factor)
  
  SS = S1%*%t(S2)*factor
  
  factor = NULL
  
  SS = apply(SS,2,cummin)
  SS = t(apply(SS,1,cummin))
  
  S1.mat =  matrix(rep(S1,npts2),npts1,npts2)
  S2.mat =  matrix(rep(S2,each=npts1),npts1,npts2)
  
  F.mat = (1-S1.mat-S2.mat+SS)/(1-S2.mat)
  F.mat = ifelse(F.mat>=0 & F.mat<=1,F.mat, ifelse(F.mat<0,0,1))
  F.mat = apply(F.mat,2,cummax)
  
  G.mat = (S2.mat - SS)/S2.mat
  G.mat = ifelse(G.mat>=0 & G.mat<=1,G.mat, ifelse(G.mat<0,0,1))
  G.mat = apply(G.mat,2,cummax)
  
  for(i in 1:numtime)
  { loc1= max((uniquetime1<=timept1[i])*(1:npts1))
  loc2= max((uniquetime2<=timept2[i])*(1:npts2))
  out[i] = F.mat[loc1,loc2]
  out[numtime+i] = G.mat[loc1,loc2]
  }
  
  return(out)
  
  }
  
  
  Shi_Gerds = function(dt, t0)
  {
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
      dt=survival
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
    
    
    npts=length(t0)
    
    cifs=cif(t0)
    if(is.null(cifs)) return(NULL)
    p1.hat.m1=cifs[[1]]
    p2.hat.m1=cifs[[2]]
    p3.hat.m1=cifs[[3]]
    p1.hat.m2=cifs[[4]]
    p2.hat.m2=cifs[[5]]
    p3.hat.m2=cifs[[6]]
    delta=dt[['event']]
    X=dt[['time']]
    
    IDI1=IDI2=IAUC1=IAUC2=NULL
    for(point in 1:length(t0)){
      s=sort(unique(p1.hat.m2[,point]-p1.hat.m1[,point]))
      n.s=length(s)
      DDorg = cbind(p1.hat.m2[,point]-p1.hat.m1[,point],X,1,ifelse(delta==1,1,0))
      est.mean = est_Dabrowska(DDorg,s,rep(t0[point],length(s)))
      F.hat.dab = est.mean[1:n.s]
      G.hat.dab = est.mean[(n.s+1):(2*n.s)]
      
      # IDI is the area between F(s) and G(s) curves
      AUC.G.dab = sum((s[2:n.s]-s[1:(n.s-1)])*(G.hat.dab[2:n.s]+G.hat.dab[1:(n.s-1)])/2)
      AUC.F.dab = sum((s[2:n.s]-s[1:(n.s-1)])*(F.hat.dab[2:n.s]+F.hat.dab[1:(n.s-1)])/2)
      IDI1 = c(IDI1,AUC.G.dab - AUC.F.dab)	
      
      IAUC.index = max(which(s<=0))
      # IAUC is the vertical distance of F(0) and G(0)
      IAUC1 = c(IAUC1,G.hat.dab[IAUC.index]-F.hat.dab[IAUC.index])
      
      
      s=sort(unique(p2.hat.m2[,point]-p2.hat.m1[,point]))
      n.s=length(s)
      DDorg = cbind(p2.hat.m2[,point]-p2.hat.m1[,point],X,1,ifelse(delta==2,1,0))
      est.mean = est_Dabrowska(DDorg,s,rep(t0[point],length(s)))
      F.hat.dab = est.mean[1:n.s]
      G.hat.dab = est.mean[(n.s+1):(2*n.s)]
      
      # IDI is the area between F(s) and G(s) curves
      AUC.G.dab = sum((s[2:n.s]-s[1:(n.s-1)])*(G.hat.dab[2:n.s]+G.hat.dab[1:(n.s-1)])/2)
      AUC.F.dab = sum((s[2:n.s]-s[1:(n.s-1)])*(F.hat.dab[2:n.s]+F.hat.dab[1:(n.s-1)])/2)
      IDI2 = c(IDI2,AUC.G.dab - AUC.F.dab)	
      
      IAUC.index = max(which(s<=0))
      # IAUC is the vertical distance of F(0) and G(0)
      IAUC2 = c(IAUC2,G.hat.dab[IAUC.index]-F.hat.dab[IAUC.index])
    }
    
    return(c(IDI1,IDI2,IAUC1,IAUC2))
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
      
      Shi=try(Shi_Gerds(data.frame(time=X,event=delta,Z1=z1,Z2=z2,Z3=z3), t0),silent=TRUE)
      if(class(Shi)=="try-error") return(NULL)
      
      
      res=rbind(res,Shi)
    }
    
  }
  return(res)
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
estShi_Gerds400_Alternative <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(estShi_Gerds400_Alternative,file = 'estShi_Gerds400_Alternative.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF