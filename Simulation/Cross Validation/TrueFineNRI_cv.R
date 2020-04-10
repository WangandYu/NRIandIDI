#!/bin/bash
#SBATCH --job-name=NRICF400
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
sfLibrary(survival)
sfLibrary(cmprsk)

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
  ndim = 1000
  nsim = 1000
  
  ### true parameters ###
  beta11=0.2
  beta12=-0.5
  beta13=1
  p=0.65
  
  T0 = c(20,21,22)
  
  
  z10=rnorm(ndim)
  z10=pmax(pmin(z10,3.5),-3.5)
  z20=ifelse(runif(ndim)<0.7,1,0)
  z30=rnorm(ndim)
  z30=pmax(pmin(z30,3.5),-3.5)
  
  F10=1-(1-p)^(exp(z10*beta11+z20*beta12+z30*beta13))
  
  W0=runif(ndim)
  epsilon0=ifelse(W0<F10,1,2)
  TT0=rep(0,ndim)
  TT0[W0<F10]=-log(1-(1-(1-W0[W0<F10])^(1/exp(z10[W0<F10]*beta11+z20[W0<F10]*beta12+z30[W0<F10]*beta13)))/p)
  zz0=exp(z10[W0>=F10]*beta11/10+z20[W0>=F10]*beta12/10+z30[W0>=F10]*beta13/10)
  TT0[W0>=F10]=sapply(zz0,function(o) return(rexp(1,o)))
  TT0=TT0^0.2
  TT0=pmin(TT0*20,100)
  
  Ns=NULL
  for(t0 in T0){
    n1 = sum(ifelse(TT0<=t0 & epsilon0==1,1,0))
    n2 = sum(ifelse(TT0<=t0 & epsilon0==2,1,0))
    n3 = ndim - n1 - n2
    
    Ns=rbind(Ns,t(c(n1,n2,n3)))
  }
  cv=split(1:ndim,sort((1:ndim)%%5))
  p1.hat.m1=p2.hat.m1=p3.hat.m1=p1.hat.m2=p2.hat.m2=p3.hat.m2=NULL
  for(layer in cv){
    z1=z10[-layer]
    z2=z20[-layer]
    z3=z30[-layer]
    TT=TT0[-layer]
    epsilon=epsilon0[-layer]
    #############################################################################
    ### Fit the cox proportional hazard regression model without z3 (model 1) ###
    #############################################################################
    
    # calcualte the cumulative baseline hazard function for cause 1 #
    cov = cbind(z1,z2)
    crr1 = crr(TT,epsilon,cov)
    pred1 = predict(crr1,cbind(z10[layer],z20[layer]))
    # summary(crr1)
    
    # indecies for the predified t0 #
    t0.index= c(max(which(pred1[,1]<T0[1])),max(which(pred1[,1]<T0[2])),max(which(pred1[,1]<T0[3])))
    
    # estimate the overall survival function #
    p1.hat.m1 = cbind(p1.hat.m1,(pred1[t0.index,-1]))
    
    crr2 = crr(TT,epsilon,cov,failcode=2)
    pred2 = predict(crr2,cbind(z10[layer],z20[layer]))
    # summary(crr2)
    
    # indecies for the predified t0 #
    t0.index= c(max(which(pred2[,1]<T0[1])),max(which(pred2[,1]<T0[2])),max(which(pred2[,1]<T0[3])))
    
    # estimate the overall survival function #
    p2.hat.m1 = cbind(p2.hat.m1, (pred2[t0.index,-1]))
    
    
    
    #######################################################################################
    ### Fit the cox proportional hazard regression model using all covariates (model 2) ###
    #######################################################################################
    
    # calcualte the cumulative baseline hazard function for cause 1 #
    cov = cbind(z1,z2,z3)
    crr1 = crr(TT,epsilon,cov,failcode=1)
    pred1 = predict(crr1,cbind(z10[layer],z20[layer],z30[layer]))
    # summary(crr1)
    
    # indecies for the predified t0 #
    t0.index= c(max(which(pred1[,1]<T0[1])),max(which(pred1[,1]<T0[2])),max(which(pred1[,1]<T0[3])))
    
    # estimate the overall survival function #
    p1.hat.m2 = cbind(p1.hat.m2,pred1[t0.index,-1])
    
    crr2 = crr(TT,epsilon,cov,failcode=2)
    pred2 = predict(crr2,cbind(z10[layer],z20[layer],z30[layer]))
    # summary(crr2)
    
    # indecies for the predified t0 #
    t0.index= c(max(which(pred2[,1]<T0[1])),max(which(pred2[,1]<T0[2])),max(which(pred2[,1]<T0[3])))
    
    # estimate the overall survival function #
    p2.hat.m2 = cbind(p2.hat.m2,pred2[t0.index,-1])
  } 
  p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
  
  p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
  
  
  
  pind.m1=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m1[i,],p2.hat.m1[i,], p3.hat.m1[i,])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m1 = cbind(pind.m1,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  pind.m2=NULL
  for(i in 1:3){
    pp = cbind(p1.hat.m2[i,],p2.hat.m2[i,], p3.hat.m2[i,])
    pmax = apply(pp,1,max)
    pmax.matrix = cbind(pmax,pmax,pmax)
    ind.matrix = cbind(rep(1,ndim),rep(2,ndim),rep(3,ndim))
    pind.m2 = cbind(pind.m2,apply((pp>=pmax.matrix)*ind.matrix,1,max))
  }
  
  sums=NULL
  for(i in 1:3){
    t0=T0[i]
    sum1 = sum(ifelse(pind.m2[,i]==1 & pind.m1[,i]!=1 & TT0<=t0 & epsilon0 ==1, 1, 0)) - sum(ifelse(pind.m2[,i]!=1 & pind.m1[,i]==1 & TT0<=t0 & epsilon0 ==1, 1, 0))
    sum2 = sum(ifelse(pind.m2[,i]==2 & pind.m1[,i]!=2 & TT0<=t0 & epsilon0 ==2, 1, 0)) - sum(ifelse(pind.m2[,i]!=2 & pind.m1[,i]==2 & TT0<=t0 & epsilon0 ==2, 1, 0))
    sum3 = sum(ifelse(pind.m2[,i]==3 & pind.m1[,i]!=3 & TT0>t0, 1, 0)) - sum(ifelse(pind.m2[,i]!=3 & pind.m1[,i]==3 & TT0>t0, 1, 0))
    
    sums=rbind(sums,t(c(sum1,sum2,sum3)))
  }
  return(rbind(Ns,sums))
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
TrueNRI_Fine_CV <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(TrueNRI_Fine_CV,file = 'TrueNRI_Fine_CV.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF
