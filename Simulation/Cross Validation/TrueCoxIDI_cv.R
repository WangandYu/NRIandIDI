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
sfLibrary(purrr)

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
  
  
  npts=3
  ### generate the covariates ###
  z10 = rnorm(ndim)					# continuous 
  z10 = pmax(pmin(z10,3.5),-3.5)
  z20 = ifelse(runif(ndim,0,1)<=0.7,1,0)		# 0,1 variable
  z30 = rnorm(ndim)					# continuous with large coefficient
  z30 = pmax(pmin(z30,3.5),-3.5)
  
  ## generate event time using Weibull model ##
  z0 = cbind(z10,z20,z30)					# covariates
  U0 = runif(ndim, 0, 1)				# uniform(0,1)
  w0 = log(-log(1-U0))					# inverse cdf of extreme value dist
  TT0 = exp(beta0 + z0%*%beta + sigma*w0)		# Weibull model
  
  
  epsilon0 = ifelse(runif(ndim)<0.5,1,2)
  
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
    data1 = data.frame(cbind(TT, event=ifelse(epsilon==1,1,0), z1=z1, z2=z2))
    cox1 = coxph(Surv(TT, event) ~ z1 + z2, data1)	
    # summary(cox1)
    cox.basehaz1 = basehaz(cox1, centered=FALSE)
    cox.basehaz1 = cox.basehaz1[cox.basehaz1[['time']]<=t0,]	# only keep the baseline hazard until time t0
    
    # Then calculate the comulative baseline hazard fuction for cause 2 #
    data2 = data.frame(cbind(TT, event=ifelse(epsilon==2,1,0), z1=z1, z2=z2))
    cox2 = coxph(Surv(TT, event) ~ z1 + z2, data2)	
    # summary(cox2)
    cox.basehaz2 = basehaz(cox2, centered=FALSE)
    cox.basehaz2 = cox.basehaz2[cox.basehaz2[['time']]<=t0,]	# only keep the baseline hazard until time t0[2]
    
    # indecies for the predified t0 #
    
    t0.index1 = max(which(cox.basehaz1[['time']]<=T0[1]))
    t0.index2 = max(which(cox.basehaz1[['time']]<=T0[2]))
    t0.index3 = max(which(cox.basehaz1[['time']]<=T0[3]))
    t0.index = c(t0.index1, t0.index2, t0.index3)
    
    dim.uniquetime = length(which(cox.basehaz1[['time']]<=T0[npts]))	
    
    
    # estimate the cumulative cause specific hazard for cause 1 #
    cox.haz1 = matrix(rep(cox.basehaz1[['hazard']],each=length(layer)),length(layer),dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
    cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z10[layer],z20[layer])%*%coef(cox1)),dim.uniquetime),length(layer),dim.uniquetime)
    # calcualte the cumulative cause specific hazard for cause 2 #
    cox.haz2 = matrix(rep(cox.basehaz2[['hazard']],each=length(layer)),length(layer),dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
    cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z10[layer],z20[layer])%*%coef(cox2)),dim.uniquetime),length(layer),dim.uniquetime)
    
    # estimate the overall survival function #
    S = exp(-cox.haz1-cox.haz2)
    S = cbind(1,S[,1:(dim(S)[2]-1)])		# S(t-;z)
    # dLemda1(t;z)
    dLamda1 = cox.haz1 - cbind(0,cox.haz1[,1:(dim(cox.haz1)[2]-1)])
    # estimtate cumulative incidence function for cause1 #
    F1 = t(apply(S*dLamda1,1,cumsum))
    # estimates of F1 at t0 #
    p1.hat.m1 = rbind(p1.hat.m1, F1[,t0.index])
    
    dLamda2 = cox.haz2 - cbind(0,cox.haz2[,1:(dim(cox.haz2)[2]-1)])
    # estimtate cumulative incidence function for cause1 #
    F2 = t(apply(S*dLamda2,1,cumsum))
    # estimates of F1 at t0 #
    p2.hat.m1 = rbind(p2.hat.m1, F2[,t0.index])
    
    
    
    #######################################################################################
    ### Fit the cox proportional hazard regression model using all covariates (model 2) ###
    #######################################################################################
    
    # calcualte the cumulative baseline hazard function for cause 1 #
    data1 = data.frame(cbind(TT, event=ifelse(epsilon==1,1,0), z1=z1, z2=z2, z3=z3))
    cox1 = coxph(Surv(TT, event) ~ z1 + z2 + z3, data1)	
    #summary(cox1)
    cox.basehaz1 = basehaz(cox1, centered=FALSE)
    cox.basehaz1 = cox.basehaz1[cox.basehaz1[['time']]<=t0,]	# only keep the baseline hazard until time t0
    
    # Then calculate the comulative baseline hazard fuction for cause 2 #
    data2 = data.frame(cbind(TT, event=ifelse(epsilon==2,1,0), z1=z1, z2=z2, z3=z3))
    cox2 = coxph(Surv(TT, event) ~ z1 + z2 + z3, data2)	
    #summary(cox2)
    cox.basehaz2 = basehaz(cox2, centered=FALSE)
    cox.basehaz2 = cox.basehaz2[cox.basehaz2[['time']]<=t0,]	# only keep the baseline hazard until time t0
    
    # index for cause 1 and 2 events at the predified t0 #
    t0.index1 = max(which(cox.basehaz1[['time']]<=T0[1]))
    t0.index2 = max(which(cox.basehaz1[['time']]<=T0[2]))
    t0.index3 = max(which(cox.basehaz1[['time']]<=T0[3]))
    t0.index = c(t0.index1, t0.index2, t0.index3)
    
    dim.uniquetime = length(which(cox.basehaz1[['time']]<=T0[npts]))	
    
    # estimate the cumulative cause specific hazard for cause 1 #
    cox.haz1 = matrix(rep(cox.basehaz1[['hazard']],each=length(layer)),length(layer),dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
    cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z10[layer],z20[layer],z30[layer])%*%coef(cox1)),dim.uniquetime),length(layer),dim.uniquetime)
    
    # calcualte the cumulative cause specific hazard for cause 2 #
    cox.haz2 = matrix(rep(cox.basehaz2[['hazard']],each=length(layer)),length(layer),dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
    cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z10[layer],z20[layer],z30[layer])%*%coef(cox2)),dim.uniquetime),length(layer),dim.uniquetime)
    
    # estimate the overall survial fuction #
    S = exp(-cox.haz1-cox.haz2)
    S = cbind(1,S[,1:(dim(S)[2]-1)])		# S(t-;z)
    # dLamda1(t;z)
    dLamda1 = cox.haz1 - cbind(0,cox.haz1[,1:(dim(cox.haz1)[2]-1)])
    # esimtate cumulative incidence function for cause1 #
    F1 = t(apply(S*dLamda1,1,cumsum))
    # estimates of F1 at t0 #
    p1.hat.m2 = rbind(p1.hat.m2, F1[,t0.index])
    
    # dLamda2(t;z)
    dLamda2 = cox.haz2 - cbind(0,cox.haz2[,1:(dim(cox.haz2)[2]-1)])
    # esimtate cumulative incidence function for cause1 #
    F2 = t(apply(S*dLamda2,1,cumsum))
    
    # estimates of F2 at t0 #
    p2.hat.m2 = rbind(p2.hat.m2, F2[,t0.index])
  }
  
  p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
  
  p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
  
  
  sums=NULL
  for(i in 1:3){
    sum1 = sum((p1.hat.m2[,i]-mean(p1.hat.m2[,i]))^2 - (p1.hat.m1[,i]-mean(p1.hat.m1[,i]))^2)
    sum2 = sum((p2.hat.m2[,i]-mean(p2.hat.m2[,i]))^2 - (p2.hat.m1[,i]-mean(p2.hat.m1[,i]))^2)
    sum3 = sum((p3.hat.m2[,i]-mean(p3.hat.m2[,i]))^2 - (p3.hat.m1[,i]-mean(p3.hat.m1[,i]))^2)
    
    sums=rbind(sums,t(c(sum1,sum2,sum3)))
  }
  return(rbind(Ns,sums))
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
TrueIDI_WB_CV <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(TrueIDI_WB_CV,file = 'TrueIDI_WB_CV.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF
