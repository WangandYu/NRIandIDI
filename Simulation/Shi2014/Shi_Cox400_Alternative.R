#!/bin/bash
#SBATCH --job-name=CoxShi
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
sfLibrary(survival)

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
  
  Shi_Cox = function(X,delta, z1, z2,z3, t0)
  {ndim = length(X)
  npts = length(t0)
  
  #############################################################################
  ### Fit the cox proportional hazard regression model without z3 (model 1) ###
  #############################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  data1 = data.frame(cbind(X, event=ifelse(delta==1,1,0), z1=z1,z2=z2))
  names(data1) = c("X","event","z1","z2")
  cox1 = coxph(Surv(X, event) ~ z1+z2, data1)	
  # summary(cox1)
  cox.basehaz1 = basehaz(cox1, centered=FALSE)
  cox.basehaz1 = cox.basehaz1[cox.basehaz1[['time']]<=t0[npts],]	# only keep the baseline hazard until the last point of t0
  
  # Then calculate the comulative baseline hazard fuction for cause 2 #
  data2 = data.frame(cbind(X, event=ifelse(delta==2,1,0), z1=z1,z2=z2))
  names(data2) = c("X","event","z1","z2")
  cox2 = coxph(Surv(X, event) ~ z1+z2, data2)	
  # summary(cox2)
  cox.basehaz2 = basehaz(cox2, centered=FALSE)
  cox.basehaz2 = cox.basehaz2[cox.basehaz2[['time']]<=t0[npts],]	# only keep the baseline hazard until the last point of t0
  
  # indecies for the predified t0 #
  t0.index1 = max(which(cox.basehaz1[['time']]<=t0[1]))
  t0.index2 = max(which(cox.basehaz1[['time']]<=t0[2]))
  t0.index3 = max(which(cox.basehaz1[['time']]<=t0[3]))
  t0.index = c(t0.index1, t0.index2, t0.index3)
  
  
  
  dim.uniquetime = length(which(cox.basehaz1[['time']]<=t0[npts]))	
  
  
  # estimate the cumulative cause specific hazard for cause 1 #
  cox.haz1 = matrix(rep(cox.basehaz1[['hazard']],each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
  cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z1,z2)%*%coef(cox1)),dim.uniquetime),ndim,dim.uniquetime)
  # calcualte the cumulative cause specific hazard for cause 2 #
  cox.haz2 = matrix(rep(cox.basehaz2[['hazard']],each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
  cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z1,z2)%*%coef(cox2)),dim.uniquetime),ndim,dim.uniquetime)
  
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
  
  p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
  
  
  
  
  
  #######################################################################################
  ### Fit the cox proportional hazard regression model using all covariates (model 2) ###
  #######################################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  data1 = data.frame(cbind(X, event=ifelse(delta==1,1,0), z1=z1, z2=z2,z3=z3))
  names(data1) = c("X","event","z1","z2","z3")
  cox1 = coxph(Surv(X, event) ~ z1 + z2+z3, data1)	
  #summary(cox1)
  cox.basehaz1 = basehaz(cox1, centered=FALSE)
  cox.basehaz1 = cox.basehaz1[cox.basehaz1[['time']]<=t0[npts],]	# only keep the baseline hazard until time t0
  
  # Then calculate the comulative baseline hazard fuction for cause 2 #
  data2 = data.frame(cbind(X, event=ifelse(delta==2,1,0), z1=z1, z2=z2,z3=z3))
  names(data2) = c("X","event","z1","z2","z3")
  cox2 = coxph(Surv(X, event) ~ z1 + z2+z3, data2)	
  #summary(cox2)
  cox.basehaz2 = basehaz(cox2, centered=FALSE)
  cox.basehaz2 = cox.basehaz2[cox.basehaz2[['time']]<=t0[npts],]	# only keep the baseline hazard until time t0
  
  # indecies for the predified t0 #
  t0.index1 = max(which(cox.basehaz1[['time']]<=t0[1]))
  t0.index2 = max(which(cox.basehaz1[['time']]<=t0[2]))
  t0.index3 = max(which(cox.basehaz1[['time']]<=t0[3]))
  t0.index = c(t0.index1, t0.index2, t0.index3)
  
  dim.uniquetime = length(which(cox.basehaz1[['time']]<=t0[npts]))
  
  # estimate the cumulative cause specific hazard for cause 1 #
  cox.haz1 = matrix(rep(cox.basehaz1[['hazard']],each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
  cox.haz1 = cox.haz1*matrix(rep(exp(cbind(z1,z2,z3)%*%coef(cox1)),dim.uniquetime),ndim,dim.uniquetime)
  
  # calcualte the cumulative cause specific hazard for cause 2 #
  cox.haz2 = matrix(rep(cox.basehaz2[['hazard']],each=ndim),ndim,dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
  cox.haz2 = cox.haz2*matrix(rep(exp(cbind(z1,z2,z3)%*%coef(cox2)),dim.uniquetime),ndim,dim.uniquetime)
  
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
  p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
  
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
  beta0 = 2.5
  beta1 = 0.05
  beta2 = -0.05
  beta3 = 0.15
  sigma = 0.2
  
  beta = c(beta1,beta2,beta3)			# regression coefficients
  
  t0 = c(11,12,13) 	
  
  
  LC=c(9,7,2,3,1)
  RC=c(38,30,31,25,21)
  per=c(10,20,30,40,50)
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
      
      ## generate event time using Weibull model ##
      z = cbind(z1,z2,z3)					# covariates
      U = runif(ndim, 0, 1)				# uniform(0,1)
      w = log(-log(1-U))					# inverse cdf of extreme value dist
      TT = exp(beta0 + z%*%beta + sigma*w)		# Weibull model
      
      
      epsilon = ifelse(runif(ndim)<0.5,1,2)
      
      
      ## uniformly generate censoring time ##
      
      C = runif(ndim, LC[ww], RC[ww])		
      
      X = as.numeric(pmin(TT,C))
      delta = as.numeric(TT<=C)*epsilon
      
      Shi=try(Shi_Cox(X,delta,z1, z2,z3, t0),silent=TRUE)
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
estShi_Cox400_Alternative <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(estShi_Cox400_Alternative,file = 'estShi_Cox400_Alternative.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF