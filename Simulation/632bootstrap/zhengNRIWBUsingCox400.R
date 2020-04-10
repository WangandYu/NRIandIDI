#!/bin/bash
#SBATCH --job-name=NRICC400
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
  
  estNRI_Weibull = function(Time,Event, Z1, Z2, Z3, t0, to_boot=FALSE)
  {
    ndim = length(Time)
    npts = length(t0)
    
    if(to_boot){
      idx_train=sample(1:ndim, ndim, replace=TRUE)
      idx_test=(1:ndim)[!(1:ndim)%in%idx_train]
    } else {
      idx_train=idx_test=1:ndim
    }
    
    #############################################################################
    ### Fit the cox proportional hazard regression model without z3 (model 1) ###
    #############################################################################
    X=Time[idx_train]
    delta=Event[idx_train]
    z1=Z1[idx_train]
    z2=Z2[idx_train]
    z3=Z3[idx_train]
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
    cox.haz1 = matrix(rep(cox.basehaz1[['hazard']],each=length(idx_test)),length(idx_test),dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
    cox.haz1 = cox.haz1*matrix(rep(exp(cbind(Z1[idx_test],Z2[idx_test])%*%coef(cox1)),dim.uniquetime),length(idx_test),dim.uniquetime)
    # calcualte the cumulative cause specific hazard for cause 2 #
    cox.haz2 = matrix(rep(cox.basehaz2[['hazard']],each=length(idx_test)),length(idx_test),dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
    cox.haz2 = cox.haz2*matrix(rep(exp(cbind(Z1[idx_test],Z2[idx_test])%*%coef(cox2)),dim.uniquetime),length(idx_test),dim.uniquetime)
    
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
    cox.haz1 = matrix(rep(cox.basehaz1[['hazard']],each=length(idx_test)),length(idx_test),dim.uniquetime)	# baseline cumulative hazard matrix for cause 1
    cox.haz1 = cox.haz1*matrix(rep(exp(cbind(Z1[idx_test],Z2[idx_test],Z3[idx_test])%*%coef(cox1)),dim.uniquetime),length(idx_test),dim.uniquetime)
    
    # calcualte the cumulative cause specific hazard for cause 2 #
    cox.haz2 = matrix(rep(cox.basehaz2[['hazard']],each=length(idx_test)),length(idx_test),dim.uniquetime)	# baseline cumulative hazard matrix for cause 2
    cox.haz2 = cox.haz2*matrix(rep(exp(cbind(Z1[idx_test],Z2[idx_test],Z3[idx_test])%*%coef(cox2)),dim.uniquetime),length(idx_test),dim.uniquetime)
    
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
    
    p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
    
    p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
    
    
    
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
    
    return(estNRI)
    
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
      
      ests=try(estNRI_Weibull(X,delta,z1, z2,z3, t0),silent=TRUE)
      if(class(ests)=="try-error") return(NULL)
      NRI = t(ests)
      
      
      boot.NRI = NULL
      for(i in 1:nboot)
      {
        index = sample(1:ndim, ndim, replace=TRUE)
        X.boot = X[index]
        delta.boot = delta[index]
        z1.boot = z1[index]
        z2.boot = z2[index]
        z3.boot = z3[index]
        zanshi=try(estNRI_Weibull(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0),silent=TRUE)
        if(class(zanshi)=="try-error") next
        boot.NRI = rbind(boot.NRI,t(zanshi))
      }
      
      zheng.NRI = NULL
      for(i in 1:nboot)
      {
        zanshi=try(estNRI_Weibull(X,delta,z1, z2,z3, t0, TRUE),silent=TRUE)
        if(class(zanshi)=="try-error") next
        zheng.NRI = rbind(zheng.NRI,t(zanshi))
      }
      
      
      
      boot.NRI = na.omit(boot.NRI)
      zheng.NRI = na.omit(zheng.NRI)
      if(nrow(boot.NRI)<=nboot/2 | nrow(zheng.NRI)<=nboot/2) return(NULL)
      
      sd.NRI = t(apply(boot.NRI,2,sd))
      
      NRI_zheng=0.632*colMeans(zheng.NRI)+(1-0.632)*NRI
      bias=NRI-NRI_zheng
      covs=(NRI+qnorm(0.025)*sd.NRI-bias < c(0.1140708,0.09986274,0.1257358))&(NRI+qnorm(0.975)*sd.NRI-bias > c(0.1140708,0.09986274,0.1257358))
      
      res=rbind(res,rbind(NRI,sd.NRI,NRI_zheng,covs))
    }
    
  }
  return(res)
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
NRIWBUsingCox400_Alternative_zheng <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(NRIWBUsingCox400_Alternative_zheng,file = 'NRIWBUsingCox400_Alternative_zheng.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF