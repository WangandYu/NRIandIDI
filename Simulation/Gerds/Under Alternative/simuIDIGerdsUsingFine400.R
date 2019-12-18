#!/bin/bash
#SBATCH --job-name=IDIGF400
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
  
  estIDI_Fine = function(X,delta, z1, z2,z3, t0, withsd)
  {ndim = length(X)
  npts = length(t0)
  
  #############################################################################
  ### Fit the Fine Gray model without z3 (model 1) ###
  #############################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  cov = cbind(z1,z2)
  crr1 = crr(X,delta,cov)
  pred1 = predict(crr1,cov)
  # summary(crr1)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred1[,1]<t0[1])),max(which(pred1[,1]<t0[2])),max(which(pred1[,1]<t0[3])))
  
  # estimate the overall survival function #
  p1.hat.m1 = t(pred1[t0.index,-1])
  
  crr2 = crr(X,delta,cov,failcode=2)
  pred2 = predict(crr2,cov)
  # summary(crr1)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred2[,1]<t0[1])),max(which(pred2[,1]<t0[2])),max(which(pred2[,1]<t0[3])))
  
  # estimate the overall survival function #
  p2.hat.m1 = t(pred2[t0.index,-1])
  
  p3.hat.m1 = 1-p1.hat.m1-p2.hat.m1
  
  
  
  #######################################################################################
  ### Fit the Fine Gray model using all covariates (model 2) ###
  #######################################################################################
  # calcualte the cumulative baseline hazard function for cause 1 #
  cov = cbind(z1,z2,z3)
  crr1 = crr(X,delta,cov)
  pred1 = predict(crr1,cov)
  # summary(crr1)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred1[,1]<t0[1])),max(which(pred1[,1]<t0[2])),max(which(pred1[,1]<t0[3])))
  
  # estimate the overall survival function #
  p1.hat.m2 = t(pred1[t0.index,-1])
  
  crr2 = crr(X,delta,cov,failcode=2)
  pred2 = predict(crr2,cov)
  # summary(crr1)
  
  # indecies for the predified t0 #
  t0.index= c(max(which(pred2[,1]<t0[1])),max(which(pred2[,1]<t0[2])),max(which(pred2[,1]<t0[3])))
  
  # estimate the overall survival function #
  p2.hat.m2 = t(pred2[t0.index,-1])
  
  p3.hat.m2 = 1-p1.hat.m2-p2.hat.m2
  
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
  
  
  return(c(estIDI.t1, estIDI.t2, estIDI.t3))
  
  }
  # number of simulation
  nsim = 1
  # sample size
  ndim = 400
  nboot= 1000
  
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
      
      ests=try(estIDI_Fine(X,delta,z1, z2,z3, t0, FALSE),silent=TRUE)
      if(class(ests)=="try-error") return(NULL)
      IDI = t(ests)
      
      boot.IDI = NULL
      for(i in 1:nboot)
      {
        index = sample(1:ndim, ndim, replace=TRUE)
        X.boot = X[index]
        delta.boot = delta[index]
        z1.boot = z1[index]
        z2.boot = z2[index]
        z3.boot = z3[index]
        zanshi=try(estIDI_Fine(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0, FALSE),silent=TRUE)
        if(class(zanshi)=="try-error") next
        boot.IDI = rbind(boot.IDI,t(zanshi))
      }
      jack.IDI = NULL
      for(i in 1:ndim)
      {
        X.boot = X[-i]
        delta.boot = delta[-i]
        z1.boot = z1[-i]
        z2.boot = z2[-i]
        z3.boot = z3[-i]
        zanshi=try(estIDI_Fine(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0, FALSE),silent=TRUE)
        if(class(zanshi)=="try-error") next
        jack.IDI = rbind(jack.IDI,t(zanshi))
      }
      
      
      
      boot.IDI = na.omit(boot.IDI)
      jack.IDI = na.omit(jack.IDI)
      if(nrow(boot.IDI)<=nboot/2 | nrow(jack.IDI)<=ndim/2) return(NULL)
      sd.IDI = t(apply(boot.IDI,2,sd))
      
      z0.hat=qnorm(colMeans(boot.IDI<matrix(rep(IDI,each=nrow(boot.IDI)),ncol=3)))
      temp=matrix(rep(colMeans(jack.IDI),each=nrow(jack.IDI)),ncol=3)-jack.IDI
      a.hat=colSums(temp^3)/6/colSums(temp^2)^(3/2)
      alpha1=pnorm(z0.hat+(z0.hat+qnorm(0.025))/(1-a.hat*(z0.hat+qnorm(0.025))))
      alpha2=pnorm(z0.hat+(z0.hat+qnorm(0.975))/(1-a.hat*(z0.hat+qnorm(0.975))))
      cov1=ifelse(quantile(boot.IDI[,1],alpha1[1])<0.2875484 & quantile(boot.IDI[,1],alpha2[1])>0.2875484,1,0)
      cov2=ifelse(quantile(boot.IDI[,2],alpha1[2])<0.2713428 & quantile(boot.IDI[,2],alpha2[2])>0.2713428,1,0)
      cov3=ifelse(quantile(boot.IDI[,3],alpha1[3])<0.2509805 & quantile(boot.IDI[,3],alpha2[3])>0.2509805,1,0)
      cov10=ifelse(quantile(boot.IDI[,1],0.025)<0.2875484 & quantile(boot.IDI[,1],0.975)>0.2875484,1,0)
      cov20=ifelse(quantile(boot.IDI[,2],0.025)<0.2713428 & quantile(boot.IDI[,2],0.975)>0.2713428,1,0)
      cov30=ifelse(quantile(boot.IDI[,3],0.025)<0.2509805 & quantile(boot.IDI[,3],0.975)>0.2509805,1,0)
      
      res=rbind(res,rbind(IDI,sd.IDI,c(cov1,cov2,cov3),c(cov10,cov20,cov30)))
    }}
return(res)
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
IDIGerdsUsingFine400_Alternative <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(IDIGerdsUsingFine400_Alternative,file = 'IDIGerdsUsingFine400_Alternative.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF