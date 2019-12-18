#!/bin/bash
#SBATCH --job-name=IDIFC200
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
  
  estIDI_Weibull = function(X,delta, z1, z2,z3, t0, withsd)
  {ndim = length(X)
  npts = length(t0)
  
  #############################################################################
  ### Fit the cox proportional hazard regression model without z3 (model 1) ###
  #############################################################################
  
  # calcualte the cumulative baseline hazard function for cause 1 #
  data1 = data.frame(cbind(X, event=ifelse(delta==1,1,0), z1=z1,z2=z2))
  names(data1) = c("X","event","z1","z2")
  cox1 = coxph(Surv(X, event) ~ z1+z2 , data1)	
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
  ndim = 200
  nboot= 1000
  
  
  ### true parameters ###
  beta11=0.2
  beta12=-0.5
  beta13=1
  p=0.65
  
  t0 =  c(20,21,22)		
  
  LC=c(9,7,10,3,7)
  RC=c(38,30,37.5,25,29.5)
  per=c(10,20,30,40,50)
  res=NULL
  for(ww in c(3,5)){
    
    for ( m in 1:nsim)	# number of simulations
    {
      ### generate the covariates ###
      
      z1=rnorm(ndim)
      z1=pmax(pmin(z1,3.5),-3.5)
      z2=ifelse(runif(ndim)<0.7,1,0)
      z3=rnorm(ndim)
      z3=pmax(pmin(z3,3.5),-3.5)
      
      F1=1-(1-p)^(exp(z1*beta11+z2*beta12+z3*beta13))
      
      W=runif(ndim)
      epsilon=ifelse(W<F1,1,2)
      TT=rep(0,ndim)
      TT[W<F1]=-log(1-(1-(1-W[W<F1])^(1/exp(z1[W<F1]*beta11+z2[W<F1]*beta12+z3[W<F1]*beta13)))/p)
      zz=exp(z1[W>=F1]*beta11/10+z2[W>=F1]*beta12/10+z3[W>=F1]*beta13/10)
      TT[W>=F1]=sapply(zz,function(o) return(rexp(1,o)))
      TT=TT^0.2
      TT=pmin(TT*20,100)
      ## uniformly generate censoring time ##
      
      C = runif(ndim, LC[ww], RC[ww])		
      
      X = as.numeric(pmin(TT,C))
      delta = as.numeric(TT<=C)*epsilon
      
      ests=try(estIDI_Weibull(X,delta,z1, z2,z3, t0, FALSE),silent=TRUE)
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
        zanshi=try(estIDI_Weibull(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0, FALSE),silent=TRUE)
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
        zanshi=try(estIDI_Weibull(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0, FALSE),silent=TRUE)
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
      cov1=ifelse(quantile(boot.IDI[,1],alpha1[1])<0.1513984 & quantile(boot.IDI[,1],alpha2[1])>0.1513984,1,0)
      cov2=ifelse(quantile(boot.IDI[,2],alpha1[2])<0.1584317 & quantile(boot.IDI[,2],alpha2[2])>0.1584317,1,0)
      cov3=ifelse(quantile(boot.IDI[,3],alpha1[3])<0.1649372 & quantile(boot.IDI[,3],alpha2[3])>0.1649372,1,0)
      cov10=ifelse(quantile(boot.IDI[,1],0.025)<0.1513984 & quantile(boot.IDI[,1],0.975)>0.1513984,1,0)
      cov20=ifelse(quantile(boot.IDI[,2],0.025)<0.1584317 & quantile(boot.IDI[,2],0.975)>0.1584317,1,0)
      cov30=ifelse(quantile(boot.IDI[,3],0.025)<0.1649372 & quantile(boot.IDI[,3],0.975)>0.1649372,1,0)
      
      res=rbind(res,rbind(IDI,sd.IDI,c(cov1,cov2,cov3),c(cov10,cov20,cov30)))
    }}
  return(res)
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
IDIFineUsingCox200_Alternative <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(IDIFineUsingCox200_Alternative,file = 'IDIFineUsingCox200_Alternative.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF
