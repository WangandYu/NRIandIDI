#!/bin/bash
#SBATCH --job-name=NRIFF400
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
  
estNRI_Fine = function(X,delta, z1, z2,z3, t0, withsd)
{ndim = length(X)
npts = length(t0)

#############################################################################
### Fit the fine gray model without z3 (model 1) ###
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
### Fit the fine gray model using all covariates (model 2) ###
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

return(cbind(estNRI,estsd))

}

# number of simulation
nsim = 1
# sample size
ndim = 400
nboot=1000



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
    
    
    ests=try(estNRI_Fine(X,delta,z1, z2,z3, t0, TRUE),silent=TRUE)
    if(class(ests)=="try-error") return(NULL)
    NRI = t(ests[,1])
    NRI_sd = t(ests[,2])
    
    cov = as.numeric((NRI-qnorm(0.975)*NRI_sd/sqrt(ndim))<c(0.109168,0.1268454,0.1225807) & (NRI+qnorm(0.975)*NRI_sd/sqrt(ndim))>c(0.109168,0.1268454,0.1225807))
    
    
    
    boot.NRI = NULL
    for(i in 1:nboot)
    {
      index = sample(1:ndim, ndim, replace=TRUE)
      X.boot = X[index]
      delta.boot = delta[index]
      z1.boot = z1[index]
      z2.boot = z2[index]
      z3.boot = z3[index]
      zanshi=try(estNRI_Fine(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0, FALSE),silent=TRUE)
      if(class(zanshi)=="try-error") next
      boot.NRI = rbind(boot.NRI,t(zanshi))
    }
    jack.NRI = NULL
    for(i in 1:ndim)
    {
      X.boot = X[-i]
      delta.boot = delta[-i]
      z1.boot = z1[-i]
      z2.boot = z2[-i]
      z3.boot = z3[-i]
      zanshi=try(estNRI_Fine(X.boot,delta.boot,z1.boot, z2.boot,z3.boot, t0, FALSE),silent=TRUE)
      if(class(zanshi)=="try-error") next
      jack.NRI = rbind(jack.NRI,t(zanshi))
    }
    
    
    
    boot.NRI2 = na.omit(boot.NRI)
    sd.NRI = t(apply(boot.NRI2,2,sd))
    
    z0.hat=qnorm(colMeans(boot.NRI<matrix(rep(NRI,each=nrow(boot.NRI)),ncol=3)))
    temp=matrix(rep(colMeans(jack.NRI),each=nrow(jack.NRI)),ncol=3)-jack.NRI
    a.hat=colSums(temp^3)/6/colSums(temp^2)^(3/2)
    alpha1=pnorm(z0.hat+(z0.hat+qnorm(0.025))/(1-a.hat*(z0.hat+qnorm(0.025))))
    alpha2=pnorm(z0.hat+(z0.hat+qnorm(0.975))/(1-a.hat*(z0.hat+qnorm(0.975))))
    cov1=ifelse(quantile(boot.NRI[,1],alpha1[1])<0.109168 & quantile(boot.NRI[,1],alpha2[1])>0.109168,1,0)
    cov2=ifelse(quantile(boot.NRI[,2],alpha1[2])<0.1268454 & quantile(boot.NRI[,2],alpha2[2])>0.1268454,1,0)
    cov3=ifelse(quantile(boot.NRI[,3],alpha1[3])<0.1225807 & quantile(boot.NRI[,3],alpha2[3])>0.1225807,1,0)
    cov10=ifelse(quantile(boot.NRI[,1],0.025)<0.109168 & quantile(boot.NRI[,1],0.975)>0.109168,1,0)
    cov20=ifelse(quantile(boot.NRI[,2],0.025)<0.1268454 & quantile(boot.NRI[,2],0.975)>0.1268454,1,0)
    cov30=ifelse(quantile(boot.NRI[,3],0.025)<0.1225807 & quantile(boot.NRI[,3],0.975)>0.1225807,1,0)
    
    res=rbind(res,rbind(NRI,NRI_sd,cov,sd.NRI,c(cov1,cov2,cov3),c(cov10,cov20,cov30)))
  }
  
}
return(res)
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
NRIFineUsingFine400_Alternative <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(NRIFineUsingFine400_Alternative,file = 'NRIFineUsingFine400_Alternative.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF