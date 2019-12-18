####Final the data to be analyzed
dt=Final
X=dt$t
delta=dt$indicator
ndim=nrow(dt)
library(survival)
library(cmprsk)
library(riskRegression)

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
pos2=c(sum(X<=8),sum(X<=10),sum(X<=12))
Cen2=censor.S1[pos2]


##Brier score is calculated from fomula

Final$ind1=ifelse(Final$indicator==1,1,0)
Final$ind2=ifelse(Final$indicator==2,1,0)
M1=CSC(Hist(t,indicator)~AGE+CESD+LEU3N+COHORT+LEU3N2,data=Final,cause=1)
M2=FGR(Hist(t,indicator)~AGE+CESD+LEU3N+COHORT,data=Final,cause=1)
R1=predictRisk(M1,Final[,c('AGE','CESD','LEU3N','COHORT','LEU3N2')],times = c(8,10,12),cause=1)
R2=predictRisk(M2,Final[,c('AGE','CESD','LEU3N','COHORT')],times = c(8,10,12),cause=1)
w1=ifelse((Final$t<=8)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>8,1,0)/Cen2[1]
mean((ifelse((Final$t<=8)&Final$ind1,1,0)-R1[,1])^2*w1)
w1=ifelse((Final$t<=10)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>10,1,0)/Cen2[1]
mean((ifelse((Final$t<=10)&Final$ind1,1,0)-R1[,2])^2*w1)
w1=ifelse((Final$t<=12)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>12,1,0)/Cen2[1]
mean((ifelse((Final$t<=12)&Final$ind1,1,0)-R1[,3])^2*w1)

w1=ifelse((Final$t<=8)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>8,1,0)/Cen2[1]
mean((ifelse((Final$t<=8)&Final$ind1,1,0)-R2[,1])^2*w1)
w1=ifelse((Final$t<=10)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>10,1,0)/Cen2[1]
mean((ifelse((Final$t<=10)&Final$ind1,1,0)-R2[,2])^2*w1)
w1=ifelse((Final$t<=12)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>12,1,0)/Cen2[1]
mean((ifelse((Final$t<=12)&Final$ind1,1,0)-R2[,3])^2*w1)



Final$co=ifelse(Final$COHORT==2,0,1)
Final$ind1=ifelse(Final$indicator==1,1,0)
Final$ind2=ifelse(Final$indicator==2,1,0)
M1=CSC(Hist(t,indicator)~AGE+CESD+LEU3N+COHORT+LEU3N2,data=Final,cause=2)
M2=FGR(Hist(t,indicator)~AGE+CESD+LEU3N+COHORT,data=Final,cause=2)
R1=predictRisk(M1,Final[,c('AGE','CESD','LEU3N','COHORT','LEU3N2')],times = c(8,10,12),cause=2)
R2=predictRisk(M2,Final[,c('AGE','CESD','LEU3N','COHORT')],times = c(8,10,12),cause=2)
w1=ifelse((Final$t<=8)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>8,1,0)/Cen2[1]
mean((ifelse((Final$t<=8)&Final$ind2,1,0)-R1[,1])^2*w1)
w1=ifelse((Final$t<=10)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>10,1,0)/Cen2[1]
mean((ifelse((Final$t<=10)&Final$ind2,1,0)-R1[,2])^2*w1)
w1=ifelse((Final$t<=12)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>12,1,0)/Cen2[1]
mean((ifelse((Final$t<=12)&Final$ind2,1,0)-R1[,3])^2*w1)

w1=ifelse((Final$t<=8)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>8,1,0)/Cen2[1]
mean((ifelse((Final$t<=8)&Final$ind2,1,0)-R2[,1])^2*w1)
w1=ifelse((Final$t<=10)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>10,1,0)/Cen2[1]
mean((ifelse((Final$t<=10)&Final$ind2,1,0)-R2[,2])^2*w1)
w1=ifelse((Final$t<=12)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>12,1,0)/Cen2[1]
mean((ifelse((Final$t<=12)&Final$ind2,1,0)-R2[,3])^2*w1)

##F31 & F32 are predicted probabilities from Gerds model implemented
load('E:Legacy_Gerds_probs.rda')
F31=res[[4]]
F32=res[[5]]
w1=ifelse((Final$t<=8)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>8,1,0)/Cen2[1]
mean((ifelse((Final$t<=8)&Final$ind1,1,0)-F31[,1])^2*w1)
w1=ifelse((Final$t<=10)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>10,1,0)/Cen2[1]
mean((ifelse((Final$t<=10)&Final$ind1,1,0)-F31[,2])^2*w1)
w1=ifelse((Final$t<=12)&Final$ind1,1,0)/est.S.censoring+ifelse(Final$ind1==0 & Final$t>12,1,0)/Cen2[1]
mean((ifelse((Final$t<=12)&Final$ind1,1,0)-F31[,3])^2*w1)

w1=ifelse((Final$t<=8)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>8,1,0)/Cen2[1]
mean((ifelse((Final$t<=8)&Final$ind2,1,0)-F32[,1])^2*w1)
w1=ifelse((Final$t<=10)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>10,1,0)/Cen2[1]
mean((ifelse((Final$t<=10)&Final$ind2,1,0)-F32[,2])^2*w1)
w1=ifelse((Final$t<=12)&Final$ind2,1,0)/est.S.censoring+ifelse(Final$ind2==0 & Final$t>12,1,0)/Cen2[1]
mean((ifelse((Final$t<=12)&Final$ind2,1,0)-F32[,3])^2*w1)


fit1=coxph(Surv(t,ind1)~AGE+CESD+LEU3N+COHORT+LEU3N2,data=Final,method='breslow')
coxsnellres=Final$ind1-resid(fit1,type="martingale")

## Then using NA method to estimate the cumulative hazard function for residuals;
fitres=survfit(coxph(Surv(coxsnellres,Final$ind1)~1,method='breslow'),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab='Cox-Snell Residuals', 
     ylab='Estimated Cumulative Hazard Function',
     main='Cox Regression with Event of Cognitive Impairment')
abline(0,1,col='red',lty=2)

fit1=coxph(Surv(t,ind2)~AGE+CESD+LEU3N+COHORT,data=Final,method='breslow')
coxsnellres=Final$ind2-resid(fit1,type="martingale")

## Then using NA method to estimate the cumulative hazard function for residuals;
fitres=survfit(coxph(Surv(coxsnellres,Final$ind2)~1,method='breslow'),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab='Cox-Snell Residuals', 
     ylab='Estimated Cumulative Hazard Function',
     main='Cox Regression with Event of Death')
abline(0,1,col='red',lty=2)

