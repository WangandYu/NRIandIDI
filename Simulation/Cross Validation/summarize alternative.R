load("E:\\Research\\NRIIDI revisions\\Cross Validation\\IDIFineUsingFine400_Alternative_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\NRIFineUsingFine400_Alternative_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\IDIWBUsingCox400_Alternative_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\NRIWBUsingCox400_Alternative_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\NRIIDIGerdsUsingGerds400_Alternative_CV.rda")

#######################################################################################
### 30% NRI 400

NRI=NRI_boot_sd=cov_bca=cov_boot=cov_norm=cov_basic=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_Alternative_CV[[i]][1,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_Alternative_CV[[i]][2,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_Alternative_CV[[i]][3,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_Alternative_CV[[i]][4,])
  cov_norm=rbind(cov_norm,NRIWBUsingCox400_Alternative_CV[[i]][5,])
  cov_basic=rbind(cov_basic,NRIWBUsingCox400_Alternative_CV[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot,na.rm=TRUE)
colMeans(cov_norm,na.rm=TRUE)
colMeans(cov_basic,na.rm=TRUE)

NRI=NRI_boot_sd=cov_bca=cov_boot=cov_norm=cov_basic=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingFine400_Alternative_CV[[i]][1,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_Alternative_CV[[i]][2,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_Alternative_CV[[i]][3,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_Alternative_CV[[i]][4,])
  cov_norm=rbind(cov_norm,NRIFineUsingFine400_Alternative_CV[[i]][5,])
  cov_basic=rbind(cov_basic,NRIFineUsingFine400_Alternative_CV[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot,na.rm=TRUE)
colMeans(cov_norm,na.rm=TRUE)
colMeans(cov_basic,na.rm=TRUE)

NRI=NRI_boot_sd=cov_bca=cov_boot=cov_norm=cov_basic=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][1,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][2,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][3,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][4,])
  cov_norm=rbind(cov_norm,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][5,])
  cov_basic=rbind(cov_basic,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot,na.rm=TRUE)
colMeans(cov_norm,na.rm=TRUE)
colMeans(cov_basic,na.rm=TRUE)


#######################################################################################
### 50% NRI 400

NRI=NRI_boot_sd=cov_bca=cov_boot=cov_norm=cov_basic=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_Alternative_CV[[i]][7,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_Alternative_CV[[i]][8,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_Alternative_CV[[i]][9,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_Alternative_CV[[i]][10,])
  cov_norm=rbind(cov_norm,NRIWBUsingCox400_Alternative_CV[[i]][11,])
  cov_basic=rbind(cov_basic,NRIWBUsingCox400_Alternative_CV[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot,na.rm=TRUE)
colMeans(cov_norm,na.rm=TRUE)
colMeans(cov_basic,na.rm=TRUE)

NRI=NRI_boot_sd=cov_bca=cov_boot=cov_norm=cov_basic=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingFine400_Alternative_CV[[i]][7,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_Alternative_CV[[i]][8,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_Alternative_CV[[i]][9,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_Alternative_CV[[i]][10,])
  cov_norm=rbind(cov_norm,NRIFineUsingFine400_Alternative_CV[[i]][11,])
  cov_basic=rbind(cov_basic,NRIFineUsingFine400_Alternative_CV[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot,na.rm=TRUE)
colMeans(cov_norm,na.rm=TRUE)
colMeans(cov_basic,na.rm=TRUE)



NRI=NRI_boot_sd=cov_bca=cov_boot=cov_norm=cov_basic=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][11,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][12,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][13,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][14,])
  cov_norm=rbind(cov_norm,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][15,])
  cov_basic=rbind(cov_basic,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][16,])
}

colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot,na.rm=TRUE)
colMeans(cov_norm,na.rm=TRUE)
colMeans(cov_basic,na.rm=TRUE)



#######################################################################################
### 30% IDI 400

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_Alternative_CV[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_Alternative_CV[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_Alternative_CV[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_Alternative_CV[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)



IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_Alternative_CV[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_Alternative_CV[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_Alternative_CV[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_Alternative_CV[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)



IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][10,])
}

colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


#######################################################################################
### 50% IDI 400

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_Alternative_CV[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_Alternative_CV[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_Alternative_CV[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_Alternative_CV[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_Alternative_CV[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_Alternative_CV[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_Alternative_CV[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_Alternative_CV[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative_CV[[i]][20,])
}

colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

