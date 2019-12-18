load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingCox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingFine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingCox200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingFine200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingCox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingFine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingCox200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingFine200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIIDIWBUsingGerds200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIIDIWBUsingGerds400_Alternative.rda")

#######################################################################################
### 30% NRI 400
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine400_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

#######################################################################################
### 30% IDI 400

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine400_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine400_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine400_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine400_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds400_Alternative[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds400_Alternative[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_Alternative[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_Alternative[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)




#######################################################################################
### 30% NRI 200
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox200_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine200_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

#######################################################################################
### 30% IDI 200

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox200_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox200_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox200_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox200_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine200_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine200_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine200_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine200_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds200_Alternative[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds200_Alternative[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_Alternative[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_Alternative[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)



#######################################################################################
### 50% NRI 400
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox400_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox400_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_Alternative[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine400_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine400_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine400_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine400_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine400_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine400_Alternative[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds400_Alternative[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds400_Alternative[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds400_Alternative[[i]][13,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds400_Alternative[[i]][14,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_Alternative[[i]][15,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_Alternative[[i]][16,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

#######################################################################################
### 50% IDI 400

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine400_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine400_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine400_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine400_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds400_Alternative[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds400_Alternative[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_Alternative[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_Alternative[[i]][20,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)





#######################################################################################
### 50% NRI 200
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox200_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox200_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox200_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox200_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox200_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox200_Alternative[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine200_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine200_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine200_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine200_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine200_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine200_Alternative[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds200_Alternative[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds200_Alternative[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds200_Alternative[[i]][13,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds200_Alternative[[i]][14,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_Alternative[[i]][15,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_Alternative[[i]][16,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

#######################################################################################
### 50% IDI 200

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox200_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox200_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox200_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox200_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine200_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine200_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine200_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine200_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds200_Alternative[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds200_Alternative[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_Alternative[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_Alternative[[i]][20,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)



