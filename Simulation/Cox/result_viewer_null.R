load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingCox400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingFine400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingCox200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIWBUsingFine200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingCox400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingFine400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingCox200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\IDIWBUsingFine200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIIDIWBUsingGerds200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Cox\\NRIIDIWBUsingGerds400_NULL.rda")


#######################################################################################
### 30% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_NULL[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% NRI 400 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox400_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox400_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_NULL[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine400_NULL[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% NRI 400 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine400_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine400_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine400_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine400_NULL[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_NULL[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% NRI 400 after
NRI=NRI_IF_sd=cov_IF=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds400_NULL[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds400_NULL[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds400_NULL[[i]][13,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)


#######################################################################################
### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine400_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine400_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine400_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine400_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine400_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine400_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine400_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine400_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds400_NULL[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds400_NULL[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_NULL[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_NULL[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds400_NULL[[i]][14,])
}
colMeans(IDI)
apply(IDI,2,sd)


#######################################################################################
### 30% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox200_NULL[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% NRI 200 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox200_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox200_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox200_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox200_NULL[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine200_NULL[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% NRI 200 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine200_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine200_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine200_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine200_NULL[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_NULL[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% NRI 200 after
NRI=NRI_IF_sd=cov_IF=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds200_NULL[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds200_NULL[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds200_NULL[[i]][13,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)



#######################################################################################
### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox200_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox200_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox200_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox200_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox200_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox200_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox200_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox200_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine200_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine200_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine200_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine200_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine200_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine200_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine200_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine200_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds200_NULL[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds200_NULL[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_NULL[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_NULL[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds200_NULL[[i]][14,])
}
colMeans(IDI)
apply(IDI,2,sd)



#######################################################################################
### 50% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox400_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox400_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_NULL[[i]][18,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 400 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox400_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox400_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_NULL[[i]][24,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine400_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine400_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine400_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine400_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine400_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine400_NULL[[i]][18,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 400 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine400_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine400_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine400_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine400_NULL[[i]][24,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds400_NULL[[i]][15,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds400_NULL[[i]][16,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds400_NULL[[i]][17,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds400_NULL[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_NULL[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_NULL[[i]][20,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 400 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds400_NULL[[i]][25,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds400_NULL[[i]][26,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds400_NULL[[i]][27,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)


#######################################################################################
### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox400_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox400_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox400_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox400_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine400_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine400_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine400_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine400_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine400_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine400_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds400_NULL[[i]][21,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds400_NULL[[i]][24,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds400_NULL[[i]][28,])
}
colMeans(IDI)
apply(IDI,2,sd)



#######################################################################################
### 50% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox200_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox200_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox200_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox200_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox200_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox200_NULL[[i]][18,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 200 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox200_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingCox200_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIWBUsingCox200_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIWBUsingCox200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox200_NULL[[i]][24,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine200_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine200_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine200_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine200_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine200_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine200_NULL[[i]][18,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 200 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingFine200_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIWBUsingFine200_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIWBUsingFine200_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingFine200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIWBUsingFine200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIWBUsingFine200_NULL[[i]][24,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds200_NULL[[i]][15,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds200_NULL[[i]][16,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds200_NULL[[i]][17,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIWBUsingGerds200_NULL[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_NULL[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_NULL[[i]][20,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% NRI 200 after
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIWBUsingGerds200_NULL[[i]][25,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIWBUsingGerds200_NULL[[i]][26,])
  cov_IF=rbind(cov_IF,NRIIDIWBUsingGerds200_NULL[[i]][27,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)

#######################################################################################
### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox200_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox200_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
apply(IDI_boot_sd,2,function(v) mean(v[v<0.05]))
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingCox200_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingCox200_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIWBUsingCox200_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIWBUsingCox200_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine200_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine200_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIWBUsingFine200_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIWBUsingFine200_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIWBUsingFine200_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIWBUsingFine200_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds200_NULL[[i]][21,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIWBUsingGerds200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIIDIWBUsingGerds200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIIDIWBUsingGerds200_NULL[[i]][24,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIWBUsingGerds200_NULL[[i]][28,])
}
colMeans(IDI)
apply(IDI,2,sd)





