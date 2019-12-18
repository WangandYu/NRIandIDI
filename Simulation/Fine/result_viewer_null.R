load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingCox400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingFine400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingCox200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingFine200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingCox400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingFine400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingCox200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingFine200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIIDIFineUsingGerds200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIIDIFineUsingGerds400_NULL.rda")


#######################################################################################
### 30% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingCox400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox400_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIFineUsingCox400_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox400_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox400_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox400_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIFineUsingFine400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIFineUsingFine400_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine400_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine400_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds400_NULL[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds400_NULL[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds400_NULL[[i]][13,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)


#######################################################################################
### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox400_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox400_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox400_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox400_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox400_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox400_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox400_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox400_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds400_NULL[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds400_NULL[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_NULL[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_NULL[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds400_NULL[[i]][14,])
}
colMeans(IDI)
apply(IDI,2,sd)


#######################################################################################
### 30% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingCox200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox200_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIFineUsingCox200_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox200_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox200_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox200_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIFineUsingFine200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine200_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIFineUsingFine200_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine200_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine200_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine200_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds200_NULL[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds200_NULL[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds200_NULL[[i]][13,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)



#######################################################################################
### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox200_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox200_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox200_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox200_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox200_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox200_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox200_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox200_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine200_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine200_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine200_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine200_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine200_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine200_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine200_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine200_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds200_NULL[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds200_NULL[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_NULL[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_NULL[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds200_NULL[[i]][14,])
}
colMeans(IDI)
apply(IDI,2,sd)



#######################################################################################
### 50% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingCox400_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox400_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox400_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox400_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox400_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox400_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIFineUsingCox400_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox400_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox400_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox400_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIFineUsingFine400_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine400_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine400_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIFineUsingFine400_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine400_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine400_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds400_NULL[[i]][15,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds400_NULL[[i]][16,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds400_NULL[[i]][17,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds400_NULL[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_NULL[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_NULL[[i]][20,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds400_NULL[[i]][25,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds400_NULL[[i]][26,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds400_NULL[[i]][27,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)


#######################################################################################
### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox400_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox400_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox400_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox400_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox400_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox400_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds400_NULL[[i]][21,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_NULL[[i]][24,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds400_NULL[[i]][28,])
}
colMeans(IDI)
apply(IDI,2,sd)



#######################################################################################
### 50% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingCox200_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox200_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox200_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox200_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox200_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox200_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIFineUsingCox200_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox200_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox200_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox200_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIFineUsingFine200_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine200_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine200_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine200_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine200_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine200_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIFineUsingFine200_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine200_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine200_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine200_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds200_NULL[[i]][15,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds200_NULL[[i]][16,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds200_NULL[[i]][17,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds200_NULL[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_NULL[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_NULL[[i]][20,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds200_NULL[[i]][25,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds200_NULL[[i]][26,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds200_NULL[[i]][27,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)

#######################################################################################
### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox200_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox200_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
apply(IDI_boot_sd,2,function(v) mean(v[v<0.05]))
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingCox200_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox200_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox200_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox200_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine200_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine200_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine200_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine200_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine200_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine200_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds200_NULL[[i]][21,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_NULL[[i]][24,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds200_NULL[[i]][28,])
}
colMeans(IDI)
apply(IDI,2,sd)





