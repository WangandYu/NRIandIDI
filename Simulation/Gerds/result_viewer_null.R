load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingCox400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingFine400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingCox200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingFine200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingCox400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingFine400_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingCox200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingFine200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIIDIGerdsUsingGerds200_NULL.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIIDIGerdsUsingGerds400_NULL.rda")


#######################################################################################
### 30% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingCox400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox400_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox400_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox400_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox400_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox400_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine400_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine400_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine400_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine400_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine400_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds400_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_NULL[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds400_NULL[[i]][13,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)


#######################################################################################
### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox400_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox400_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox400_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox400_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox400_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox400_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox400_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox400_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine400_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine400_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine400_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine400_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine400_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine400_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine400_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine400_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_NULL[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_NULL[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_NULL[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 400 after

IDI=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_NULL[[i]][14,])
}
colMeans(IDI)
apply(IDI,2,sd)


#######################################################################################
### 30% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingCox200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox200_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox200_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox200_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox200_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox200_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine200_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine200_NULL[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine200_NULL[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine200_NULL[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine200_NULL[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds200_NULL[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds200_NULL[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_NULL[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_NULL[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds200_NULL[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds200_NULL[[i]][13,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)



#######################################################################################
### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox200_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox200_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox200_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox200_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox200_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox200_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox200_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox200_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine200_NULL[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine200_NULL[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine200_NULL[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine200_NULL[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine200_NULL[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine200_NULL[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine200_NULL[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine200_NULL[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 30% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds200_NULL[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_NULL[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_NULL[[i]][10,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 30% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds200_NULL[[i]][14,])
}
colMeans(IDI)
apply(IDI,2,sd)



#######################################################################################
### 50% NRI 400 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingCox400_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox400_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox400_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox400_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox400_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox400_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox400_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox400_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox400_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox400_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine400_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine400_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine400_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine400_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine400_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine400_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine400_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine400_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine400_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine400_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_NULL[[i]][15,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][16,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds400_NULL[[i]][17,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_NULL[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_NULL[[i]][20,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_NULL[[i]][25,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][26,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds400_NULL[[i]][27,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/20
colMeans(cov_IF)


#######################################################################################
### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox400_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox400_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox400_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox400_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox400_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox400_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine400_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine400_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine400_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine400_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine400_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine400_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine400_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine400_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 400 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_NULL[[i]][21,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds400_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_NULL[[i]][24,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 400 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_NULL[[i]][28,])
}
colMeans(IDI)
apply(IDI,2,sd)



#######################################################################################
### 50% NRI 200 before
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingCox200_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox200_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox200_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox200_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox200_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox200_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox200_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox200_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox200_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox200_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine200_NULL[[i]][13,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine200_NULL[[i]][14,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine200_NULL[[i]][15,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine200_NULL[[i]][16,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine200_NULL[[i]][17,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine200_NULL[[i]][18,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine200_NULL[[i]][19,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine200_NULL[[i]][20,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine200_NULL[[i]][21,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine200_NULL[[i]][24,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds200_NULL[[i]][15,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][16,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds200_NULL[[i]][17,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_NULL[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_NULL[[i]][20,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds200_NULL[[i]][25,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][26,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds200_NULL[[i]][27,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(200)
colMeans(cov_IF)

#######################################################################################
### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox200_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox200_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
apply(IDI_boot_sd,2,function(v) mean(v[v<0.05]))
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox200_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox200_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox200_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox200_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine200_NULL[[i]][9,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine200_NULL[[i]][10,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine200_NULL[[i]][11,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine200_NULL[[i]][12,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine200_NULL[[i]][13,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine200_NULL[[i]][14,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine200_NULL[[i]][15,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine200_NULL[[i]][16,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


### 50% IDI 200 before

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds200_NULL[[i]][21,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds200_NULL[[i]][22,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_NULL[[i]][23,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_NULL[[i]][24,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

### 50% IDI 200 after

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds200_NULL[[i]][28,])
}
colMeans(IDI)
apply(IDI,2,sd)





