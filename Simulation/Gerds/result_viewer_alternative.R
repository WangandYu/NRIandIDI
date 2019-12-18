load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingCox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingFine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingCox200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIGerdsUsingFine200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingCox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingFine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingCox200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\IDIGerdsUsingFine200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIIDIGerdsUsingGerds200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Gerds\\NRIIDIGerdsUsingGerds400_Alternative.rda")

#######################################################################################
### 30% NRI 400
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingCox400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox400_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(400)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingFine400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine400_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(400)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative[[i]][6,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(400)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

#######################################################################################
### 30% IDI 400

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox400_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox400_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox400_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox400_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine400_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine400_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine400_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine400_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_Alternative[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative[[i]][10,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox200_Alternative[[i]][6,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine200_Alternative[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_Alternative[[i]][6,])
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
  IDI=rbind(IDI,IDIGerdsUsingCox200_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox200_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox200_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox200_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd[IDI_boot_sd[,1]<.1,])
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine200_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine200_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine200_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine200_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds200_Alternative[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds200_Alternative[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_Alternative[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_Alternative[[i]][10,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox400_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox400_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox400_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox400_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox400_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox400_Alternative[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(400)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIGerdsUsingFine400_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine400_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine400_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine400_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine400_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine400_Alternative[[i]][12,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(400)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds400_Alternative[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds400_Alternative[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds400_Alternative[[i]][13,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative[[i]][14,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative[[i]][15,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative[[i]][16,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_IF_sd)/sqrt(400)
colMeans(cov_IF)
colMeans(NRI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)

#######################################################################################
### 50% IDI 400

IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingCox400_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox400_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox400_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox400_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine400_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine400_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine400_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine400_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds400_Alternative[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds400_Alternative[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds400_Alternative[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds400_Alternative[[i]][20,])
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
  NRI=rbind(NRI,NRIGerdsUsingCox200_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingCox200_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingCox200_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingCox200_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingCox200_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingCox200_Alternative[[i]][12,])
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
  NRI=rbind(NRI,NRIGerdsUsingFine200_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIGerdsUsingFine200_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIGerdsUsingFine200_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIGerdsUsingFine200_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIGerdsUsingFine200_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIGerdsUsingFine200_Alternative[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIGerdsUsingGerds200_Alternative[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIGerdsUsingGerds200_Alternative[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIGerdsUsingGerds200_Alternative[[i]][13,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIGerdsUsingGerds200_Alternative[[i]][14,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_Alternative[[i]][15,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_Alternative[[i]][16,])
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
  IDI=rbind(IDI,IDIGerdsUsingCox200_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingCox200_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingCox200_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingCox200_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIGerdsUsingFine200_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIGerdsUsingFine200_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIGerdsUsingFine200_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIGerdsUsingFine200_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIGerdsUsingGerds200_Alternative[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIGerdsUsingGerds200_Alternative[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIGerdsUsingGerds200_Alternative[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIGerdsUsingGerds200_Alternative[[i]][20,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)



