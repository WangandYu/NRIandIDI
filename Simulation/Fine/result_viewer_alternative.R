load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingCox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingFine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingCox200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIFineUsingFine200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingCox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingFine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingCox200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\IDIFineUsingFine200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIIDIFineUsingGerds200_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Fine\\NRIIDIFineUsingGerds400_Alternative.rda")

#######################################################################################
### 30% NRI 400
NRI=NRI_IF_sd=cov_IF=NRI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingCox400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox400_Alternative[[i]][6,])
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
  NRI=rbind(NRI,NRIFineUsingFine400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_Alternative[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds400_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds400_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds400_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds400_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_Alternative[[i]][6,])
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
  IDI=rbind(IDI,IDIFineUsingCox400_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox400_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox400_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox400_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds400_Alternative[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds400_Alternative[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_Alternative[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_Alternative[[i]][10,])
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
  NRI=rbind(NRI,NRIFineUsingCox200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox200_Alternative[[i]][6,])
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
  NRI=rbind(NRI,NRIFineUsingFine200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine200_Alternative[[i]][6,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds200_Alternative[[i]][1,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds200_Alternative[[i]][2,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds200_Alternative[[i]][3,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds200_Alternative[[i]][4,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_Alternative[[i]][5,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_Alternative[[i]][6,])
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
  IDI=rbind(IDI,IDIFineUsingCox200_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox200_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox200_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox200_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine200_Alternative[[i]][1,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine200_Alternative[[i]][2,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine200_Alternative[[i]][3,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine200_Alternative[[i]][4,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds200_Alternative[[i]][7,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds200_Alternative[[i]][8,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_Alternative[[i]][9,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_Alternative[[i]][10,])
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
  NRI=rbind(NRI,NRIFineUsingCox400_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox400_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox400_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox400_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox400_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox400_Alternative[[i]][12,])
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
  NRI=rbind(NRI,NRIFineUsingFine400_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine400_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine400_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine400_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_Alternative[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds400_Alternative[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds400_Alternative[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds400_Alternative[[i]][13,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds400_Alternative[[i]][14,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_Alternative[[i]][15,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_Alternative[[i]][16,])
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
  IDI=rbind(IDI,IDIFineUsingCox400_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox400_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox400_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox400_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine400_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine400_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine400_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine400_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd)
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds400_Alternative[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds400_Alternative[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds400_Alternative[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds400_Alternative[[i]][20,])
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
  NRI=rbind(NRI,NRIFineUsingCox200_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingCox200_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingCox200_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingCox200_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingCox200_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingCox200_Alternative[[i]][12,])
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
  NRI=rbind(NRI,NRIFineUsingFine200_Alternative[[i]][7,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIFineUsingFine200_Alternative[[i]][8,])
  cov_IF=rbind(cov_IF,NRIFineUsingFine200_Alternative[[i]][9,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine200_Alternative[[i]][10,])
  cov_bca=rbind(cov_bca,NRIFineUsingFine200_Alternative[[i]][11,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine200_Alternative[[i]][12,])
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
  NRI=rbind(NRI,NRIIDIFineUsingGerds200_Alternative[[i]][11,])
  NRI_IF_sd=rbind(NRI_IF_sd,NRIIDIFineUsingGerds200_Alternative[[i]][12,])
  cov_IF=rbind(cov_IF,NRIIDIFineUsingGerds200_Alternative[[i]][13,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIIDIFineUsingGerds200_Alternative[[i]][14,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_Alternative[[i]][15,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_Alternative[[i]][16,])
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
  IDI=rbind(IDI,IDIFineUsingCox200_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingCox200_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingCox200_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingCox200_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
colMeans(IDI_boot_sd[IDI_boot_sd[,3]<0.1,])
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,IDIFineUsingFine200_Alternative[[i]][5,])
  IDI_boot_sd=rbind(IDI_boot_sd,IDIFineUsingFine200_Alternative[[i]][6,])
  cov_bca=rbind(cov_bca,IDIFineUsingFine200_Alternative[[i]][7,])
  cov_boot=rbind(cov_boot,IDIFineUsingFine200_Alternative[[i]][8,])
}
colMeans(IDI)
apply(IDI,2,sd)
apply(IDI_boot_sd,2,function(v) mean(v[v<0.1]))
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)


IDI=IDI_boot_sd=cov_bca=cov_boot=NULL
for(i in 1:1000){
  IDI=rbind(IDI,NRIIDIFineUsingGerds200_Alternative[[i]][17,])
  IDI_boot_sd=rbind(IDI_boot_sd,NRIIDIFineUsingGerds200_Alternative[[i]][18,])
  cov_bca=rbind(cov_bca,NRIIDIFineUsingGerds200_Alternative[[i]][19,])
  cov_boot=rbind(cov_boot,NRIIDIFineUsingGerds200_Alternative[[i]][20,])
}
colMeans(IDI)
apply(IDI,2,sd)
apply(IDI_boot_sd,2,function(v) mean(v[v<0.1]))
colMeans(cov_bca,na.rm=TRUE)
colMeans(cov_boot)



