load("E:\\Research\\NRIIDI revisions\\632bootstrap\\IDIFineUsingFine400_Alternative_zheng.rda")
load("E:\\Research\\NRIIDI revisions\\632bootstrap\\NRIFineUsingFine400_Alternative_zheng.rda")
load("E:\\Research\\NRIIDI revisions\\632bootstrap\\IDIWBUsingCox400_Alternative_zheng.rda")
load("E:\\Research\\NRIIDI revisions\\632bootstrap\\NRIWBUsingCox400_Alternative_zheng.rda")
load("E:\\Research\\NRIIDI revisions\\632bootstrap\\NRIIDIGerdsUsingGerds400_Alternative_zheng.rda")

#######################################################################################
### 30% NRI 400

NRI=NRI_zheng=cov_boot=NRI_boot_sd=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_Alternative_zheng[[i]][1,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_Alternative_zheng[[i]][2,])
  NRI_zheng=rbind(NRI_zheng,NRIWBUsingCox400_Alternative_zheng[[i]][3,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_Alternative_zheng[[i]][4,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(NRI_zheng)
apply(NRI_zheng,2,sd)
colMeans(cov_boot,na.rm=TRUE)

NRI=NRI_zheng=cov_boot=NRI_boot_sd=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingFine400_Alternative_zheng[[i]][1,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_Alternative_zheng[[i]][2,])
  NRI_zheng=rbind(NRI_zheng,NRIFineUsingFine400_Alternative_zheng[[i]][3,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_Alternative_zheng[[i]][4,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(NRI_zheng)
apply(NRI_zheng,2,sd)
colMeans(cov_boot,na.rm=TRUE)


#######################################################################################
### 50% NRI 400

NRI=NRI_zheng=cov_boot=NRI_boot_sd=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIWBUsingCox400_Alternative_zheng[[i]][5,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIWBUsingCox400_Alternative_zheng[[i]][6,])
  NRI_zheng=rbind(NRI_zheng,NRIWBUsingCox400_Alternative_zheng[[i]][7,])
  cov_boot=rbind(cov_boot,NRIWBUsingCox400_Alternative_zheng[[i]][8,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(NRI_zheng)
apply(NRI_zheng,2,sd)
colMeans(cov_boot,na.rm=TRUE)

NRI=NRI_zheng=cov_boot=NRI_boot_sd=NULL
for(i in 1:1000){
  NRI=rbind(NRI,NRIFineUsingFine400_Alternative_zheng[[i]][5,])
  NRI_boot_sd=rbind(NRI_boot_sd,NRIFineUsingFine400_Alternative_zheng[[i]][6,])
  NRI_zheng=rbind(NRI_zheng,NRIFineUsingFine400_Alternative_zheng[[i]][7,])
  cov_boot=rbind(cov_boot,NRIFineUsingFine400_Alternative_zheng[[i]][8,])
}
colMeans(NRI)
apply(NRI,2,sd)
colMeans(NRI_boot_sd)
colMeans(NRI_zheng)
apply(NRI_zheng,2,sd)
colMeans(cov_boot,na.rm=TRUE)


