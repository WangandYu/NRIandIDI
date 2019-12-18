load("E:\\Research\\NRIIDI revisions\\Shi2014\\estShi_Cox400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Shi2014\\estShi_Fine400_Alternative.rda")
load("E:\\Research\\NRIIDI revisions\\Shi2014\\estShi_Gerds400_Alternative.rda")

Shi=sd_Shi=NULL
for(i in 1:1000){
  Shi=rbind(Shi,estShi_Cox400_Alternative[[i]][1,])
}
colMeans(Shi)
apply(Shi,2,sd)


Shi=sd_Shi=NULL
for(i in 1:1000){
  Shi=rbind(Shi,estShi_Fine400_Alternative[[i]][1,])
}
colMeans(Shi,na.rm = TRUE)
apply(Shi,2,function(v) sd(na.omit(v)))




Shi=sd_Shi=NULL
for(i in 1:1000){
  Shi=rbind(Shi,estShi_Gerds400_Alternative[[i]][1,])
}
colMeans(Shi)
apply(Shi,2,sd)

Shi=sd_Shi=NULL
for(i in 1:1000){
  Shi=rbind(Shi,estShi_Cox400_Alternative[[i]][2,])
}
colMeans(Shi)
apply(Shi,2,sd)


Shi=sd_Shi=NULL
for(i in 1:1000){
  Shi=rbind(Shi,estShi_Fine400_Alternative[[i]][2,])
}
colMeans(Shi,na.rm = TRUE)
apply(Shi,2,function(v) sd(na.omit(v)))




Shi=sd_Shi=NULL
for(i in 1:1000){
  Shi=rbind(Shi,estShi_Gerds400_Alternative[[i]][2,])
}
colMeans(Shi)
apply(Shi,2,sd)