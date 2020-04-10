load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueNRI_Fine_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueIDI_Fine_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueNRI_WB_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueIDI_WB_CV.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueNRI_Fine_testing.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueIDI_Fine_testing.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueNRI_WB_testing.rda")
load("E:\\Research\\NRIIDI revisions\\Cross Validation\\TrueIDI_WB_testing.rda")

N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueNRI_WB_CV[[i]][1,]))
  N2=rbind(N2,t(TrueNRI_WB_CV[[i]][2,]))
  N3=rbind(N3,t(TrueNRI_WB_CV[[i]][3,]))
  sum1=rbind(sum1,t(TrueNRI_WB_CV[[i]][4,]))
  sum2=rbind(sum2,t(TrueNRI_WB_CV[[i]][5,]))
  sum3=rbind(sum3,t(TrueNRI_WB_CV[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum))/3 ,
  sum(apply(sum2,2,sum)/apply(N2,2,sum))/3 ,
  sum(apply(sum3,2,sum)/apply(N3,2,sum))/3 )

N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueNRI_WB_testing[[i]][1,]))
  N2=rbind(N2,t(TrueNRI_WB_testing[[i]][2,]))
  N3=rbind(N3,t(TrueNRI_WB_testing[[i]][3,]))
  sum1=rbind(sum1,t(TrueNRI_WB_testing[[i]][4,]))
  sum2=rbind(sum2,t(TrueNRI_WB_testing[[i]][5,]))
  sum3=rbind(sum3,t(TrueNRI_WB_testing[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum))/3 ,
  sum(apply(sum2,2,sum)/apply(N2,2,sum))/3 ,
  sum(apply(sum3,2,sum)/apply(N3,2,sum))/3 )


N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueNRI_Fine_CV[[i]][1,]))
  N2=rbind(N2,t(TrueNRI_Fine_CV[[i]][2,]))
  N3=rbind(N3,t(TrueNRI_Fine_CV[[i]][3,]))
  sum1=rbind(sum1,t(TrueNRI_Fine_CV[[i]][4,]))
  sum2=rbind(sum2,t(TrueNRI_Fine_CV[[i]][5,]))
  sum3=rbind(sum3,t(TrueNRI_Fine_CV[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum))/3 ,
  sum(apply(sum2,2,sum)/apply(N2,2,sum))/3 ,
  sum(apply(sum3,2,sum)/apply(N3,2,sum))/3 )

N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueNRI_Fine_testing[[i]][1,]))
  N2=rbind(N2,t(TrueNRI_Fine_testing[[i]][2,]))
  N3=rbind(N3,t(TrueNRI_Fine_testing[[i]][3,]))
  sum1=rbind(sum1,t(TrueNRI_Fine_testing[[i]][4,]))
  sum2=rbind(sum2,t(TrueNRI_Fine_testing[[i]][5,]))
  sum3=rbind(sum3,t(TrueNRI_Fine_testing[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum))/3 ,
  sum(apply(sum2,2,sum)/apply(N2,2,sum))/3 ,
  sum(apply(sum3,2,sum)/apply(N3,2,sum))/3 )

N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueIDI_WB_CV[[i]][1,]))
  N2=rbind(N2,t(TrueIDI_WB_CV[[i]][2,]))
  N3=rbind(N3,t(TrueIDI_WB_CV[[i]][3,]))
  sum1=rbind(sum1,t(TrueIDI_WB_CV[[i]][4,]))
  sum2=rbind(sum2,t(TrueIDI_WB_CV[[i]][5,]))
  sum3=rbind(sum3,t(TrueIDI_WB_CV[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum)/(1-apply(N1,2,sum)/sum(N1)))/3,
  sum(apply(sum2,2,sum)/apply(N2,2,sum)/(1-apply(N2,2,sum)/sum(N2)))/3,
  sum(apply(sum3,2,sum)/apply(N3,2,sum)/(1-apply(N3,2,sum)/sum(N3)))/3)

N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueIDI_WB_testing[[i]][1,]))
  N2=rbind(N2,t(TrueIDI_WB_testing[[i]][2,]))
  N3=rbind(N3,t(TrueIDI_WB_testing[[i]][3,]))
  sum1=rbind(sum1,t(TrueIDI_WB_testing[[i]][4,]))
  sum2=rbind(sum2,t(TrueIDI_WB_testing[[i]][5,]))
  sum3=rbind(sum3,t(TrueIDI_WB_testing[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum)/(1-apply(N1,2,sum)/sum(N1)))/3,
  sum(apply(sum2,2,sum)/apply(N2,2,sum)/(1-apply(N2,2,sum)/sum(N2)))/3,
  sum(apply(sum3,2,sum)/apply(N3,2,sum)/(1-apply(N3,2,sum)/sum(N3)))/3)


N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueIDI_Fine_CV[[i]][1,]))
  N2=rbind(N2,t(TrueIDI_Fine_CV[[i]][2,]))
  N3=rbind(N3,t(TrueIDI_Fine_CV[[i]][3,]))
  sum1=rbind(sum1,t(TrueIDI_Fine_CV[[i]][4,]))
  sum2=rbind(sum2,t(TrueIDI_Fine_CV[[i]][5,]))
  sum3=rbind(sum3,t(TrueIDI_Fine_CV[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum)/(1-apply(N1,2,sum)/sum(N1)))/3,
  sum(apply(sum2,2,sum)/apply(N2,2,sum)/(1-apply(N2,2,sum)/sum(N2)))/3,
  sum(apply(sum3,2,sum)/apply(N3,2,sum)/(1-apply(N3,2,sum)/sum(N3)))/3)

N1=N2=N3=sum1=sum2=sum3=NULL
for(i in 1:1000){
  N1=rbind(N1,t(TrueIDI_Fine_testing[[i]][1,]))
  N2=rbind(N2,t(TrueIDI_Fine_testing[[i]][2,]))
  N3=rbind(N3,t(TrueIDI_Fine_testing[[i]][3,]))
  sum1=rbind(sum1,t(TrueIDI_Fine_testing[[i]][4,]))
  sum2=rbind(sum2,t(TrueIDI_Fine_testing[[i]][5,]))
  sum3=rbind(sum3,t(TrueIDI_Fine_testing[[i]][6,]))
}
c(sum(apply(sum1,2,sum)/apply(N1,2,sum)/(1-apply(N1,2,sum)/sum(N1)))/3,
  sum(apply(sum2,2,sum)/apply(N2,2,sum)/(1-apply(N2,2,sum)/sum(N2)))/3,
  sum(apply(sum3,2,sum)/apply(N3,2,sum)/(1-apply(N3,2,sum)/sum(N3)))/3)
