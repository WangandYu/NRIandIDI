library(demography)
library(latex2exp)
par(mfrow=c(3,2))
x=(1:200000)/10000
y0=1-exp(-exp((log(x)-2.5)/0.2))
y0=y0/2
y1=1-exp(-exp((log(x)-2.5+0.05)/0.2))
y1=y1/2
y2=1-exp(-exp((log(x)-2.5-0.15)/0.2))
y2=y2/2
plot(y0~x,type='l',main='Case 1 (both causes overlap)',lwd=2, xlab=TeX("t"), ylab=TeX("F_k, k=1,2"))
lines(y1~x,lty=2,lwd=2)
lines(y2~x,lty=3,lwd=2)
legend("topleft", legend=c(TeX("Z_1=0, Z_2=0, Z_3=0"),TeX("Z_1=0, Z_2=1, Z_3=0"),TeX("Z_1=0, Z_2=0, Z_3=1"),TeX('F_1 for cause 1'), TeX('F_2 for cause 2')), 
       lty=c(1,2,3,NA,NA),pch=c(NA, NA, NA, 15, 15), col=c('gray0','gray0','gray0','gray0','gray60'),lwd=2)


x=(1:200000)/10000
y0=1-exp(-exp((log(x)-2.5)/0.2))
y0=y0/2
y1=1-exp(-exp((log(x)-2.5+0.05)/0.2))
y1=y1/2
plot(y0~x,type='l',main='Case 4 (both causes overlap)',lwd=2, xlab=TeX("t"), ylab=TeX("F_k, k=1,2"))
lines(y1~x,lty=2,lwd=2)
legend("topleft", legend=c(TeX("Z_1=0, Z_2=0"),TeX("Z_1=0, Z_2=1"),TeX('F_1 for cause 1'), TeX('F_2 for cause 2')), 
       lty=c(1,2,NA,NA),pch=c(NA, NA,  15, 15), col=c('gray0','gray0','gray0','gray60'),lwd=2)



x=(1:300000)/10000
y01=1-(1-0.65*(1-exp(-(x/20)^5)))^exp(0)
y02=1-exp(-exp(0)*(x/20)^5)
y02=y02*(1-0.65)^exp(0)
y11=1-(1-0.65*(1-exp(-(x/20)^5)))^exp(-0.5)
y12=1-exp(-exp(-0.05)*(x/20)^5)
y12=y12*(1-0.65)^exp(-0.5)
y21=1-(1-0.65*(1-exp(-(x/20)^5)))^exp(1)
y22=1-exp(-exp(0.1)*(x/20)^5)
y22=y22*(1-0.65)^exp(1)
plot(y01~x,type='l',main='Case 2',lwd=2, xlab=TeX("t"),ylim=c(0,0.95), ylab=TeX("F_k, k=1,2"))
lines(y02~x,lwd=2,col='gray60')
lines(y11~x,lwd=2,lty=2)
lines(y12~x,lwd=2,lty=2,col='gray60')
lines(y21~x,lwd=2,lty=3)
lines(y22~x,lwd=2,lty=3,col='gray60')
legend("topleft", legend=c(TeX("Z_1=0, Z_2=0, Z_3=0"),TeX("Z_1=0, Z_2=1, Z_3=0"),TeX("Z_1=0, Z_2=0, Z_3=1"),TeX('F_1 for cause 1'), TeX('F_2 for cause 2')), 
       lty=c(1,2,3,NA,NA),pch=c(NA, NA, NA, 15, 15), col=c('gray0','gray0','gray0','gray0','gray60'),lwd=2)

x=(1:300000)/10000
y01=1-(1-0.65*(1-exp(-(x/20)^5)))^exp(0)
y02=1-exp(-exp(0)*(x/20)^5)
y02=y02*(1-0.65)^exp(0)
y11=1-(1-0.65*(1-exp(-(x/20)^5)))^exp(-0.5)
y12=1-exp(-exp(-0.05)*(x/20)^5)
y12=y12*(1-0.65)^exp(-0.5)
plot(y01~x,type='l',main='Case 5',lwd=2, xlab=TeX("t"), ylab=TeX("F_k, k=1,2"))
lines(y02~x,lwd=2,col='gray60')
lines(y11~x,lwd=2,lty=2)
lines(y12~x,lwd=2,lty=2,col='gray60')
legend("topleft", legend=c(TeX("Z_1=0, Z_2=0"),TeX("Z_1=0, Z_2=1"),TeX('F_1 for cause 1'), TeX('F_2 for cause 2')), 
       lty=c(1,2,NA,NA),pch=c(NA, NA,  15, 15), col=c('gray0','gray0','gray0','gray60'),lwd=2)


x=(1:150000)/10000
y0=exp(x-11)/(1+2*exp(x-11))
y1=exp(x-11+0.5)/(1+2*exp(x-11+0.5))
y2=exp(x-11+1)/(1+2*exp(x-11+1))
plot(y0~x,type='l',main='Case 3 (both causes overlap)',lwd=2, xlab=TeX("t"), ylab=TeX("F_k, k=1,2"))
lines(y1~x,lty=2,lwd=2)
lines(y2~x,lty=3,lwd=2)
legend("topleft", legend=c(TeX("Z_1=0, Z_2=0, Z_3=0"),TeX("Z_1=0, Z_2=1, Z_3=0"),TeX("Z_1=0, Z_2=0, Z_3=1"),TeX('F_1 for cause 1'), TeX('F_2 for cause 2')), 
       lty=c(1,2,3,NA,NA),pch=c(NA, NA, NA, 15, 15), col=c('gray0','gray0','gray0','gray0','gray60'),lwd=2)


x=(1:150000)/10000
y0=exp(x-11)/(1+2*exp(x-11))
y1=exp(x-11+1)/(1+2*exp(x-11+1))
plot(y0~x,type='l',main='Case 6 (both causes overlap)',lwd=2, xlab=TeX("t"), ylab=TeX("F_k, k=1,2"))
lines(y1~x,lty=2,lwd=2)
legend("topleft", legend=c(TeX("Z_1=0, Z_2=0"),TeX("Z_1=0, Z_2=1"),TeX('F_1 for cause 1'), TeX('F_2 for cause 2')), 
       lty=c(1,2,NA,NA),pch=c(NA, NA,  15, 15), col=c('gray0','gray0','gray0','gray60'),lwd=2)
