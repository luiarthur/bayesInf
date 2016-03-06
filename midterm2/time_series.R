system("mkdir -p output")
source("../R_Functions/plotPost.R",chdir=T)
source("mcmc/full_like_mcmc_2.R")
source("mcmc/conditional_like_mcmc.R")


y <- read.table("data/dataexam2th_0.txt",header=TRUE)
N <- nrow(y)

pdf("output/rawY.pdf",w=13,h=9)
  plot(ts(y),main="")
dev.off()

pdf("output/lag1.pdf",w=13,h=9)
  par(mfrow=c(3,2))
    sapply(1:ncol(y), function(i) plot(y[-N,i],y[-1,i],xlab="lag 1",ylab=paste0("y",i),pch=20))
  par(mfrow=c(1,1))
dev.off()

tau.sq <- c(.1, 1, 10)

# Conditional Likelihood
set.seed(206)
out_1 <- cond_like_gibbs(y=y,tau2=tau.sq[1],B=2000,burn=1000,printProgress=TRUE)

#plot.posts(out_1$phi_i,names=paste0("y",1:5),cex.a=1.2)
pdf("output/phiv1.pdf",w=13,h=9)
  plot.posts(cbind(out_1$phi,out_1$v),names=c(expression(phi),"v"),cex.a=2,cex.leg=1.5)
dev.off()

pdf("output/phii1.pdf",w=13,h=9)
  par(mfrow=c(2,3))
    sapply(1:ncol(y),function(i) plot.post(out_1$phi_vec[,i],cex.l=1.5,cex.axis=1.5,stay=T,main=paste0(expression(phi),i)))
  par(mfrow=c(1,1))
dev.off()


# Full Likelihood
set.seed(206)
out_2 <- full_like_gibbs(y=y,tau2=tau.sq[1],cand_s=c(.05,.05,.1,.1,.1),B=2000,burn=5000,printProgress=TRUE)

#plot.posts(out_2$phi_vec,names=paste0("y",1:5))
pdf("output/phiv2.pdf",w=13,h=9)
  plot.posts(cbind(out_2$phi,out_2$v),names=c(expression(phi),"v"),cex.a=2,cex.leg=1.5)
dev.off()

pdf("output/phii2.pdf",w=13,h=9)
  par(mfrow=c(2,3))
    sapply(1:ncol(y),function(i) plot.post(out_2$phi_vec[,i],cex.l=1.5,cex.axis=1.5,stay=T,main=paste0(expression(phi),i)))
  par(mfrow=c(1,1))
dev.off()


## Conditional model
ll.cond <- log.like.cond(y,out_1)
plot.post(ll.cond,stay=T)

## Full model
ll.full <- log.like.full(y,out_2)
plot.post(ll.full,stay=T)

## Difference in BIC
plot(density(ll.cond),col="red",xlim=range(c(ll.full,ll.cond)),ylim=c(0,.3))
lines(density(ll.full),col="blue")

pdf("output/bic.pdf",w=13,h=9)
  plot.post(-2*(ll.full-ll.cond),stay=T,main="BIC (Full - Conditional Model)")
dev.off()


# Posterior Predictive: #################################3
#set.seed(206)
#pp.cond <- post.pred.cond(y,out_1$phi,out_1$phi_vec,out_1$v,tau.sq[1])
#
#par(mfrow=c(3,2))
#for (i in 1:5) {
#  bd <- range(apply(pp.cond[[i]],2,function(x) quantile(x,c(.025,.975))))
#  bd <- range(bd,pp.cond[[i]][1,])
#  plot(y[,i],type='l',col="grey",ylab=paste0("y",i),ylim=bd,lwd=1)
#  lines(pp.cond[[i]][1,],col="black",lwd=1)
#  lines(apply(pp.cond[[i]],2,mean),col="blue",lwd=1)
#  lines(apply(pp.cond[[i]],2,function(x)quantile(x,.025)),col="blue",lwd=1)
#  lines(apply(pp.cond[[i]],2,function(x)quantile(x,.975)),col="blue",lwd=1)
#}
#par(mfrow=c(1,1))


set.seed(206)
# tau.sq = 1
out_1.2 <- cond_like_gibbs(y=y,tau2=tau.sq[2],B=2000,burn=1000,printProgress=TRUE)
out_2.2 <- cond_like_gibbs(y=y,tau2=tau.sq[2],B=2000,burn=1000,printProgress=TRUE)
# tau.sq = 10
out_1.3 <- cond_like_gibbs(y=y,tau2=tau.sq[3],B=2000,burn=1000,printProgress=TRUE)
out_2.3 <- cond_like_gibbs(y=y,tau2=tau.sq[3],B=2000,burn=1000,printProgress=TRUE)

pdf("output/priorSensitivity.pdf",w=13,h=9)
par(mfrow=c(2,3))
  plot.post(out_1$phi,main=expression(phi~" for "~tau^2~"=.1, conditional model"),stay=TRUE)
  plot.post(out_1.2$phi,main=expression(phi~" for "~tau^2~"=1, conditional model"),stay=TRUE)
  plot.post(out_1.3$phi,main=expression(phi~" for "~tau^2~"=10 conditional model"),stay=TRUE)
  plot.post(out_2$phi,main=expression(phi~" for "~tau^2~"=.1, full model"),stay=TRUE)
  plot.post(out_2.2$phi,main=expression(phi~" for "~tau^2~"=1, full model"),stay=TRUE)
  plot.post(out_2.3$phi,main=expression(phi~" for "~tau^2~"=10, full model"),stay=TRUE)
par(mfrow=c(1,1))
dev.off()


## Conditional model
ll.cond.2 <- log.like.cond(y,out_1.2)
ll.full.2 <- log.like.full(y,out_2.2)
ll.cond.3 <- log.like.cond(y,out_1.3)
ll.full.3 <- log.like.full(y,out_2.3)

ll.diff.1 <- ll.full-ll.cond
ll.diff.2 <- ll.full.2-ll.cond.2
ll.diff.3 <- ll.full.3-ll.cond.3


# BIC's not different
#plot.post(-2*(ll.diff.1[!is.na(ll.diff.1)]),stay=T,main=expression("BIC (Full - Conditional Model), for "~tau^2~"=.1"))
#plot.post(-2*(ll.diff.2[!is.na(ll.diff.2)]),stay=T,main=expression("BIC (Full - Conditional Model), for "~tau^2~"=1"))
#plot.post(-2*(ll.diff.3[!is.na(ll.diff.3)]),stay=T,main=expression("BIC (Full - Conditional Model), for "~tau^2~"=10"))


