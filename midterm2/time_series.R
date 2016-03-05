system("mkdir -p output")
source("../R_Functions/plotPost.R",chdir=T)
source("mcmc/full_like_mcmc_2.R")
source("mcmc/conditional_like_mcmc.R")


y <- read.table("data/dataexam2th_0.txt",header=TRUE)

tau.sq <- c(.1, 1, 10)

# Conditional Likelihood
set.seed(206)
out_1 <- cond_like_gibbs(y=y,tau2=tau.sq[3],B=2000,burn=1000,printProgress=TRUE)

#plot.posts(out_1$phi_i,names=paste0("y",1:5),cex.a=1.2)
pdf("output/phiv1.pdf",w=13,h=9)
  plot.posts(cbind(out_1$phi,out_1$v),names=c(expression(phi),"v"),cex.a=2,cex.leg=1.5)
dev.off()

pdf("output/phii1.pdf",w=13,h=9)
  par(mfrow=c(2,3))
    sapply(1:ncol(y),function(i) plot.post(out_1$phi_i[,i],cex.l=1.5,cex.axis=1.5,stay=T,main=paste0(expression(phi),i)))
  par(mfrow=c(1,1))
dev.off()


# Full Likelihood
set.seed(206)
out_2 <- full_like_gibbs(y=y,tau2=tau.sq[1],cand_s=c(.05,.05,.1,.1,.1),B=2000,burn=5000,printProgress=TRUE)

#plot.posts(out_2$phi_vec,names=paste0("y",1:5))
pdf("output/phiv2.pdf",w=13,h=9)
  plot.posts(cbind(out_2$phi,out_1$v),names=c(expression(phi),"v"),cex.a=2,cex.leg=1.5)
dev.off()

pdf("output/phii2.pdf",w=13,h=9)
  par(mfrow=c(2,3))
    sapply(1:ncol(y),function(i) plot.post(out_2$phi_vec[,i],cex.l=1.5,cex.axis=1.5,stay=T,main=paste0(expression(phi),i)))
  par(mfrow=c(1,1))
dev.off()

