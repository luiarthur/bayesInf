source("../R_Functions/plotPost.R",chdir=T)
sapply(paste0("mcmc/", system("ls mcmc",intern=TRUE)),source)
y <- read.table("data/dataexam2th_0.txt",header=TRUE)

tau.sq <- c(.1, 1, 10)
prior.1 <- list("phi","v")

out_1 <- cond_like_gibbs(y=y,tau2=tau.sq[3],B=2000,burn=1000,printProgress=TRUE)

plot.post(out_1$phi,main=expression(phi),stay=T)
plot.post(out_1$v,main=expression(nu),stay=T)
plot.posts(out_1$phi_i,names=paste0("y",1:5))

par(mfrow=c(2,3))
  apply(out_1$phi_i,2,function(x) plot.post(x,stay=T))
par(mfrow=c(1,1))

source("mcmc/full_like_mcmc.R")
out_2 <- full_like_gibbs(y=y,tau2=tau.sq[1],cand_S=diag(7)*.005,B=3000,burn=3000,printProgress=TRUE)

plot.post(out_2$phi,main=expression(phi),stay=T)
plot.post(out_2$v,main=expression(nu),stay=T)
plot.posts(out_2$phi_vec,names=paste0("y",1:5))

par(mfrow=c(2,3))
  apply(out_2$phi_vec,2,function(x) plot.post(x,stay=T))
par(mfrow=c(1,1))
