sapply(paste0("mcmc/", system("ls mcmc",intern=TRUE)),source)
y <- read.table("data/dataexam2th_0.txt",header=TRUE)

tau.sq <- c(.1, 1, 10)
prior.1 <- list("phi","v")

out_1 <- cond_like_gibbs(y=y,tau2=tau.sq[1],B=2000,burn=1000,printProgress=TRUE)

plot(out_1$phi)
plot(out_1$v)
plot(out_1$phi_i[,2])

