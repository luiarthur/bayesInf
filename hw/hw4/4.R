#xi | v,t ~ Gamma(v,t), shape and rate
#       v ~ Gamma(3,1)
#       t ~ Gamma(2,2)

source("../../R_Functions/plotPost.R",chdir=T)
source("mh/mh.R")
x <- c( read.table("dat/dat4.dat")$V1 )


#a) 
out <- mh_in_gibbs(x,B=2000,burn=10000,cs_w=.4)
out$acc_v
par(mfrow=c(3,1))
  plot(log(out$v),type='l')
  plot(out$v,type='l')
  plot(out$theta,type='l')
par(mfrow=c(1,1))

# Plot Posterior Predictive with data
plot.posts(cbind(out$v,out$theta),names=c("v","theta"))
plot.post(rgamma(length(out$v),out$v,out$theta),main="Posterior Predictive",stay=T)
lines(density(x),lwd=3,col="grey");

#b)
