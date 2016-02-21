# Problems I need to do:
# 1a, 1b, 1c, 2a, 4, 5, 6, 7

# idk: 1b

source("../../myfunc/func.R")

######################################################################

#1) 
params_1 <- matrix(c(10,1,100,10,1000,100),3,2,byrow=T)
colnames(params_1) <- c("n","x")

#1a) Beta(1.5,9.5)
ci_1a <- qbeta(c(.025,.975), params_1[1,2]+.5,params_1[1,1]-params_1[1,2]+.5) # (0.01101167, 0.38131477)
ci_1a

#1c) 
samps_1c <- rbeta(1000000,1.5,9.5)
quantile(samps_1c,c(.025,.975)) # (0.01101609 0.38192279)


onesim <- function(param,B=10^6) {
  ci <- matrix(0,3,2)
  a <- param[2]+.5
  b <- param[1]-param[2]+.5
  ci[1,] <- ci_exact <- qbeta(c(.025,.975),a,b)
  # ci[2,] <- laplace...
  ci[3,] <- ci_mcmc <- quantile( rbeta(B,a,b), c(.025,.975) )
  ci
}

#1d) Need to comment
onesim(params_1[1,])
onesim(params_1[2,])
onesim(params_1[3,])

#2a #######################
dat1 <- c( read.table("dat/dat1.dat")$V1 )
dat2 <- c( read.table("dat/dat2.dat")$V1 )

like <- function(x,m) 
  apply(matrix(1:length(m)),1,function(j) 1/sqrt(2*pi)*exp(-sum((x-m[j])^2)/(2)) )

out1 <- importance.sampler(p.den=function(p) like(dat1,p) * 1,
                          h=function(x) x,
                          q.den=function(x) 1,
                          q.sampler=function(n) rcauchy(n,0,1),B=1000,ret_w=T)
plot(out1$s,out1$w,xlim=c(-2,2),pch=20)
abline(v=out1$mean,lwd=3,col="grey") # -.247
out2 <- importance.sampler(p.den=function(p) like(dat2,p) * 1,
                          h=function(x) x,
                          q.den=function(x) 1,
                          q.sampler=function(n) rcauchy(n,0,1),B=1000,ret_w=T)
plot(out2$s,out2$w,xlim=c(0,20),pch=20)
abline(v=out2$mean,lwd=3,col="grey") #7.96


