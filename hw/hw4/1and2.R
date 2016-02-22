# Problems I need to do:
# 1a, 1b, 1c, 2a, 4, 5, 6, 7

# idk: 1b

source("../../myfunc/func.R")

######################################################################

#1a) Beta(1.5,9.5)
x <- 1; n <- 10
ci_1a <- qbeta(c(.025,.975), x+.5, n-x+.5) # (0.01101167, 0.38131477)
ci_1a

#1b) Robert 6.2.9
a <- x+.5
b <- n-x+.5
post.mean <- a / (a+b)
post.var <- a*b / ( (a+b)^2 * (a+b+1) )
ci_1b <- qnorm(c(.025,.975), post.mean, sqrt(post.var)) # (-0.05780193, 0.33052920)
ci_1b

rm(list=c("a","b")) # removes the variables a, b

#http://www.sumsar.net/blog/2013/11/easy-laplace-approximation/
#(1.5*.5*.1^-.5  + 8.5*7.5*.9^6.5 - 2*1.5*8.5*.1^.5*.9^7.5) * .1^2 /
#(.5*-.5*.1^-1.5 + 8.5*7.5*.9^6.5 - 2* .5*8.5*.1^-.5*.9^7.5) # .1961632
#
#pp <- seq(0,1,len=100)
#curve(pp * function(p) dbinom(1,10,p)
#
#laplace.approx <- function(bn,bd,n,tn,td,s2_n,s2_d,hn,hd) {
#  bn(tn) / bd(td) * s2_n(tn,n) / s2_d(td,n) * exp(-(hn(tn,n)-hd(td,n)))
#}
#
#fd <- function(p) dbinom(1,10,p) * dbeta(p,.5,.5)
#fn <- function(p) p * fd(p)
#curve(fn,from=0,to=1,col='blue',ylim=c(0,.5))
#curve(fd,from=0,to=1,add=T,col='red')
#pp <- seq(0,1,len=500)
#td <- pp[which.max(fd(pp))]
#tn <- pp[which.max(fn(pp))]
#
#laplace.approx(bn=function(p) 1,bd=function(p) 1,n=10,td=td,tn=tn,
#               s2_n=function(p,n) 1/ ((1.5/p^2 + 8.5/(1-p)^2)),
#               s2_d=function(p,n) 1/ ((0.5/p^2 + 8.5/(1-p)^2)),
#               hn=function(p,n) -( dbinom(1,n,p,log=T) + dbeta(p,1+.5+1,n-1+.5,log=T) ),
#               hd=function(p,n) -( dbinom(1,n,p,log=T) + dbeta(p,1+.5+0,n-1+.5,log=T) ))

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


