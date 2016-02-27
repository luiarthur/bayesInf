#xi | v,t ~ Gamma(v,t), shape and rate
#       v ~ Gamma(3,1)
#       t ~ Gamma(2,2)

source("../../R_Functions/plotPost.R",chdir=T)
source("mh/mh.R")
x <- c( read.table("dat/dat4.dat")$V1 )


#a) 
priors <- list("a_v"=3,"b_v"=1,"a_theta"=2,"b_theta"=2)
out <- mh_in_gibbs(x, priors, cs_w=.4, B=2000,burn=10000)
out$acc_v

# Plot Posterior Predictive with data
plot.posts(cbind(out$v,out$theta),names=c("v","theta"))
plot.post(rgamma(length(out$v),out$v,out$theta),main="Posterior Predictive",stay=T)
lines(density(x),lwd=3,col="grey");

#b)
prelim <- mh_multivariate(x, priors, cand_S=.05*diag(2), init=c(1,1), B=2000, burn=10000)
out_b <- mh_multivariate(x, priors, cand_S=cov(prelim$post), init=tail(prelim$post,1), B=2000, burn=1000)
plot.posts(out_b$post,names=c("v","theta"))

#c)
priors <- list("a_v"=3,"b_v"=1,"a_theta"=2,"b_theta"=2)
out_c <- mh_multivariate_laplace_proposal(x, priors, cand_S=cov(prelim$post), init=tail(prelim$post,1), B=2000, burn=50000)
plot.posts(out_c$post)

#d)
priors <- list("a_v"=3,"b_v"=1,"a_theta"=2,"b_theta"=2)
out_d <- mh_multivariate_cauchy_laplace_proposal(x, priors, cand_S=cov(prelim$post), init=tail(prelim$post,1), B=2000, burn=50000)
plot.posts(out_d$post)
