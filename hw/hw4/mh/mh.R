mh_in_gibbs <- function(x, cs_w=2, B=10000, burn=1000) {
  n <- length(x)
  sum_log_x <- sum(log(x))
  sum_x <- sum(x)
  acc_v <- 0

  # Initialize parameters
  w <- 0
  v <- rep(exp(w),B+burn)
  theta <- rep(1,B+burn)

  # log likelihood * prior for w, w = log(v)
  llpw <- function(w,t,a=3,b=1) n*exp(w)*log(t) - n*lgamma(exp(w)) + exp(w)*sum_log_x + a*w - exp(w)*b

  # MCMC Chain
  for (it in 2:(B+burn)) {
    # Update w
    cand <- rnorm(1,w,cs_w)
    lg <- llpw(cand,theta[it-1]) - llpw(w,theta[it-1])
    if (lg > log(runif(1))) {
      w <- cand
      v[it] <- exp(w)
      if (it > burn) acc_v <- acc_v + 1
    } else {
      v[it] <- v[it-1]
    }

    # Update theta
    theta[it] <- rgamma(1,n*v[it]+2,rate=(sum_x+2))
    
    cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\n")

  # Return the chain after burn
  list("v"=tail(v,B),"theta"=tail(theta,B),"acc_v"=acc_v/B)
}
