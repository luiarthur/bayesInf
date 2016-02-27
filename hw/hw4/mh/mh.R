mvrnorm <- function(M,S,n=nrow(S)) M + t(chol(S)) %*% rnorm(n)

# http://mathoverflow.net/questions/66307/generating-samples-of-a-multivariate-cauchy-distribution
mvrcauchy <- function(m,S) {
  z <- rchisq(1,1)
  mvrnorm(m,S/z)
}

# a) Regular Metropolis for v
mh_in_gibbs <- function(x, priors, cs_w=2, B=10000, burn=1000) {
  n <- length(x)
  sum_log_x <- sum(log(x))
  sum_x <- sum(x)
  acc_v <- 0

  # Initialize parameters
  w <- 0
  v <- rep(exp(w),B+burn)
  theta <- rep(1,B+burn)

  # log likelihood * prior for w, w = log(v)
  llpw <- function(w, t, a=priors$a_v, b=priors$b_v) 
    n*exp(w)*log(t) - n*lgamma(exp(w)) + exp(w)*sum_log_x + a*w - exp(w)*b

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
    theta[it] <- rgamma(1,n*v[it]+priors$a_theta,rate=(sum_x+priors$b_theta))
    
    cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\n")

  # Return the chain after burn
  list("v"=tail(v,B),"theta"=tail(theta,B),"acc_v"=acc_v/B)
}

# b) Multivariate proposal for log(v), log(theta) --- easy
mh_multivariate <- function(x, priors, cand_S, init, B=10000, burn=1000) {

  # Precomputes
  n <- length(x)
  sum_log_x <- sum(log(x))
  sum_x <- sum(x)
  a_v <- priors$a_v
  b_v <- priors$b_v
  a_theta <- priors$a_theta
  b_theta <- priors$b_theta
  p <- ncol(cand_S)

  # Proposal Sampler
  propose <- function(current_param) mvrnorm(current_param,cand_S)

  # log likelihood + log prior
  log_lik_plus_prior <- function(param) {
    log_v <- param[1] 
    log_theta <- param[2] 
    v <- exp(log_v)
    theta <- exp(log_theta)

    log_dprior <- function(w,a,b) a*w - b*exp(w)
    
    log_like <- n*v*log(theta) - n*lgamma(v) + v*sum_log_x - theta*sum_x
    log_prior <- log_dprior(log_v,a_v,b_v) + log_dprior(log_theta,a_theta,b_theta)

    log_like + log_prior
  }

  # Initialize
  posterior_draws <- matrix(0,B+burn,p)
  colnames(posterior_draws) <- c("v","theta")
  posterior_draws[1,] <- init
  acceptance_rate <- 0

  for (it in 2:(B+burn)) {
    cand <- propose(posterior_draws[it-1,])
    acceptance_ratio <- log_lik_plus_prior(cand) - 
                        log_lik_plus_prior(posterior_draws[it-1,])
    if (acceptance_ratio > log( runif(1) )) {
      posterior_draws[it,] <- cand
      if (it > burn) acceptance_rate <- acceptance_rate + 1
    } else {
      posterior_draws[it,] <- posterior_draws[it-1,]
    }
    cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\nAcceptance Ratio: ", acceptance_rate/B,"\n")
  

  settings <- list("priors"=priors,"cand_S"=cand_S,"init"=init)
  list("posterior"=exp(tail(posterior_draws,B)), "accept"=acceptance_rate/B, "settings"=settings)
}


# c) Multivariate proposal for log(v), log(theta) with Laplace Approx
mh_multivariate_laplace_proposal <- function(x, priors, cand_S, init, B=10000, burn=1000) {

  # Precomputes
  n <- length(x)
  sum_log_x <- sum(log(x))
  sum_x <- sum(x)
  a_v <- priors$a_v
  b_v <- priors$b_v
  a_theta <- priors$a_theta
  b_theta <- priors$b_theta
  p <- ncol(cand_S)

  # log likelihood + log prior
  log_lik_plus_prior <- function(param) {
    log_v <- param[1] 
    log_theta <- param[2] 
    v <- exp(log_v)
    theta <- exp(log_theta)

    log_dprior <- function(w,a,b) a*w - b*exp(w)
    
    log_like <- n*v*log(theta) - n*lgamma(v) + v*sum_log_x - theta*sum_x
    log_prior <- log_dprior(log_v,a_v,b_v) + log_dprior(log_theta,a_theta,b_theta)

    log_like + log_prior
  }

  # Proposal Sampler
  phi <- optim(c(0,0),log_lik_plus_prior,control=list("fnscale"=-1))$par #maximize
  H <- optimHess(c(0,0),log_lik_plus_prior)
  I <- solve(-H)
  propose <- function() mvrnorm(phi,I)

  # Initialize
  posterior_draws <- matrix(0,B+burn,p)
  colnames(posterior_draws) <- c("v","theta")
  posterior_draws[1,] <- init
  acceptance_rate <- 0

  for (it in 2:(B+burn)) {
    cand <- propose()
    acceptance_ratio <- log_lik_plus_prior(cand) - 
                        log_lik_plus_prior(posterior_draws[it-1,])
    if (acceptance_ratio > log( runif(1) )) {
      posterior_draws[it,] <- cand
      if (it > burn) acceptance_rate <- acceptance_rate + 1
    } else {
      posterior_draws[it,] <- posterior_draws[it-1,]
    }
    cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\nAcceptance Ratio: ", acceptance_rate/B,"\n")
  

  settings <- list("priors"=priors,"cand_S"=cand_S,"init"=init)
  list("posterior"=exp(tail(posterior_draws,B)), "accept"=acceptance_rate/B, "settings"=settings)
}

# d) Multivariate proposal for log(v), log(theta) with Cauchy Laplace Approx
mh_multivariate_cauchy_laplace_proposal <- function(x, priors, cand_S, init, B=10000, burn=1000,printProgress=T) {

  # Precomputes
  n <- length(x)
  sum_log_x <- sum(log(x))
  sum_x <- sum(x)
  a_v <- priors$a_v
  b_v <- priors$b_v
  a_theta <- priors$a_theta
  b_theta <- priors$b_theta
  p <- ncol(cand_S)

  # log likelihood + log prior
  log_lik_plus_prior <- function(param) {
    log_v <- param[1] 
    log_theta <- param[2] 
    v <- exp(log_v)
    theta <- exp(log_theta)

    log_dprior <- function(w,a,b) a*w - b*exp(w)
    
    log_like <- n*v*log(theta) - n*lgamma(v) + v*sum_log_x - theta*sum_x
    log_prior <- log_dprior(log_v,a_v,b_v) + log_dprior(log_theta,a_theta,b_theta)

    log_like + log_prior
  }

  # Proposal Sampler

  phi <- optim(c(0,0),log_lik_plus_prior,control=list("fnscale"=-1))$par #maximize
  H <- optimHess(c(0,0),log_lik_plus_prior)
  I <- solve(-H)
  propose <- function() { 
    out <- c( mvrcauchy(phi,I) )
    ifelse(abs(out) > 100, 100, out)
  }

  # Initialize
  posterior_draws <- matrix(0,B+burn,p)
  colnames(posterior_draws) <- c("v","theta")
  posterior_draws[1,] <- init
  acceptance_rate <- 0

  for (it in 2:(B+burn)) {
    cand <- propose()
    acceptance_ratio <- log_lik_plus_prior(cand) - 
                        log_lik_plus_prior(posterior_draws[it-1,])
    if (acceptance_ratio > log( runif(1) )) {
      posterior_draws[it,] <- cand
      if (it > burn) acceptance_rate <- acceptance_rate + 1
    } else {
      posterior_draws[it,] <- posterior_draws[it-1,]
    }
    if (printProgress) cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\nAcceptance Ratio: ", acceptance_rate/B,"\n")
  

  settings <- list("priors"=priors,"cand_S"=cand_S,"init"=init)
  list("posterior"=exp(tail(posterior_draws,B)), "accept"=acceptance_rate/B, "settings"=settings)
}
