cat("Sourcing full_like_mcmc.R\n")


full_like_gibbs <- function(y,tau2,cand_S=diag(7)*.01,B=2000,burn=10000,printProgress=TRUE) {

  logit <- function(x,a,b) log(x-a) - log(b-x)
  logit_inv <- function(x,a,b) (b*exp(x) + a) / (1 + exp(x))
  log_logit_prior <- function(x) x - 2*log(1+exp(x))

  mvrnorm <- function(M,S,n=nrow(S)) M + t(chol(S)) %*% rnorm(n)
  rinvgamma <- function(n,a,b) 1 / rgamma(n,a,rate=b) # a,b are the parameters in wikipedia

  N <- nrow(y)
  I <- ncol(y)

  params <- matrix(0,B+burn,2+I)
  acceptance_rate <- 0

  log_like_plus_log_prior <- function(param_curr) {
    ### Functions for computing log_like_plus_log_prior
    log_like <- function(phi_vec, v) {

      Q_star <- function(i,phi_vec) y[1,i]^2 * (1-phi_vec[i]^2) + sum( (y[-1,i] - phi_vec[i]*y[-N,i])^2 )

      log_like_i <- function(i,phi_vec,v) {
        .5*log(1-phi_vec[i]^2) - N/2*log(v) - Q_star(i,phi_vec) / (2*v)
      }

      sum(  sapply(1:5,function(i) log_like_i(i, phi_vec, v))  )
    }

    #log_prior <- function(log_v,logit_phi,phi_vec) { # normal phi_i
    log_prior <- function(log_v,logit_phi,logit_phi_vec) { # unif phi_i
      log_prior_log_v <- -2 * log_v

      #log_prior_logit_phi <- log_logit_prior(logit_phi) # normal phi_i
      log_prior_logit_phi <- 0 # unif phi_i

      #log_prior_phi_vec <- sum( -(phi_vec-phi)^2 / (2*tau2) ) # normal phi_i
      log_prior_phi_vec <- sum( log_logit_prior(logit_phi_vec) ) # unif phi_i

      log_prior_log_v + log_prior_logit_phi + log_prior_phi_vec
    }

    ### main of log_like_plus_log_prior
    log_v <- param_curr[1]
    v <- exp(log_v)
    logit_phi <- param_curr[2]
    phi <- logit_inv(logit_phi,-1,1)

    #phi_vec <- tail(param_curr,I) # normal phi_i
    logit_phi_vec <- param_curr[-c(1:2)] # unif phi_i
    phi_vec <- logit_inv(logit_phi_vec,-1,1) # unif phi_i

    #log_like(phi_vec,v) + log_prior(log_v, logit_phi, phi_vec) # normal phi_i
    log_like(phi_vec,v) + log_prior(log_v, logit_phi, logit_phi_vec) # unif phi_i
  }
 

  for (it in 2:(B+burn)) {
    params_cand <- mvrnorm(params[it-1,], cand_S)

    accept_ratio <- log_like_plus_log_prior(params_cand) - 
                    log_like_plus_log_prior(params[it-1,])

    if (accept_ratio > log(runif(1))) {
      params[it,] <- params_cand
      if (it > burn) acceptance_rate <- acceptance_rate + 1
    } else {
      params[it,] <- params[it-1,]
    }
    
    if (printProgress) cat("\rProgress: ",it,"/",B+burn)
  }
  
  cat("\nAcceptance: ",acceptance_rate/B,"\n")

  #list("phi_vec"=tail(params[,-c(1:2)],B), # normal phi_i
  list("phi_vec"=logit_inv(  tail(params[,-c(1:2)],B),  -1,1), # unif phi_i
       "phi"=logit_inv(  tail(params[,2],B),  -1,1),
       "v"=exp(  tail(params[,1],B)  ),
       "acceptance_rate"=acceptance_rate)
}
