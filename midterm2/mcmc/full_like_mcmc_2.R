cat("Sourcing full_like_mcmc.R\n")


full_like_gibbs <- function(y,tau2,cand_s=rep(.1,ncol(y)),B=2000,burn=10000,printProgress=TRUE) {

  rinvgamma <- function(n,a,b) 1 / rgamma(n,a,rate=b) # a,b are the parameters in wikipedia

  # Precomputes:
  N <- nrow(y)
  I <- ncol(y)
  sqrt_tau2_div_I <- sqrt( tau2 / I )

  # Likelihood and Prior for phi_vec
  Q_star <- function(i,phi_i) y[1,i]^2*(1-phi_i^2) + sum( (y[-1,i] - phi_i*y[-N,i])^2 )
  log_like_plus_log_prior_phi_vec <- function(i,phi_i,v,phi) {
    log_like <- .5*log(1-phi_i^2) - Q_star(i,phi_i)/(2*v)
    log_prior <- -(phi_i-phi)^2/(2*tau2)
    log_like + log_prior
  }

  # Initialize and Preallocate
  phi_vec <- matrix(0,B+burn,I)
  phi <- rep(0,B+burn)
  v <- rep(1,B+burn)
  acc_phi_vec <- rep(0,I)

  for (it in 2:(B+burn)) {
    # Update phi_vec
    for (i in 1:I) {
      curr <- phi_vec[it,i] <- phi_vec[it-1,i]
      cand <- rnorm(1,curr,cand_s[i])
      if (abs(cand) < 1) {
        accept_ratio <- log_like_plus_log_prior_phi_vec(i,cand,v[it-1],phi[it-1]) - 
                        log_like_plus_log_prior_phi_vec(i,curr,v[it-1],phi[it-1])
                        
        if (accept_ratio > log(runif(1))) {
          phi_vec[it,i] <- cand
          if (it > burn) acc_phi_vec[i] <- acc_phi_vec[i] + 1
        }
      }
    }
    # End of Update phi_vec

    # Update phi
    phi[it] <- rnorm(1, mean(phi_vec[it,]), sqrt_tau2_div_I)

    # Update v
    num_v <- sum(  sapply(1:I,function(i) Q_star(i,phi_vec[it,i]))  )
    v[it] <- rinvgamma(1, I*N/2, num_v / 2)
    # End of Update v

    if (printProgress) cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\nAcceptance: ",acc_phi_vec/B,"\n")

  list("phi_vec"=tail(phi_vec,B),"phi"=tail(phi,B),"v"=tail(v,B),"acc_phi_vec"=acc_phi_vec/B)

}
