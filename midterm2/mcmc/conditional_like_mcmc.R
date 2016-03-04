cat("Sourcing conditional_like_mcmc.R\n")

rinvgamma <- function(n,a,b) 1 / rgamma(n,a,rate=b) # a,b are the parameters in wikipedia
logit(x,a,b) <- log(x-a) - log(b-x)
logit_inverse <- function(x,a,b) (b*exp(x) + a) / (1+exp(x))

cond_like_gibbs <- function(y,tau2,B=2000,burn=10000,printProgress=TRUE) {

  N <- nrow(y)
  I <- ncol(y)
  sqrt_tau2_div_I <- sqrt( tau2 / I )

  phi_i <- matrix(0,B+burn,I)
  phi <- rep(0,B+burn)
  v <- rep(1,B+burn)
 
  for (it in 2:(B+burn)) {
    # Update phi_i
    for (i in 1:I) {
      phi_i_denom <- v[it-1] + tau2*sum(y[-N,i]^2)
      phi_i_num_m <- phi[it-1]*v[it-1] + tau2*sum(y[-1,i] * y[-N,i])
      phi_i_num_s <- v[it-1] * tau2
      m <- phi_i_num_m / phi_i_denom
      s <- sqrt( phi_i_num_s / phi_i_denom )
      phi_i[it,i] <- rnorm(1,m,s)
    }

    # Update phi
    phi[it] <- rnorm(1, mean(phi_i[it,]), sqrt_tau2_div_I)

    # update v
    num_v <- 0
    for (i in 1:I) {
      #for (t in 2:N) {
      #  num_v <- num_v + (y[t,i] - phi_i[it,i] * y[t-1,i])^2
      #}
      num_v <- num_v + sum( (y[-1,i] - phi_i[it,i] * y[-N,i])^2 )
    }
    v[it] <- rinvgamma(1, I*N/2, num_v / 2)
    
    if (printProgress) cat("\rProgress: ",it,"/",B+burn)
  }
  cat("\n")

  list("phi_i"=tail(phi_i,B),"phi"=tail(phi,B),"v"=tail(v,B))
}
