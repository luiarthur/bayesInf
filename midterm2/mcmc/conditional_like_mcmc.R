cat("Sourcing conditional_like_mcmc.R\n")

cond_like_gibbs <- function(y,tau2,B=2000,burn=10000,printProgress=TRUE) {
  rinvgamma <- function(n,a,b) 1 / rgamma(n,a,rate=b) # a,b are the parameters in wikipedia
  logit <- function(x,a,b) log(x-a) - log(b-x)
  logit_inverse <- function(x,a,b) (b*exp(x) + a) / (1+exp(x))

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

  list("phi_vec"=tail(phi_i,B),"phi"=tail(phi,B),"v"=tail(v,B))
}

post.pred.cond <- function(y,phi,phi_vec,v,tau2) {
  B <- length(phi)
  N <- nrow(y)
  I <- ncol(y)

  pred_chain_i <- function(i,b) {
    pred_y <- NULL
    pred_y[1] <- y[1,i]
    for (t in 2:N) {
      pred_y[t] <- rnorm(1, pred_y[t-1] * phi_vec[b,i], sqrt(v[b]))
    }
    pred_y
  }
  
  out <- as.list(1:I)
  for (i in 1:I) {
    out[[i]] <- matrix(0,B,N)
    for (b in 1:B) {
      out[[i]][b,] <- pred_chain_i(i,b)
      cat("\rb: ",b,";  i:",i,"     ")
    }
  }
  cat("\n")

  out
}

log.like.cond <- function(y,out) {
  phi <- out$phi
  v <- out$v
  phi_vec <- out$phi_vec

  B <- length(phi)
  N <- nrow(y)
  I <- ncol(y)

  Q <- function(i,phi_i,b) {
    sum( (y[-1,i] - phi_vec[b,i] * y[-N,i])^2 )
  }

  one_like <- function(b)
    -I*N/2*log(v[b])-sum(sapply(1:I,function(i) Q(i,phi_vec[b,i],b))) /(2*v[b])

  sapply(1:B, one_like)
}
