# Rejection Sampler #############################################################
#
# A function taking the following arguments:
#     log.density.target: a function with one argument x giving the log of the target density at x.
#     log.density.envelope: a function with one argument x giving the log of the envelope density at x.
#     sample.envelope: a function taking no arguments and returning a sample from the envelope distribution.
#     alpha: a numeric such that the density of envelope distribution divided by alpha is greater than
#        the density of the target distribution for all x.
#     sample.size: a numeric given the number of samples to draw.
# The function returns sample.size realizations from the target distribution.

rejection.sampler <- function(log.density.target,log.density.envelope,sample.envelope,alpha,sample.size) {

  # NEW CODE: Aboout 1 second faster for 2000 draws
  n <- sample.size
  x <- matrix(0,n,1)
  if ( is.vector(sample.envelope()) ) {
    x <- matrix(0,n,length(sample.envelope()))
  }

  my.sampler <- function(x){
    out <- sample.envelope()
    while ( log.density.target(out)-(log.density.envelope(out)-log(alpha)) <= log(runif(1)) ){
      out <- sample.envelope()
    }
    out
  }

  t( apply(x,1,my.sampler) )
}
# End of Rejection Sampler ######################################################


# Importance Sampler ############################################################
importance.sampler <- function(p.den,h,q.den,q.sampler,B=1e5,log=F) {
  samps <- q.sampler(B)
  out <- NULL
  if (log) {
    w <- function(x) lp.den(x) - lq.den(x)
    log_sum_w <- log(sum( exp(w(samps)) ))
    out <- exp( sum(log( h(samps) ) + w(samps) - log_sum_w) )
    print("Not Done")
  } else {
    w <- function(x) p.den(x) / q.den(x) 
    sum_w <- sum( w(samps) )
    out <- sum(h(samps) * w(samps)) / sum_w
  }
  out
}
