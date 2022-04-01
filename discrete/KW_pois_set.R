# set of distribution-specific functions for the Poisson distribution 

source("KW_common.R")


Hbound <- function(l0, l1, th0, th1, th) {
  floor(
    (log(l0) / (log(th) - log(th0)) + log(l1) /  (log(th1) - log(th))) /
    ((th1 - th)/(log(th1) - log(th)) - (th - th0)/(log(th) - log(th0)))
  )
}


Ubound <- function(n, l0, th0, th) {
  # an upper bound for the continuation region at step n
  floor(log(l0) / log(th / th0) + n * (th - th0) / log(th / th0))
}


Lbound <- function(n, l1, th1, th) {
  # a lower bound for  the continuation region at step n
  ceiling(max(0, -log(l1) / log(th1 / th) + n * (th1 - th) / log(th1 / th)))
}


quant <- function(p, n, th) {
  qpois(p, th * n)
}


cdf <- function(s, n, th) {
  ppois(s, th * n)
}


pmf <- function(s, n, th) {
  dpois(s, th * n)
}


rngen <- function(n, th) {
  rpois(n, th)
}


d <- function(n, s, x) {
  if(n == 1)
    return(1)
  choose(x + s, s) * (1 - 1/n)^s/n^x
}
