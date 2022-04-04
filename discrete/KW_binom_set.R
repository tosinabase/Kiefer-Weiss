# set of distribution-specific functions for the binomial distribution

source("KW_common.R")

size = 1 # parameter of the binomial distribution, the default is size=1 (Bernoulli)


Hbound <- function(l0, l1, th0, th1, th) {
  #for binom distribution

  floor(
    (log(l1)/log(th1/th * (1 - th)/(1 - th1)) - log(l0)/log(th0/th * (1 - th)/(1 - th0))) /
    (log((1 - th0)/(1 - th))/log(th0/th * (1 - th)/(1 - th0)) -
      log((1 - th1)/(1 - th))/log(th1/th * (1 - th)/(1 - th1))
    ) / size
  )
}


Ubound <- function(n, l0, th0, th) {

  floor(
    (-log(l0) - n * size * log((1 - th0)/(1 - th)))/(log(th0/th * (1 - th)/(1 - th0)))
  )
}


Lbound <- function(n, l1, th1, th) {

  max(
    0,
    ceiling(
      (-log(l1) - n * size * log((1 - th1)/(1 - th)))/(log(th1/th * (1 - th)/(1 - th1)))
    )
  )
}


quant <- function(p, n, th) {
  qbinom(p, size * n, th)
}


cdf <- function(s, n, th) {
  pbinom(s, size * n, th)
}


pmf <- function(s, n, th) {
  dbinom(s, size * n, th)
}


rngen <- function(n, th) {
  rbinom(n, size, th)
}


d <- function(n, s, x) {

  if(s > (n - 1) * size)
    return(0)

  if(x > size)
    return(0)

  res = choose(size, x)
  a = n * size - size - s + 1
  b = a + s

  if(x < size)
    for(i in 1:(size - x)) {
      res = res * a / b
      a = a + 1
      b = b + 1
    }

  a = s + 1

  if(x > 0) {
    for(i in 1:x) {
      res = res * a / b
      a = a + 1
      b = b + 1
    }
  }

  res
}

