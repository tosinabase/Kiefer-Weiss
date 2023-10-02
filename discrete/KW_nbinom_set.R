# set of distribution-specific functions for negative binomial model
#  parametrized by the mean value

source("KW_common.R")

size = 1  # default value of parameter of negative binomial distribution


Hbound <- function(l0, l1, th0, th1, th) {
    # theoretical estimation of the horizon H

  floor(
    (log(l0) * log(th1/th * (th + 1)/(th1 + 1)) + log(l1) * log(th/th0 * (th0 + 1)/(th + 1))) /
    (log((th1 + 1)/(th + 1)) * log(th/th0 * (th0 + 1)/(th + 1))
      - log((th + 1)/(th0 + 1)) * log(th1/th * (th + 1)/(th1 + 1))) / size
  )
}


Lbound <- function(n, l1, th1, th) {
  # a lower bound for  the continuation region at step n

  max(
    ceiling((log(l1) + n * size * log((th + 1)/(th1 + 1))) / (log((th)/(th1) * (th1 + 1)/(th + 1)))),
     0
  )
}


Ubound <- function(n, l0, th0, th) {
  # an upper bound for the continuation region at step n

  floor(
    (log(l0) + n * size * log((th + 1)/(th0 + 1))) / (log((th) / (th0) * (th0 + 1) / (th + 1)))
  )
}


rngen <- function(n, th) {
  rnbinom(n, size, prob = 1 / (1 + th))
}


cdf <- function(x, n, th) {
  pnbinom(x, size * n, prob=1 / (1 + th))
}


pmf <- function(x, n, th) {
  dnbinom(x, size * n, prob=1 / (1 + th))
}

#########################
quant <- function(p, n, th) {
  qnbinom(p, size * n, prob=1 / (1 + th))
}


d <- function(n, s, x) {

  if(n == 1)
    return(1)

  res = choose(size + x - 1, x)
  a = n * size - size
  b = a + s

  for(i in 1:size) {
    res = res * a/b
    a = a + 1
    b = b + 1
  }

  a = s + 1

  if(x > 0) {

    for(i in 1:x){
      res = res * a/b
      a = a + 1
      b = b + 1
    }
  }
  res
}

