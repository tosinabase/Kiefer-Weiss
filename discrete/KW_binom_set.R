# set of distribution-specific functions for the binomial distribution

source("KW_common.R")

size = 1 # parameter of the binomial distribution, the default is size=1 (Bernoulli)


Hbound <- function(l0, l1, th0, th1, th) {
    # theoretical estimation of the horizon H

  floor(
    (log(l1)/log(th1/th * (1 - th)/(1 - th1)) - log(l0)/log(th0/th * (1 - th)/(1 - th0))) /
    (log((1 - th0)/(1 - th))/log(th0/th * (1 - th)/(1 - th0)) -
      log((1 - th1)/(1 - th))/log(th1/th * (1 - th)/(1 - th1))
    ) / size
  )
}


Ubound <- function(n, l0, th0, th) {
    # an upper bound for the continuation region at step n

  floor(
    (-log(l0) - n * size * log((1 - th0)/(1 - th)))/(log(th0/th * (1 - th)/(1 - th0)))
  )
}


Lbound <- function(n, l1, th1, th) {
  # a lower bound for  the continuation region at step n

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

# input parameters can be set here
l0 = 449.995105192065
l1 = 489.748951643705

th0 = 0.05
th1 = 0.08
th = 0.6

H = 30
test = modified_kw(l0, l1, th0, th1, th, H=H)
print("asd")

for (i in 1:H){
  test[[i]]
}