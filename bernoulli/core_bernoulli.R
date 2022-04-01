

# number of replicas for Monte Carlo
nMC = 10000


pmf <- function(n, th) {
  # binomial probability mass function

  s = seq(0, n)
  return(dbinom(s, n, th))
}


modified_kw <- function(H, lam0, lam1, th0, th1, th) {
  # Сomputation of the optimal test in the modified Kiefer-Weiss problem.
  # Returns the optimal test (in modified Kiefer-Weiss problem))

  # Parameters:
  # lam0, lam1 Lagrange Multipliers
  # th0, th1 hipothesized values
  # th the point for minimization
  # H horizon

  lagr = list()
  accept = list()
  cont = list()
  z0 = pmf(H, th0)
  z1 = pmf(H, th1)
  lagr[[H]] = pmin(lam0 * z0, lam1 * z1)
  accept[[H]] = lam0 * z0 >= lam1 * z1
  cont[[H]] = rep(FALSE, H + 1)

  if(H > 1)
    for(n in (H - 1):1){
      z0 = pmf(n, th0)
      z1 = pmf(n, th1)
      tmp = pmin(lam0 * z0, lam1 * z1)
      x = seq(0, n + 1)

      h = head(lagr[[n + 1]]/(n + 1) * (n + 1 - x), n + 1)
      t = tail(lagr[[n + 1]]/(n + 1) * (x), n + 1)

      lagr[[n]] = pmin(pmf(n, th) + h + t, tmp)

      accept[[n]] = lam0 * z0 >= lam1 * z1
      cont[[n]] = lagr[[n]] < tmp
    }

  return(rbind(cont, accept))
}


average_sample_number <- function(test, th) {
  #average_sample_number of a test given theta = th

  cont = test[1,]
  H = length(cont)
  asn = list()
  asn[[H]] = rep(0, H + 1)

  if(H == 1)
    return(1)

  for(n in (H - 1):1) {
    x = seq(0, n + 1)


    h = head(asn[[n + 1]]/(n + 1) * (n + 1 - x), n + 1)
    t = tail(asn[[n + 1]]/(n + 1) * (x), n + 1)

    asn[[n]] = ifelse(cont[[n]], (pmf(n, th) + h + t), 0)
  }

  return(1 + head(asn[[1]], 1) + tail(asn[[1]], 1))
}


prob_not_to_stop_before <- function(test, th, k) {
  #probability to stop at  k or thereafter of a test, given theta = th 

  cont = test[1,]
  H = length(cont)

  if(k > H)
    return(0)

  dsn = list()
  dsn[[k]] = rep(0, k + 1)

  if(k > 1){
    for(n in (k - 1):1) {
      x = seq(0, n + 1)
      if(n == k - 1) {

        h = head(dsn[[n + 1]]/(n + 1) * (n + 1 - x), n + 1)
        t = tail(dsn[[n + 1]]/(n + 1) * (x), n + 1)

        dsn[[n]] = ifelse(cont[[n]], (pmf(n, th) + h + t), 0)
      }
      else {

        h = head(dsn[[n + 1]]/(n + 1) * (n + 1 - x), n + 1)
        t = tail(dsn[[n + 1]]/(n + 1) * (x), n + 1)

        dsn[[n]] = ifelse(cont[[n]], (h + t), 0)

      }
    }

    return(head(dsn[[1]], 1) + tail(dsn[[1]], 1))
  }
  else
    return(1)
}


operating_characteristic <- function(test, th) {
  # operating characteristic of a test given theta = th

  cont = test[1,]
  accept = test[2,]

  H = length(cont)
  a = list()
  a[[H]] = ifelse(accept[[H]], pmf(H, th),0)

  if(H > 1)
    for(n in (H - 1):1){
      x = seq(0, n + 1)

      h = head(a[[n + 1]]/(n + 1) * (n + 1 - x), n + 1)
      t = tail(a[[n + 1]]/(n + 1) * (x), n + 1)

      a[[n]] = ifelse(cont[[n]], h + t, ifelse(accept[[n]], pmf(n, th), 0))
    }

  return(head(a[[1]], 1) + tail(a[[1]], 1))
}


monte_carlo_simulation <- function(test, hyp) {
  # optimal test simulation
  # hyp = true success probability
  # returns rate of accepting
  # the average run length  and its standard deviation

  cont = test[1,]
  accept = test[2,]

  H = length(cont)
  ss = 0
  totaccepted = 0
  totn = 0

  for(i in 1 : nMC) {
    s = 0 # accumulated sum for test run
    n = 0 # number of steps in current run
    accepted = 0 # number of accepting in the run (0 or 1)

    for (stage in 1:H){
      # run starts

      # generate
      summand = rbinom(1, 1, hyp) # 1 bernoulli
      s = s + summand # accumulated sum
      n = n + 1   # one step more

      if(stage == H) {
        # the last stage of the run
        if(accept[[stage]][s + 1])
          # if s is below of the decision constant then accept
          accepted = accepted + 1
      }
      else {
        if(cont[[stage]][s + 1] == FALSE){
          # stop by the optimal stopping rule
          if (accept[[stage]][s + 1])
            accepted = accepted + 1

          break # accepted or rejected, stop the run; no more stages
        }

      }

    }
    totaccepted = totaccepted + accepted
    totn = totn + n
    ss = ss + n^2

  }

  OC = (totaccepted)/as.double(nMC)
  ASN = totn/as.double(nMC)
  stdev = sqrt((ss - totn^2 / as.double(nMC)) / as.double(nMC - 1))

  return(c(OC, ASN, stdev))
}


original_kw <- function(H = 1, option, lam0, lam1, th0, th1, tol=.Machine$double.eps^0.25) {
  # minimizing  maximum average_sample_number cycle
  # H only used when option != 1

  delta <- function(x){

    if(option == 1)
      H=Hbound(lam0,lam1,th0,th1,x)

    test = modified_kw(H, lam0, lam1, th0, th1, x)
    ASNold = average_sample_number(test, x)
    a = optimize(average_sample_number, c(th0, th1), maximum=TRUE, test=test, tol=tol)
    return(a$objective - ASNold)
  }

  res = optimize(delta, c(th0, th1), maximum=FALSE, tol=tol)
  if(option == 1)
    H = Hbound(lam0, lam1, th0, th1, res$minimum)

  return(c(res$minimum, res$objective))
}


max_number <- function(test) {
  H = length(test[1,])
  n = 1
  repeat {
    if (all(!test[1,][[n]]))
      break
    if(n == H)
      break

    n = n+1
  }
  return(n)
}


Hbound <- function(lam0, lam1, th0, th1, th) {
  # Upper bound for the maximum sample number
  # for the optimal test in the modified Kiefer-Weiss problem, see Hawix and Schmitz

  y = c(log(th / th0 / (1 - th) * (1 - th0)), log(th / th1 / (1 - th) * (1 - th1)))
  z = c(log((1 - th) / (1 - th0)), log((1 - th) / (1 - th1)))
  a = solve(rbind(y, z))
  b = a%*%(c(0, 1))
  c = sum(ifelse(dbinom(c(0, 1), 1, th0) < dbinom(c(0, 1), 1, th1), dbinom(c(0, 1), 1, th0), 0))
  d = sum(ifelse(dbinom(c(0, 1), 1, th0) >= dbinom(c(0, 1), 1, th1), dbinom(c(0, 1), 1, th1), 0))

  return(ceiling(b[1] * log(lam0) + b[2] * log(lam1) - (b[1] + b[2]) * (-log(1 - c - d))))
}


sample_number_quantile <- function(test, th, q) {
  # q-quantile of the sample number at point th

  H = length(test[1,])

  k = H
  repeat{
    if(prob_not_to_stop_before(test, th, k) > 1 - q)
      break
    k = k - 1
  }

  return(k)
}
