#  KW_common.R contains common functions for any Kiefer-Weiss problem
#  model-specific functions are placed in a separate module, like KW_pois_set.R, KW_binom_set.R, etc.


back_step_int <- function(stepdata, n, s, l0, l1, th0, th1) {

  incorp <- function(x) {
    if(!stepdata$laststep){

      if(x >= stepdata$from & x-stepdata$from + 1 <= stepdata$length)
        return(stepdata$val[x - stepdata$from + 1])
      else
        return(min(l0 * pmf(x, n, th0), l1 * pmf(x, n, th1)))

    }
    else
      return(min(l0 * pmf(x, n, th0), l1 * pmf(x, n, th1)))
  }

  sum = 0
  k = 0
  repeat{
    sumold = sum
    sum = sum + incorp(s + k) * d(n, s, k)
    if(sum == sumold)
      break
    k = k + 1
  }
  sum
}


back_step_int_oc <- function(stepdata, s, th) {
  # backward induction step summation
  if(!stepdata$laststep) {

    sum = cdf(min(stepdata$from - 1 - s,stepdata$acceptAt - s), 1, th)
    for(k in seq(stepdata$from - s, length.out=stepdata$length))
      if(k >= 0)
        sum = sum + stepdata$val[k + s - stepdata$from + 1] * pmf(k, 1, th)

  }
  else
    sum = cdf(stepdata$acceptAt - s, 1, th)

  sum
}


back_step_int_asn <- function(stepdata, s, th) {
  # backward induction step summation

  if(!stepdata$laststep) {
    k = 0
    sum = 0
    repeat{
      if(k + s >= stepdata$from & k + s <= stepdata$from + stepdata$length - 1) {
        sum = sum + stepdata$val[k + s - stepdata$from + 1] * pmf(k, 1, th)
      }

      if(k + s >= stepdata$from + stepdata$length - 1)
        break

      k = k + 1
    }

  }
  else
    sum = 0
  sum
}


step_effect <- function(stepdata, n, s, l0, l1, th0, th1, th) {
  #will continue if negative
  back_step_int(stepdata, n + 1, s, l0, l1, th0, th1) +
    pmf(s, n, th) - min(l0 * pmf(s, n, th0), l1 * pmf(s, n, th1))
}


operating_characteristic <- function(test, th) {

  H = length(test)
  stepdata = test[[H]]
  H = H - 1

  repeat{

    for(i in seq(test[[H]]$from,length.out=test[[H]]$length)){
      test[[H]]$val[i - test[[H]]$from + 1] = back_step_int_oc(stepdata, i, th)
    }
    stepdata = test[[H]]
    if(H == 1)
      break
    H = H - 1
  }
  back_step_int_oc(stepdata, 0, th)
}


average_sample_number <- function(test, th) {
  H = length(test)
  stepdata = test[[H]]

  H=H-1
  repeat{
    for(i in seq(test[[H]]$from,length.out=test[[H]]$length)) {
      test[[H]]$val[i - test[[H]]$from + 1]= back_step_int_asn(stepdata, i, th) + 1
    }
    stepdata = test[[H]]
    if(H == 1)
      break
    H = H - 1
  }
  back_step_int_asn(stepdata, 0, th) + 1
}


Lagr <- function(test, l0, l1, th0, th1, th) {

  H = length(test)
  stepdata = test[[H]]
  H = H - 1

  repeat{
    for(i in seq(test[[H]]$from, length.out=test[[H]]$length)) {
     test[[H]]$val[i - test[[H]]$from + 1] = back_step_int(stepdata, H + 1, i, l0, l1, th0, th1) + pmf(i, H, th)
    }
    stepdata=test[[H]]

    if(H == 1)
      break
    H = H - 1

  }
  back_step_int(stepdata, 1, 0, l0, l1, th0, th1) + 1
}


monte_carlo_simulation <- function(test, hyp, nMC=10000) {
  # optimal test simulation
  # hyp = true parameter value
  # nMC number of Monte Carlo replications 
  # returns relative frequency of acceptations
  # the average run length  and its standard deviation

  H = length(test)
  const = test[[H]]$acceptAt
  H = H - 1
  ss = 0
  totaccepted = 0
  totn = 0
  for(i in 1:nMC) {
    s = 0  # accumulated sum for test run
    n = 0  # number of steps in current run
    accepted = 0  #number of accepting in the run (0 or 1)
    for (stage in 1:H){
      # run starts
      # generate
      summand = sum(rngen(1, hyp))
      s = s + summand  # accumulated sum
      n = n + 1  # one step more
      if(stage == H){
        # the last stage of the run
        if(s <= const)
          accepted = accepted + 1  # if s is below of the decision constant then accept
      }
      else {
        if(s < test[[stage]]$from| s > test[[stage]]$from + test[[stage]]$length - 1){
           # stop by the optimal stopping rule
          if (s <= test[[stage]]$acceptAt)
            accepted = accepted + 1  # if s is in the acceptation
          # region then accept
          break #  accepted or rejected, stop the run; no more stages
        }
      }

    }
    totaccepted=totaccepted+accepted
    totn=totn+n
    ss=ss+n^2
  }
  OC = (totaccepted)/as.double(nMC)
  ASN = totn/as.double(nMC)
  stdev = sqrt((ss - totn^2/as.double(nMC))/as.double(nMC - 1))
  return(c(OC, ASN, stdev))
}


modified_kw <- function(l0, l1, th0, th1, th, H=0) {


  #  returns optimal test in modified Kiefer Weiss problem
  #  truncation level given by H if H is not zero
  fill_in <- function(stepdata, a, b, H) {
    new = list(H=H, from=a, length=b - a + 1, val=array(dim=b - a + 1), laststep=FALSE)
    for(i in seq(new$from, length.out=new$length)){
      new$val[i - new$from + 1] = min(
        back_step_int(stepdata, H + 1, i, l0, l1, th0, th1) + pmf(i, H, th),
        l0 * pmf(i, H, th0),
        l1 * pmf(i, H, th1)
      )
    }
    new$acceptAt = Ubound(H, l0/l1, th0, th1)
    return(new)
  }

  if(th < th1 & th > th0){
    if(H==0)
      h = Hbound(l0, l1, th0, th1, th)
    else
      h = H
    while(Lbound(h, l1, th1, th) > Ubound(h, l0, th0, th))
      h = h - 1

    test = list()
    stepdata = list(H=h + 1, acceptAt=Ubound(h + 1, l0/l1, th0, th1), laststep=TRUE)
    test[[h + 1]] = stepdata
    repeat{
      nocont = TRUE
      for(i in seq(Lbound(h, l1, th1, th), Ubound(h, l0, th0, th))) {
        if(step_effect(stepdata, h, i, l0, l1, th0, th1, th) < 0) {
          nocont = FALSE;
          a = i;
          break
        }
      }
      if(nocont){
        h = h - 1
        if(h == 0)
          stop("stop after one observation")
        rm(test);
        test = list()
        stepdata = list(H=h + 1, acceptAt=Ubound(h + 1, l0/l1, th0, th1), laststep=TRUE)
        test[[h + 1]] = stepdata
        next
      }

      s = Ubound(h, l0, th0, th)
      repeat{
        if(step_effect(stepdata, h, s, l0, l1, th0, th1, th) < 0)
          break
        s = s - 1
      }
      b = s
      stepdata = fill_in(stepdata, a, b, h)
      test[[h]] = stepdata
      if(h == 1){
        break
      }
      h = h - 1
    }
  }

  if(th >= th1){
    if(H <= 0)
      stop("need to specify horizon")

    test = list()
    h = H
    stepdata = list(H=h + 1,acceptAt=Ubound(h + 1, l0/l1, th0, th1), laststep=TRUE)
    test[[h + 1]] = stepdata
    repeat{

      s = Ubound(h, l0, th0, th)
      if(th > th1) s = min(s, Lbound(h, l1, th1, th))
      nocont = TRUE

      repeat{
        if(step_effect(stepdata, h, s, l0, l1, th0, th1, th) < 0){
          nocont = FALSE
          b = s
          break
        }
        if(s == 0)
          break
        s = s - 1
      }
      if(nocont){
        h = h - 1;
        if(h == 0)
          stop("stop after one observation")
        rm(test);
        test = list()
        stepdata = list(H=h + 1,acceptAt=Ubound(h + 1, l0/l1, th0, th1),laststep=TRUE)
        test[[h + 1]] = stepdata
        next
      }

      if(th>th1) {
        s = 0
        repeat{
           if(step_effect(stepdata, h, s, l0, l1, th0, th1, th) < 0) {
             a = s;
             break
           }
        s = s + 1
        }
      }
      else {
        s = b - 1
        repeat {
          if(step_effect(stepdata, h, s, l0, l1, th0, th1, th) >= 0) {
            a = s + 1;
            break
          }
          if(s == 0) {
            a = 0;
            break
          }
          s = s - 1
        }
      }
      stepdata = fill_in(stepdata, a, b, h)
      test[[h]] = stepdata
      if(h == 1)
        break
      h = h - 1
    }
  }

  if(th <= th0) {
    if(H <= 0)
      stop("need to specify horizon")
    test = list()
    h = H
    stepdata = list(H=h + 1,acceptAt=Ubound(h + 1, l0/l1, th0, th1), laststep=TRUE)
    test[[h + 1]] = stepdata

    repeat{
      s = Lbound(h, l1, th1, th)
      if(th < th0)
        s = max(s, Ubound(h, l0, th0, th))
      nocont = TRUE
      repeat {
        if(step_effect(stepdata, h, s, l0, l1, th0, th1, th) < 0) {
          nocont = FALSE
          a = s
          break
        }
        s = s + 1
      }

      if(nocont){
        rm(test);
        test=list();
        h = h - 1;
        stepdata = list(H=h + 1, acceptAt=Ubound(h + 1, l0/l1, th0, th1), laststep=TRUE)
        test[[h + 1]] = stepdata
        if(h == 0)
          stop("stop after one observation")
        next
      }
      s = a + 1
      repeat{
        if(step_effect(stepdata,h,s,l0,l1,th0,th1,th)>=0){b=s-1;break}
        s = s + 1
      }

      stepdata = fill_in(stepdata, a, b, h)
      test[[h]] = stepdata
      if(h == 1)
        break
      h = h - 1
    }
  }

  test
}


prob_to_stop_after <- function(test, th, k) {
  #probability to stop after step  k for a test, given theta = th

  if(k >= length(test))
    return(0)
  if(k <= 0)
    return(1)
  test[[k]]$val = rep(1, test[[k]]$length)
  stepdata = test[[k]]
  if(k > 1) {
    n = k-1
    repeat {
      for(i in seq(test[[n]]$from, length.out=test[[n]]$length)) {
        test[[n]]$val[i - test[[n]]$from + 1] = back_step_int_asn(stepdata, i, th)
      }
      stepdata = test[[n]]
      if(n == 1)
        break
      n = n - 1
    }
  }
  return(back_step_int_asn(stepdata, 0, th))
}


sample_number_quantile <- function(test, th, q) {
  # q-quantile of the sample number given theta=th
  H = length(test)
  n1 = 1
  n2 = H + 1
  repeat{
    n = floor((n1 + n2)/2)
    if(prob_to_stop_after(test, th, n - 1) > 1 - q)
      n1 = n
    else
      n2 = n

    if(n2-n1<=1)
      break
  }
  return(n1)
}


original_kw <- function(lam0, lam1, th0, th1, H=0, tol=.Machine$double.eps^0.25) {
  # minimizing  maximum average_sample_number cycle

  option = ifelse(H == 0, 1, 2)

  delta <- function(x) {

    if(option == 1)
      H = Hbound(lam0, lam1, th0, th1, x)

    test = modified_kw(lam0, lam1, th0, th1, x, H=H)
    ASNold = average_sample_number(test, x)
    a = optimize(average_sample_number, c(th0, th1), maximum=TRUE, test=test, tol=tol)

    return(a$objective - ASNold)
  }

  res = optimize(delta, c(th0, th1), maximum=FALSE, tol=tol)

  return(c(res$minimum, res$objective))
}


tsprt <- function(H, A, B, th0, th1) {

  test = list()
  test[[H + 1]] = list(H=H + 1, acceptAt=Ubound(H + 1, A * B, th0, th1), laststep=TRUE)
  for(i in seq(1, H)){
    from = Lbound(i, 1/A, th1, th0)
    to = Ubound(i, B, th0, th1)
    test[[i]] = list(H=i, from=from, length=to - from + 1, var=array(dim=to - from + 1),
                     laststep=FALSE, acceptAt=Ubound(i, A * B, th0, th1))
  }
  test
}


sprt <- function(A, B, th0, th1) {

  H = 16 * ceiling(np((1 - A)/(B - A),A * (B - 1)/(B - A), th0, th1))
  ASN0 = 0
  ASN1 = 0

  repeat{
    ASN0old = ASN0
    ASN1old = ASN1
    test = tsprt(H, A, B, th0, th1)
    ASN0 = average_sample_number(test, th0)
    ASN1 = average_sample_number(test, th1)
    if(ASN0 == ASN0old & ASN1 == ASN1old)
      break
    H = H * 2
  }
  test
}


np_test <- function(alpha, n, th0, th1) {
  # type II error probability of the Neyman-Pearson test
  #based on n observations
  #returns vector of critical constant, radndomization constant and the type II error probability

  c = quant(1 - alpha, n, th0)
  gam = (1 - alpha - cdf(c - 1, n, th0))/pmf(c, n, th0)
  return(c(c, gam, cdf(c - 1, n, th1) + gam * pmf(c, n, th1)))
}


np <- function(alpha, beta, th0, th1) {
  # minimum fixed sample size test
  n = 1
  repeat{
    n = n * 2
    if(np_test(alpha, n, th0, th1)[3] < beta)
      break
  }
  n1 = 1
  n2 = n
  repeat{
    n = floor((n1 + n2)/2)
    if(np_test(alpha, n, th0, th1)[3] < beta)
      n2 = n
    else
      n1 = n

    if(n2 - n1 <= 1)
      break
  }
  ant = np_test(alpha, n1, th0, th1)[3]
  cur = np_test(alpha, n2, th0, th1)[3]
  return(n2 - 1 + (beta - ant)/(cur - ant))
}

