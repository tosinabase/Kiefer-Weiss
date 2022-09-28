# install.packages("pracma")
library(pracma)


intgr <- function(spts, vls, t){
  # Backward induction integration
  # numerical integral of \int u(x)*dnorm(x-t) dx
  # spts is grid of x values, vls corresponding values of u
  delta = spts[2] - spts[1]
  PhiVec = pnorm(spts - t[1])
  phiVec = dnorm(spts - t[1])
  diffPhi = diff(PhiVec)
  diffphi = diff(phiVec)
  tmp = sum(head(vls, -1) * diffPhi)
  vec = diffPhi * head(spts - t[1], -1) + diffphi
  tmp = tmp - sum((diff(vls)/delta) * vec)
  return(tmp)
}

intgr_trapezoidal <- function (spts, vls){
  # numerical integral of \int f(x) dx using trapezoidal rule
  # spts is grid of x values, vls corresponding values of f(x)

  # vector withount first minus without vector last
  x_diff = tail(spts, -1) - head(spts, -1)
  y_diff = tail(vls, -1) + head(vls, -1)

  trapezes = x_diff * y_diff / 2
  sum(trapezes)
}


b_normal <- function (thx){
  thx^2 / 2
}

h_theoretical <- function (l0, l1, th0, th1, th){
  v = (log(l0) / (th - th0) + log(l1) / (th1 - th)) /
    (
      (b_normal(th1) - b_normal(th)) / (th1 - th) -
      (b_normal(th) - b_normal(th0)) / (th - th0)
    )
  return(floor(v))
}

bound <- function (l, tx0, tx1, n){
  log(l) / (tx1 - tx0) + n * (b_normal(tx1) - b_normal(tx0)) / (tx1 - tx0)
}


d <- function (s, thx, n) {
  # probability density function without common multiple
  # here the variable s can be a vector, other variables should be numbers

  v = exp(thx * s - (n * thx ^ 2) / 2 )
  return(v)
}


L <- function (s, l0, l1, th0, th1, H){
  # first (at step H - 1) integral value (found analiticaly)
  # here the variable s can be a vector, other variables should be numbers
  B = H * (th0 + th1) / 2 + log(l1 / l0) / (th0 - th1)
  rp1 = l0 * d(s, th0, H - 1) * (1 - erf((B - s - th0) / sqrt(2)))
  rp2 = l1 * d(s, th1, H - 1) * (1 + erf((B - s - th1) / sqrt(2)))
  return((rp1 + rp2) / 2)
}


calc_center_bound <- function (n, th0, th1){
  bound(l0/l1, th0, th1, n)
}

calc_left_bound <- function(n, l1, th1, th){
  # v = log(1 / l1) / (th1 - th) + n * (th1 + th) / 2
  # return(v)
  bound(1/l1, th, th1, n)
}

calc_right_bound <- function(n, l0, th0, th){
  # v = log(1 / l0) / (th0 - th) + n * (th0 + th) / 2
  # return(v)
  bound(l0, th0, th, n)
}


find_first_continuation_interval <- function (l0, l1, th0, th1, th, n){
  # exactly for step n

  current_loss <- function(s){ min(l0 * d(s,th0, n), l1 * d(s, th1, n)) }
  future_loss <- function(s){ d(s, th, n) + L(s, l0, l1, th0, th1, n + 1) }
  to_solve <- function (s){ current_loss(s) - future_loss(s) }

  c = calc_center_bound(n, th0, th1)  # center_bound
  a = calc_left_bound(n, l1, th1, th)  # left_bound
  b = calc_right_bound(n, l0, th0, th)  # right_bound

  if(sign(to_solve(c)) == sign(to_solve(b))) {
    cat("There is no continuation!")
    return(NA)
  }

  left_solution = uniroot(to_solve, c(a, c))
  right_solution = uniroot(to_solve, c(c, b))

  return(c(left_solution$root, right_solution$root))
}

calculate_integral <- function (s, stepdata, l0, l1, th0, th1){
  # s should be a number
  if (length(s) != 1) stop("s should be a number here")

  n = stepdata$n
  #x = stepdata$grid - s
  # left bound of interval
  a = head(stepdata$grid, 1)
  # right bound of interval
  b = tail(stepdata$grid, 1)

  #y = dnorm(x) * (stepdata$val)

  # numerical part of integral
  int3 = intgr(stepdata$grid,stepdata$val,s)
  # analytical parts of integral
  int1 = l0 * d(s, th0,  n - 1) * (1 - erf((b - s - th0) / sqrt(2))) / 2
  int2 = l1 * d(s, th1,  n - 1) * (1 + erf((a - s - th1) / sqrt(2))) / 2

  parts = c(int1, int2, int3)
  return(sum(parts))
}


find_continuation_interval <- function (stepdata, l0, l1, th0, th1, th){
  # all calculatuions exactly for n (previous step minus one)
  n = stepdata$n - 1

  current_loss<-function(s){ min(l0 * d(s, th0, n), l1 * d(s, th1, n)) }
  future_loss<-function(s){ d(s, th, n) + calculate_integral(s, stepdata, l0, l1, th0, th1) }
  to_solve <- function (s){ current_loss(s) - future_loss(s) }

  c = calc_center_bound(n, th0, th1)  # center_bound
  a = calc_left_bound(n, l1, th1, th)  # left_bound
  b = calc_right_bound(n, l0, th0, th)  # right_bound

  left_solution = uniroot(to_solve, c(a, c))
  right_solution = uniroot(to_solve, c(c, b))

  return(c(left_solution$root, right_solution$root))
}


step_n_has_continuation <- function (l0, l1, th0, th1, th, n){
  # checking exactly for step n

  current_loss <- function(s){ min(l0 * d(s,th0, n), l1 * d(s, th1, n)) }
  future_loss <- function(s){ d(s, th, n) + L(s, l0, l1, th0, th1, n + 1) }
  to_solve <- function (s){ current_loss(s) - future_loss(s) }

  c = calc_center_bound(n, th0, th1)  # center_bound
  b = calc_right_bound(n, l0, th0, th)  # right_bound

  if(sign(to_solve(c)) == sign(to_solve(b))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


modified_kw <- function(l0, l1, th0, th1, th, H=0, precision=0.005, print_output=TRUE){
  if (th <= th0 || th >= th1){
    stop("Only case th0 < th < th1 is supported!")
  }

  theor_h = h_theoretical(l0, l1, th0, th1, th)

  if (H > theor_h || H == 0){
    H = theor_h
  }

  test = list()

  n = H - 1
  repeat{
    if(step_n_has_continuation(l0, l1, th0, th1, th, n=n))
      break
    else {

      if (print_output) {
        print(n)
        print("There is no continuation!")
      }

      n = n - 1
    }
    if(n == 1)
      return(test)
  }

  acceptance_constant = bound(l0/l1, th0, th1, n + 1)
  test[[n + 1]] = list(n=n + 1, const=acceptance_constant, grid=acceptance_constant, val=NA)

  if (print_output) {
    print(n + 1)
    print(acceptance_constant)
  }


  continuation_interval = find_first_continuation_interval(l0, l1, th0, th1, th, n)
  a = continuation_interval[1]
  b = continuation_interval[2]
  nint = ceiling((b - a) / precision)
  grid = seq(a, b, length=nint + 1)
  val = d(grid, th, n) + L(grid, l0, l1, th0, th1, n + 1)
  stepdata = list(n=n, grid=grid, val=val)
  test[[n]] = stepdata

  if (print_output) {
    print(n)
    print(continuation_interval)
  }


  repeat{
    continuation_interval = find_continuation_interval(stepdata, l0, l1, th0, th1, th)
    n = stepdata$n - 1

    if (print_output) {
      print(n)
      print(continuation_interval)
    }

    a = continuation_interval[1]
    b = continuation_interval[2]
    nint = ceiling((b - a) / precision)
    grid = seq(a, b, length=nint + 1)
    val = d(grid, th, n) + Vectorize(calculate_integral, vectorize.args="s")(grid, stepdata, l0, l1, th0, th1)
    stepdata = list(n=n, grid=grid, val=val)
    test[[n]] = stepdata
    if(n == 1)
      break
  }
  test
}

calculate_g_values <- function (test, th){

  H = tail(test, n=1)[[1]]$n
  g_values = list()

  # Gn calculation
  grid = test[[1]]$grid
  val = dnorm(grid - th)
  stepdata = list(n=1, grid=grid, val=val)
  g_values[[1]] = stepdata

  for (n in 2:(H - 1)){
    calc_new_Gn <- function (y) intgr(g_values[[n - 1]]$grid, g_values[[n - 1]]$val, y - th)

    grid = test[[n]]$grid
    val = Vectorize(calc_new_Gn, vectorize.args="y")(grid)

    stepdata = list(n=n, grid=grid, val=val)
    g_values[[n]] = stepdata
  }

  g_values
}

operating_characteristic <- function (test, th, g_values=NULL){
  H = tail(test, n=1)[[1]]$n

  if (is.null(g_values)){
    g_values = calculate_g_values(test, th)
  }

  # Operating characteristic calculation
  oc = pnorm(test[[1]]$grid[1] - th)
  for (n in 1:(H - 1)){
    a = test[[n + 1]]$grid[1]
    grid = test[[n]]$grid
    to_int_val = g_values[[n]]$val * pnorm(a - grid - th)
    oc = oc + intgr_trapezoidal(grid, to_int_val)
  }

  oc
}

average_sample_number <- function (test, th, g_values=NULL){
  H = tail(test, n=1)[[1]]$n

  if (is.null(g_values)){
    g_values = calculate_g_values(test, th)
  }

  asn = 1
  for (n in 1:(H - 1)){
    asn = asn + intgr_trapezoidal(g_values[[n]]$grid, g_values[[n]]$val)
  }

  asn
}

sample_number_quantile <- function (test, th, q, g_values=NULL){
  H = tail(test, n=1)[[1]]$n

  if (is.null(g_values)){
    g_values = calculate_g_values(test, th)
  }

  n = 0
  p = 1
  if (p <= 1 - q) return(0)

  repeat{
    n = n + 1
    if (n == H) return(H)

    p = intgr_trapezoidal(g_values[[n]]$grid, g_values[[n]]$val)
    if (p <= 1 - q) return(n)
  }
}

