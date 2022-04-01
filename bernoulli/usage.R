source("core_bernoulli.R")

calculate <- function(l0, l1, th0, th1) {
  # main function

  a = original_kw(H=1, option=1, l0, l1, th0, th1)
  th = a[1]
  Delta = a[2]
  H = Hbound(l0, l1, th0, th1, th)

  test = modified_kw(H, l0, l1, th0, th1, th) # modified_kw(H, l1, l2, w, th0, th1, th, thp)
  maxN = max_number(test)
  ASN = average_sample_number(test, th)
  ASN0 = average_sample_number(test, th0)
  ASN1 = average_sample_number(test, th1)


  SNQ = sample_number_quantile(test, th, 0.99)


  result_list = list(
    "lambda0" = l0,
    "lambda1" = l1,
    "theta0" = th0,
    "theta1" = th1,
    "theta" = th,
    "maxNumber" = maxN,
    "alpha" = 1 - operating_characteristic(test, th0),
    "beta" = operating_characteristic(test, th1),
    "Delta" = Delta,
    "ASN" = ASN,
    "ASN0" = ASN0,
    "ASN1" = ASN1,
    "Quantile" = SNQ
  )

  return (result_list)

}

l0=157.696751972207
l1=193.349705609267

th0=0.05
th1=0.15

results = calculate(l0, l1, th0, th1)
print(results)

for (name in names(results)){
  print(paste(name, results[[name]]))

}
