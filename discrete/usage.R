# The usage example

#  One of the following options should be loaded here depending on which distribution will be used for the calculation:
#  KW_binom_set.R, KW_nbinom_set.R, KW_pois_set.R, etc...

source("KW_binom_set.R")
# source("KW_nbinom_set.R")
# source("KW_pois_set.R")



calculate <- function(l0, l1, th0, th1) {
  # main function

  a = original_kw(l0, l1, th0, th1)
  th = a[1]
  Delta = a[2]

  test = modified_kw(l0, l1, th0, th1, th)
  maxN = length(test)
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


# input parameters can be set here
l0 = 449.995105192065
l1 = 489.748951643705

th0 = 0.05
th1 = 0.08


# size is the global variable which definded in KW_binom_set.R, KW_nbinom_set.R. Can be changed here.
size = 3

results = calculate(l0, l1, th0, th1)
print(results)

for (name in names(results)){
  print(paste(name, results[[name]]))

}
