# The usage example

source("normal.R")



modified_kw_summary <- function(l0, l1, th0, th1, th, H, precision, print_output=TRUE) {
  # main function

  test = modified_kw(l0, l1, th0, th1, th, H, precision, print_output)
  maxN = length(test)

  g_values = calculate_g_values(test, th)
  g_values0 = calculate_g_values(test, th0)
  g_values1 = calculate_g_values(test, th1)

  ASN = average_sample_number(test, th, g_values)
  ASN0 = average_sample_number(test, th0, g_values0)
  ASN1 = average_sample_number(test, th1, g_values1)

  alpha_error = 1 - operating_characteristic(test, th0, g_values0)
  beta_error = operating_characteristic(test, th1, g_values1)


  SNQ = sample_number_quantile(test, th, 0.99, g_values)


  result_list = list(
    "lambda0" = l0,
    "lambda1" = l1,
    "theta0" = th0,
    "theta1" = th1,
    "theta" = th,
    "maxNumber" = maxN,
    "alpha" = alpha_error,
    "beta" = beta_error,
    # "Delta" = Delta,
    "ASN" = ASN,
    "ASN0" = ASN0,
    "ASN1" = ASN1,
    "Quantile" = SNQ
  )

  return (result_list)

}

if (sys.nframe() == 0){

  # input parameters can be set here
  l0 = 530.3
  l1 = 530.3
  th0 = -0.1
  th1 = 0.1
  th = 0
  H = 749
  precision = 0.01



  results = modified_kw_summary(l0, l1, th0, th1, th, H, precision)
  print(results)

  for (name in names(results)){
    print(paste(name, results[[name]]))

  }

}