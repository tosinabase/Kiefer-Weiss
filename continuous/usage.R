# The usage example

source("normal.R")



modified_kw_summary <- function(test, l0, l1, th0, th1, th, to_calc_delta=FALSE) {
  H = test[["info"]]$H
  g_values = calculate_g_values(test, th)
  g_values0 = calculate_g_values(test, th0)
  g_values1 = calculate_g_values(test, th1)

  result_list = list(
    "lambda0" = l0,
    "lambda1" = l1,
    "theta0" = th0,
    "theta1" = th1,
    "theta" = th,
    "H" = H,
    "alpha" = 1 - operating_characteristic(test, th0, g_values0),
    "beta" = operating_characteristic(test, th1, g_values1),
    "ASN" = average_sample_number(test, th, g_values),
    "ASN0" = average_sample_number(test, th0, g_values0),
    "ASN1" = average_sample_number(test, th1, g_values1),
    "quantile" = sample_number_quantile(test, th, 0.99, g_values),
    "acceptance_constant" = test$data[[H]]$const
  )

  if (to_calc_delta)
    result_list[["Delta"]] = calc_delta(test, c(th0, th1), asn_th=result_list[["ASN"]])

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


  test = modified_kw(l0, l1, th0, th1, th, H, precision)
  results = modified_kw_summary(test, l0, l1, th0, th1, th)
  print(results)

  for (name in names(results)){
    print(paste(name, results[[name]]))

  }

}