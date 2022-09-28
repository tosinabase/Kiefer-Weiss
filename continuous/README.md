# An R project for numerical solution of the Kiefer-Weiss problem
## Continuous Koopman-Darmois families

This  directory contains the R code for sampling from a normal population.


A detailed discussion of the problem and the description of the algorithms 
is provided in 
**Novikov, A., Novikov, A., Farkhshatov, F. Numerical solution of Kiefer-Weiss problems when sampling from 
continuous exponential families, 2022, to appear**



## Content description
* The file [normal.R](normal.R) contains all the functions providing the user interface for the Kiefer-Weiss 
and related problems. List of functions can be seen below.

* The file [usage.R](usage.R) is a usage example of these functions.


## [normal.R](normal.R) functions list

**modified_kw** calculates the optimal test  for the modified Kiefer-Weiss problem

**operating_characteristic** calculates the operating characteristic curve of a sequential test

**average_sample_number** calculates the average sample number of a sequential test

**sample_number_quantile** calculates the quantile of a given level of the sample number distribution of a sequential test

**calculate_g_values**  calculates values for functions **operating_characteristic**, **average_sample_number**, **sample_number_quantile**

## The following functions are used for calculation process of **modified_kw**

**d**

**bound**

**h_theoretical**

**b_normal**

**step_n_has_continuation**

**find_first_continuation_interval**

**find_continuation_interval**

**L**

**calculate_integral**

**calc_center_bound**

**calc_left_bound**

**calc_right_bound**

**intgr_trapezoidal**

**intgr**











