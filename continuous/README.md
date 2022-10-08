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
* The notebook [modified_kw_performance.ipynb](modified_kw_performance.ipynb) 
contains a demonstration of usage of the program code for solution of the modified Kiefer-Weiss problem 
when sampling from a normal distribution, visualization and related calculations.  


## [normal.R](normal.R) functions list

**modified_kw** calculates the optimal test for the modified Kiefer-Weiss problem.

Takes arguments `l0, l1, th0, th1, th, H=0, precision=0.01, print_output=TRUE`.

* `l0, l1, th0, th1, th` are input parameters of sequential test.  Currently, only the case `th0 < th < th1` is supported.

* Input variable `H` is the desired maximum number of steps the test can take. It  will be reduced if step `H - 1` has no continuation. 
In this case or if `H` was set to 0, it will be redefined as the biggest number `n` such that `n - 1` 
has continuation.  The final  adjusted `H` can be obtained  as `length(test)`.

* `precision` is the grid size for numerical integration.
* `print_output` (boolean parameter)  Print continuation intervals at each step of test evaluation.


Returns R list `test`. Let `H` be `length(test)`. For `n` from `1` to `H - 1` `test[[n]]` is an R list with the following fields 
* `n` - step number,
* `grid` - grid on the continuation interval at step `n`,
* `val`- stored values for  construction of optimal test, calculated at each point of `grid`.

For the case`n=H`
* `n` -  step number,
* `const` - the acceptance constant,
* `grid` - only one value which equals `const` field,
* `val`- equals `NA`.



**operating_characteristic** calculates the operating characteristic curve of a sequential test

**average_sample_number** calculates the average sample number of a sequential test

**sample_number_quantile** calculates the quantile of a given level of the sample number distribution of a sequential test

**calculate_g_values**  calculates values for functions **operating_characteristic**, **average_sample_number**, **sample_number_quantile**

**calc_delta**

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











