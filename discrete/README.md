# An R project for numerical solution of the Kiefer-Weiss problem
## Discrete Koopman-Darmois families

This  directory contains the R code for discrete exponential (Koopman-Darmois) families. For the time being, only binomial, Poisson and 
negative binomial (Pascal) distributions are covered.


A detailed discussion of the problem and the description of the algorithms is provided in
[*Novikov, A. and Farkhshatov, F. Design and performance evaluation in Kiefer-Weiss problems when sampling from discrete 
exponential families, 2022, arXiv.org preprint arXiv:2203.13957 [stat.ME]*](https://arxiv.org/abs/2203.13957)


## Content description
* The file [KW_common.R](KW_common.R) contains all the functions providing the common user interface for the Kiefer-Weiss 
and related problems. List of functions can be seen below. 
* The distribution-specific functions are placed in:

    * [KW_binom_set.R](KW_binom_set.R)  for the binomial distribution with success probability theta. 
  The number of trials is a global variable **size** (default value 1). 

    * [KW_pois_set.R](KW_pois_set.R) for the Poisson distribution with mean value theta.

    * [KW_nbinom_set](KW_nbinom_set) for the negative binomial (Pascal) distribution with mean **size***theta, 
  where **size** is a global 
variable with default value equal to 1 (which corresponds to the geometric distribution with mean theta).

* The file [usage.R](usage.R) is a usage example of these functions.


## [KW_common.R](KW_common.R) functions list

**modified_kw** calculates the optimal test  for the modified Kiefer-Weiss problem

**original_kw** calculates the "minimax" point \theta and  \Delta for construction of optimal test in the original
Kiefer-Weiss problem (see Proposition 2 in [https://arxiv.org/abs/2110.04802](https://arxiv.org/abs/2110.04802))

**operating_characteristic** calculates the operating characteristic curve of a sequential test

**average_sample_number** calculates the average sample number of a sequential test

**Lagr** calculates the Lagrangian function for the modified Kiefer-Weiss problem

**monte_carlo_simulation** performs Monte Carlo simulation for a sequential test. Returns estimated values of operating 
characteristic, average sample number and standard deviation of the sample number of the test.

**prob_to_stop_after** calculates the probability that a sequential test stops (strictly) after a given number of steps

**sample_number_quantile** calculates the quantile of a given level of the sample number distribution of a sequential test

**np_test** calculates the most powerful (Neyman-Pearson) test based on a given number of observations. Returns the 
critical constant, the randomization constant and the error probability of the second kind

**np** calculates the minimum sample number for the (one-sample) Neyman-Pearson test providing given error probabilities 
of the first and second kind (see Subsection 3.5 in [https://arxiv.org/abs/2203.13957](https://arxiv.org/abs/2203.13957) 
for details).





