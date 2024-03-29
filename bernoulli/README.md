# An R project for numerical solution of the Kiefer-Weiss problem
## R code to provide a numerical solution to the Kiefer-Weiss problem for sampling from Bernoulli populations


The Kiefer-Weiss problem consists in finding a sequential hypothesis tests minimizing the maximum average sample size 
among all the tests whose error probabilities of the first and second kind do not exceed some prescribed levels.

A detailed discussion of the problem and the description of the algorithms is provided in
[*Novikov, A., Novikov, A., Farkhshatov, F. A computational approach to  the Kiefer-Weiss problem for sampling from 
a Bernoulli population, Sequential Analysis - Volume 41, 2022 - Issue 2, Pages 198-219*](https://www.tandfonline.com/doi/full/10.1080/07474946.2022.2070212)
[(arXiv.org preprint)](https://arxiv.org/abs/2110.04802)

At the time being, only sampling from a Bernoulli population is covered here.

## Content description
* The file [core_bernoulli.R](core_bernoulli.R) contains all functions used for calculation. List of functions can be seen below. 
* The file [usage.R](usage.R) is a usage example of these functions.
* The file [results.csv](results.csv) contains numerical results of these scripts. 
You can try to use `lambda0, lambda1, th0, th1` variables as input in [usage.R](usage.R) script to verify the data in [results.csv](results.csv).
* The directory [graphs](graphs) contains interactive versions of graphs which were prepared using the 3D visualization [rgl package](https://github.com/dmurdoch/rgl). 
To see these plots you have to download html file and run on your browser.


## [core_bernoulli.R](core_bernoulli.R) functions list

* **modified_kw**: calculates the optimal  test for the modified Kiefer-Weiss problem
* **original_kw**: calculates  the optimal test for the original Kiefer-Weiss problem
* **operating_characteristic**: calculates the operating characteristic function for a test
* **average_sample_number**: calculates the average sample number of a test
* **prob_not_to_stop_before**: calculates the probability that a test does not stop before a given time k
* **sample_number_quantile**: calculates the quantile  of the distribution of the sample number
* **monte_carlo_simulation**: calculates, by Monte Carlo simulation, the operating characteristic, the average sample number, 
and  the standard deviation of the  sample number of a test
