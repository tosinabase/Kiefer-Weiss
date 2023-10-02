# An R project for numerical solution of the Kiefer-Weiss problem

The Kiefer-Weiss problem consists in finding  sequential hypothesis tests minimizing the maximum average sample size 
among all the tests whose error probabilities of the first and second kind do not exceed some prescribed levels.

This repository contains accompanying R program code for three articles:

[*Novikov, A., Novikov, A., Farkhshatov, F. 2022. 
A computational approach to  the Kiefer-Weiss problem for sampling from a Bernoulli population, 
Sequential Analysis  41 (02),  
pages 198-219*](https://www.tandfonline.com/doi/full/10.1080/07474946.2022.2070212),
arXiv.org preprint [*arXiv:2110.04802 [stat.ME]*](https://arxiv.org/abs/2110.04802)


[*Novikov, A. and Farkhshatov, F. 2022. 
Design and performance evaluation in Kiefer-Weiss problems when sampling from discrete exponential families, 
Sequential Analysis 41 (04),
pages 417 – 434*](https://doi.org/10.1080/07474946.2022.2109673), 
arXiv.org preprint [*arXiv:2203.13957 [stat.ME]*](https://arxiv.org/abs/2203.13957)

and

[*Novikov, A., Novikov, A., Farkhshatov, F. 2023. 
Numerical solution of Kiefer-Weiss problems when sampling from continuous exponential families, 
Sequential Analysis 42 (02), 
pages 189-209*](https://www.tandfonline.com/doi/abs/10.1080/07474946.2023.2193602)


The accompanying code for the first one, covering sampling from a Bernoulli population is placed in [bernoulli](bernoulli) directory. 

The accompanying code for the second paper, covering sampling from binomial, Poisson and negative binomial (Pascal)
distributions is placed in [discrete](discrete) directory.

The accompanying code for the  third paper covers  sampling from a normal population with a known variance
and is placed in [continuous](continuous) directory.
