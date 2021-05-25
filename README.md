# R package pmerge :This package contains several functions/methods to merge p-values that are either independent or arbitrarily dependent in multiple hypothesis testing.


### Table of Contents
**[Installation](#installation)**<br>
**[Introduction](#introduction)**<br>
**[Methods](#methods)**<br>
**[Usage](#usage)**<br>
**[References](#references)**<br>
## Installation
Install this package from Github with 
```r
#install.packages("devtools")
devtools::install_github("YuyuChen-UW/pmerge")
```
## Introduction
Merging p-values from different sources in multiple hypothesis testing has been a long-standing issue in many scientific investigation procedures. Many classical methods are designed for combining p-values with cerntain dependence assumption (e.g., independence). However, the validity of the test (in the sense that the probability of making a Type-I error is below the significance level) cannot be guaranteed if the dependence assumption of p-values is not satisfied. Let's take the classic Fisher's method for independent p-values as an example. P-values are obtained from one-sided normal tests and their test statistcs jointly follow a equicorrelated Gaussian normal distribution.
![\rho](https://latex.codecogs.com/svg.latex?\Large&space;\rho) | 0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 
--- | --- | --- | --- |--- |--- |--- 
Probability of making a Type-one error | 301 | 283 | 290 | 286 | 289 | 285 




## Methods
This package contains the following functions/methods to merge p-values:
1. The Generalized Mean Merging Function (`pmean`): 
   - The generalized mean methods for independent p-values (Chen et al. 2020).
   - The generalized mean methods for arbitrarily dependent p-values (Vovk and Wang (2020) and Vovk et al. (2020)).
2. The Order Statistics Merging Function (`porder`): 
   - The order statistics merging method for arbitrarily dependent p-values (Vovk et al. 2020).
3. The Harmonic Mean Merging Function (`pharmonic`): 
   - The harmonic mean method for independent p-values (Wilson 2019).
   - The harmonic mean method for arbitrarily dependent p-values (Vovk and Wang 2020).
   - The harmonic* merging method for arbitrarily dependent p-values (Vovk et al. 2020).
4. The Simes Merging Function (`pSimes`): 
   - The Hommel method for arbitrarily dependent p-values (Hommel 1983).
   - The Simes method for independent p-values (Simes 1986).
   - The grid harmonic merging method for arbitrarily dependent p-values (Vovk et al. 2020).
5. The Cauchy Merging Function (`pCauchy`): 
   - The Cauchy combination methods for independent p-values (Liu and Xie 2020).
   - The Cauchy combination methods for arbitrarily dependent p-values (Chen et al. 2020).
## Usage
- The functions `pmean`, `porder`, `pharmonic` and `pSimes` return the merged p-value (reject if the merged p-value is smaller than the significance level).
- The function `pCauchy` returns the details of the hypothesis test including the test statistic and the threshold (reject if the test statistic is smaller than the threshold).
- Examples:
```r
library("pmerge")
# Generate 1000 p-values
P <- runif(1000)

# The generalized mean method for arbitrarily dependent p-values with exponent being 0
pmean(p = P, r = 0, dependence = "A")
# The generalized mean method for independent p-values with exponent being 0
pmean(p = P, r = 0, dependence = "I")

# The order statistic method for arbitrarily dependent p-values
porder(p = P, k=1)

# The Hommel method for arbitrarily dependent p-values
pSimes(p = P, method = "H")
# The Simes method for independent p-values
pSimes(p = P, method = "S")
# The grid harmonic method for arbitrarily dependent p-values
pSimes(p = P, method = "G")

# The harmonic mean method for arbitrarily dependent p-values
pharmonic(p = P, method = "H1")
# The harmonic* merging method for arbitrarily dependent p-values
pharmonic(p = P, method = "H2")
# The harmonic mean method for independent p-values
pharmonic(p = P, method = "H3")

# The Cauchy combination method for arbitrarily dependent p-values with significance level 0.1
pCauchy(p = P, method = "A",epi = 0.1)
# The Cauchy combination method for independent p-values with significance level 0.1
pCauchy(p = P, method = "I",epi = 0.1)
```


## References
Chen Y, Liu P, Tan KS, and Wang R (2020). “Trade-off between validity and efficiency of merging p-values under arbitrary dependence.” arXiv preprint arXiv:2007.12366.

Hommel G (1983). “Tests of the overall hypothesis for arbitrary dependence structures.” Biometrical Journal, 25(5), 423–430.

Liu Y, and Xie J (2020). “Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures.” Journal of the American Statistical Association, 115(529), 393–402.

Simes RJ (1986). “An improved Bonferroni procedure for multiple tests of significance.” Biometrika, 73(3), 751–754.

Vovk V, and Wang R (2020). “Combining p-values via averaging.” Biometrika, 107(4), 791–808.

Vovk V, Wang B, and Wang R (2020). “Admissible ways of merging p-values under arbitrary dependence.” arXiv preprint arXiv:2007.14208.

Wilson DJ (2019). “The harmonic mean p-value for combining dependent tests.” Proceedings of the National Academy of Sciences, 116(4), 1195–1200.
