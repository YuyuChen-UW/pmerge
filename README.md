# pmerge R package

This package contains several functions/methods to merge p-values that are either independent or arbitrarily dependent from multiple hypothesis testing.

## Installation
Install this package from Github with 
```
#install.packages("devtools")
devtools::install_github("YuyuChen-UW/pmerge")
library("pmerge")

```
## Methods to merge p-values
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
   - The Simes method for independent p-values (Simes 1986).
   - The Hommel method for arbitrarily dependent p-values (Hommel 1983).
   - The grid harmonic merging method for arbitrarily dependent p-values (Vovk et al. 2020).
5. The Cauchy Merging Function (`pCauchy`): 
   - The Cauchy combination methods for independent p-values (Liu and Xie 2020).
   - The Cauchy combination methods for arbitrarily dependent p-values (Chen et al. 2020).
## Usage
- The functions `pmean`, `porder`, `pharmonic` and `pSimes` return the merged p-value, e.g.,
``
# Generate 1000 p-values
P <- runif(1000)
# Merge p-values using the harmonic mean method for arbitrarily dependent p-values 
``
- The function `pCauchy` returns the details of the hypothesis test including the test statistic and the threshold (reject if the test statistic is smaller than the threshold).
## References
Chen Y, Liu P, Tan KS, and Wang R (2020). “Trade-off between validity and efficiency of merging p-values under arbitrary dependence.” arXiv preprint arXiv:2007.12366.

Hommel G (1983). “Tests of the overall hypothesis for arbitrary dependence structures.” Biometri- cal Journal, 25(5), 423–430.

Liu Y, and Xie J (2020). “Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures.” Journal of the American Statistical Association, 115(529), 393–402.

Simes RJ (1986). “An improved Bonferroni procedure for multiple tests of significance.” Biometrika, 73(3), 751–754.

Vovk V, and Wang R (2020). “Combining p-values via averaging.” Biometrika, 107(4), 791–808.

Vovk V, Wang B, and Wang R (2020). “Admissible ways of merging p-values under arbitrary dependence.” arXiv preprint arXiv:2007.14208.

Wilson DJ (2019). “The harmonic mean p-value for combining dependent tests.” Proceedings of the National Academy of Sciences, 116(4), 1195–1200.
