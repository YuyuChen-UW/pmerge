# pmerge R package

This package contains several functions to merge p-values that are either independent or arbitrarily dependent from multiple hypothesis testing.

## Installation

## Methods to merge p-values
This package contains the following functions/methods to merge p-values:
* The Generalized Mean Merging Function (pmean): It contains the generalized mean methods for independent p-values (Chen et al. 2020) and arbitrarily dependent p-values (Vovk and Wang (2020) and Vovk et al. (2020)).
* The Order Statistics Merging Function (porder): It contains the order statistics merging method (Vovk et al. 2020) for arbitrarily dependent p-values.
* The Harmonic Mean Merging Function (pharmonic): It contains the harmonic mean method (Wilson 2019) for independent p-values, the harmonic mean method (Vovk and Wang 2020) and the harmonic* merging method (Vovk et al. 2020) for arbitrarily depen- dent p-values.
* The Simes Merging Function (pSimes): It contains the Simes method (Simes 1986) for independent p-values, the Hommel method (Hommel 1983) and the grid harmonic merg- ing method (Vovk et al. 2020) for arbitrarily dependent p-values.
* The Cauchy Merging Function (pCauchy): It contains the Cauchy combination methods for independent p-values (Liu and Xie 2020) and arbitrarily dependent p-values (Chen et al. 2020).
## References
Chen Y, Liu P, Tan KS, and Wang R (2020). “Trade-off between validity and efficiency of merging p-values under arbitrary dependence.” arXiv preprint arXiv:2007.12366.

Hommel G (1983). “Tests of the overall hypothesis for arbitrary dependence structures.” Biometri- cal Journal, 25(5), 423–430.

Liu Y, and Xie J (2020). “Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures.” Journal of the American Statistical Association, 115(529), 393–402.

Simes RJ (1986). “An improved Bonferroni procedure for multiple tests of significance.” Biometrika, 73(3), 751–754.

Vovk V, and Wang R (2020). “Combining p-values via averaging.” Biometrika, 107(4), 791–808.

Vovk V, Wang B, and Wang R (2020). “Admissible ways of merging p-values under arbitrary dependence.” arXiv preprint arXiv:2007.14208.

Wilson DJ (2019). “The harmonic mean p-value for combining dependent tests.” Proceedings of the National Academy of Sciences, 116(4), 1195–1200.
