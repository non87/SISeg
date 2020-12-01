
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SISeg

SISeg (Statistical Inference for Segregation indices) provides general
tools to estimate segregation indices from sample data. The package
natively calculates the segregation indices and provides confidence
intervals based on several inferential techniques.

SISeg natively supports 7 segregation indices (D, Gini, Atkinson, Mutual
information, Theil, V and Isolation). Moreover, the user can specify her
own indices and use SISeg functions to get statistical inference. At the
present stage of development, SISeg only supports segregation indices
for two groups, but support to multi-group indices will be added in the
near future.

<!-- badges: start -->

<!-- badges: end -->

<!-- ## Installation -->

<!-- You can install the released version of SISeg from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("SISeg") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

<!-- ``` r -->

<!-- # install.packages("devtools") -->

<!-- devtools::install_github("non87/SISeg") -->

<!-- ``` -->

## Basic Principles

Intuitively, segregation regards the different distributions of groups
(for example, racial groups) among different units (for example,
schools). We can represent the groups’ distributions in a matrix,
showing how many individuals from each groups are in any unit.
Segregation indices are functions that map a positive real matrix to a
positive real number, quantifying the amount of segregation in an
environment.

For SISeg, the environment whose segregation is to be quantified is
represented in the form of a 2 by k matrix, where there are two groups
and k units. If the matrix is the product of a sample from a larger
population, one can apply a segregation index to the sample matrix to
estimate the population segregation. Yet, it turns out this procedure is
unreliable and likely to produce a positively-biased estimate.

SISeg uses bootstrap and Bayesian modeling to provide statistical
inference for segregation indices on sample matrices. The methodologies
used are based on 3 simple models of the sampling procedure, named
statistical frameworks: independent groups, full multinomial, and
independent units. The choice of a framework is practically irrelevant
for mid-sized and large samples. However, all the methodologies used
work in slightly-different ways under the different frameworks.

## Example

Based on the discussion above, SISeg works on the assumption that
relevant data for segregation in an environment is organized in a
matrix. Ths example starts from there and show how to conduct
statistical inference:

``` r
library(SISeg)
set.seed(1234)
## Create artificial data drawing from a multinomial. 
## This procedure corresponds to the full multinomial framework
env1 <- rmultinom(1, 1000, c(0.1, 0.1, 0.2, 0.2, 0.05, 0.35))
## Put it into a matrix form
env1 <- matrix(env1, nrow = 2, byrow = TRUE)
## We can simply calculate an index on the data. D index:
print(d_ind(env1))
#> [1] 0.2033333

## But it is just as easy to create a full confidence interval
## Estimate the confidence interval for D and V on this data
env1_cis <- index_ci(env1, seg_index = c("D", "V"))
#> [1] "d"
#> [1] "v"
#> [1] "d"
#> [1] "v"
## This is the bootstrap studentized confidence interval for the data for D
print(env1_cis$TBoot$d)
#>     Lower     Upper  Estimate 
#> 0.1562018 0.2530208 0.2033333

## SISeg can be used to check the difference between two environments
## We create another environment
env2 <- rmultinom(1, 1000, c(0.05, 0.05, 0.3, 0.1, 0.25, 0.25))
## Put it into a matrix form
env2 <- matrix(env1, nrow = 2, byrow = TRUE)
## Let's check the difference in Theil index between environment 1 and 2
## Let's calculate the confidence intervals of the difference from samples
env_diff_cis <- index_difference_ci(env1, env2, seg_index = c("Theil"))
#> [1] "theil"
#> [1] "theil"
```
