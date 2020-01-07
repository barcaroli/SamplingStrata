

<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

# SamplingStrata

This package offers an approach for the determination of the best
stratification of a sampling frame, the one that ensures the minimum
sample cost under the condition to satisfy precision constraints in a
multivariate and multidomain case. This approach is based on the use of
the genetic algorithm: each solution (i.e. a particular partition in
strata of the sampling frame) is considered as an individual in a
population; the fitness of all individuals is evaluated applying the
Bethel-Chromy algorithm to calculate the sampling size satisfying
precision constraints on the target estimates.

Functions in the package allows to:

  - support in the preparation of required input data;
  - execute the optimisation step;
  - analyse the obtained results of the optimisation step;
  - select a sample from the new frame accordingly to the best
    allocation.

Functions for the execution of the genetic algorithm are a modified
version of the functions in the ‘genalg’ package.

A complete illustration of all features and functions can be found at
the link:

<https://barcaroli.github.io/SamplingStrata/articles/SamplingStrata.html>

## Installation

You can install SamplingStrata from github with:

``` r
# install.packages("devtools")
devtools::install_github("barcaroli/SamplingStrata")
```

## Example

This is a basic example which shows you how to solve a common problem:

```` r
library(SamplingStrata)

# Data ----------------------------------------------------------------------
data(swisserrors)
swisserrors
#    DOM  CV1  CV2  CV3  CV4 domainvalue
# 1 DOM1 0.08 0.12 0.08 0.12           1
# 2 DOM1 0.08 0.12 0.08 0.12           2
# 3 DOM1 0.08 0.12 0.08 0.12           3
# 4 DOM1 0.08 0.12 0.08 0.12           4
# 5 DOM1 0.08 0.12 0.08 0.12           5
# 6 DOM1 0.08 0.12 0.08 0.12           6
# 7 DOM1 0.08 0.12 0.08 0.12           7

data(swissstrata)
head(swissstrata)
#        STRATO   N       M1        M2        M3       M4        S1       S2        S3       S4 cost cens DOM1 X1 X2 X3 X4 X5 X6
# 1 1*1*1*1*1*1 184 48.31522  49.40217  61.44022 28.40761 26.888529 28.57606 32.719652 14.67916    1    0    1  1  1  1  1  1  1
# 2 1*1*1*1*1*2   1 98.00000 106.00000 116.00000 43.00000  0.000000  0.00000  0.000000  0.00000    1    0    1  1  1  1  1  1  2
# 3 1*1*1*2*1*1   2 57.00000  64.00000  70.00000 50.00000  5.656854  0.00000  1.414214 21.21320    1    0    1  1  1  1  2  1  1
# 4 1*1*2*1*1*1  11 77.72727  81.18182  92.36364 47.00000 15.994317 19.61029 17.862098 11.67048    1    0    1  1  1  2  1  1  1
# 5 1*2*1*1*1*1   9 58.22222  61.55556  66.77778 36.22222 27.008229 21.50065 26.409173 16.43759    1    0    1  1  2  1  1  1  1
# 6 1*2*1*2*1*1   8 61.00000  68.00000  84.62500 58.37500 26.262412 20.82581 28.172618 28.38982    1    0    1  1  2  1  2  1  1

data(swissframe)
head(swissframe)
#   progr REG X1 X2 X3 X4 X5 X6         id    Y1     Y2     Y3    Y4 domainvalue
# 1     1   4 18  3  2  1  3  3     Zurich 57324 131422 108178 66349           4
# 2     2   1 17  1  1  1  3  2     Geneve 32429  60074  57063 28398           1
# 3     3   3 17  1  1  1  3  3      Basel 28161  50349  53734 34314           3
# 4     4   2 17  2  3  1  3  3       Bern 19399  44263  39397 25575           2
# 5     5   1 17  2  2  1  3  2   Lausanne 24291  44202  35421 21000           1
# 6     6   4 16  3  3  1  3  3 Winterthur 18942  28958  27696 14887           4


# Starting solution with kmeans clustering -------------------------------------
kmean <- KmeansSolution(swissstrata, swisserrors, maxclusters = 10)
# number of strata to be obtained in each domain in final solution  
nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
                 FUN=function(x) length(unique(x)))
nstrat
# 1 2 3 4 5 6 7 
# 9 9 7 6 7 7 8 

# Optimisation step ------------------------------------------------------------
solution <- optimStrata (
  method = "atomic",
  errors = swisserrors,
  framesamp = swissframe,
  nStrata = nstrat,
  suggestions = kmean)

# optimized sampling strata with allocated units -------------------------------
outstrata <- solution$aggr_strata

# sampling frame with optimized strata labels ----------------------------------
framenew <- solution$framenew

# evaluate the current solution ------------------------------------------------
eval <- evalSolution(framenew, outstrata)
eval$coeff_var
#      CV1    CV2    CV3    CV4  dom
# 1 0.0800 0.0789 0.0811 0.0817 DOM1
# 2 0.0767 0.0842 0.0780 0.0858 DOM2
# 3 0.0746 0.0753 0.0750 0.0779 DOM3
# 4 0.0785 0.0738 0.0777 0.0748 DOM4
# 5 0.0810 0.0806 0.0809 0.0801 DOM5
# 6 0.0781 0.0763 0.0782 0.0808 DOM6
# 7 0.0680 0.0719 0.0719 0.0808 DOM7

# Select a sample --------------------------------------------------------------
s <- selectSample(framenew, outstrata)
# *** Sample has been drawn successfully ***
#   99  units have been selected from  50  strata
# 
# ==> There have been  6  take-all strata 
# from which have been selected  7 units

head(s)
#   DOMAINVALUE STRATO      STRATUM PROGR REG X1 X2 X3 X4 X5 X6                   ID   Y1   Y2   Y3   Y4 LABEL   WEIGHTS        FPC
# 1           1      1  1*1*1*1*1*1  2132   1  1  1  1  1  1  1            Bremblens   89   90  138   43     1 97.500000 0.01025641
# 2           1      1  3*1*2*1*1*1  1227   1  3  1  2  1  1  1 Granges-pres-Marnand  331  289  350  176     1 97.500000 0.01025641
# 3           1     10 11*1*1*1*2*2   119   1 11  1  1  1  2  2        Ecublens (VD) 2320 3548 3309 1050    10  3.428571 0.29166667
# 4           1     10  7*3*2*2*2*1   362   1  7  3  2  2  2  1            Le Chenit  993 1054 1340  910    10  3.428571 0.29166667
# 5           1     10 12*1*1*1*2*2    70   1 12  1  1  1  2  2               Morges 2940 4388 4408 2418    10  3.428571 0.29166667
# 6           1     10  2*1*1*2*1*1  2041   1  2  1  1  2  1  1             Grimentz  121   96  125   62    10  3.428571 0.29166667```
````
