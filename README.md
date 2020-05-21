
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/SamplingStrata)](https://cran.r-project.org/package=SamplingStrata)
[![Downloads](http://cranlogs.r-pkg.org/badges/SamplingStrata?color=brightgreen)](http://www.r-pkg.org/pkg/SamplingStrata)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)

# SamplingStrata <img src="pkgdown/favicon/apple-touch-icon-152x152.png" align="right" />

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

<https://barcaroli.github.io/SamplingStrata/>

## Installation

You can install SamplingStrata from github with:

``` r
install.packages("devtools")
devtools::install_github("barcaroli/SamplingStrata")
```

<img src="cheat_sheet_page1.png" /> <img src="cheat_sheet_page2.png" />

## Three different methods for the optimization step

The optimization can be run by indicating three different methods, on
the basis of the following: A. if stratification variables are
categorical (or reduced to) then the method is the “atomic”; B. if
stratification variables are continuous, then the method is the
“continuous”; C. if stratification variables are continuous, and there
is spatial correlation among units in the sampling frame, then the
required method is the “spatial”.

## Example with the “continuous” method

``` r
library(SamplingStrata)

# Load data ---------------------------------------------------------------------------------
data("swissmunicipalities")
head(swissmunicipalities[,c(2:6,9,22)])
#   REG  COM        Nom HApoly Surfacesbois Airbat POPTOT
# 1   4  261     Zurich   8781         2326   2884 363273
# 2   1 6621     Geneve   1593           67    773 177964
# 3   3 2701      Basel   2391           97   1023 166558
# 4   2  351       Bern   5162         1726   1070 128634
# 5   1 5586   Lausanne   4136         1635    856 124914
# 6   4  230 Winterthur   6787         2807    972  90483

# Define the sampling frame -----------------------------------------------------------------
frame <-buildFrameDF(df= swissmunicipalities,
                     id = "COM",                    # unique identifier of sampling units
                     domainvalue= "REG",            # domain variable (region)
                     X = c("POPTOT","HApoly"),      # stratification variables
                     Y =c("Surfacesbois","Airbat")) # target variables
head(frame)
#     id     X1   X2   Y1   Y2 domainvalue
# 1  261 363273 8781 2326 2884           4
# 2 6621 177964 1593   67  773           1
# 3 2701 166558 2391   97 1023           3
# 4  351 128634 5162 1726 1070           2
# 5 5586 124914 4136 1635  856           1
# 6  230  90483 6787 2807  972           4

# Define precision constraints ------------------------------------------------------------
ndom <- length(unique(frame$domainvalue))
cv <- as.data.frame(list(DOM = rep("DOM1",ndom),
                         CV1 = rep(0.10,ndom),      # precision (cv=10%) for 'Surfacesbois'
                         CV2 = rep(0.10,ndom),      # precision (cv=10%) for 'Airind'
                         domainvalue= c(1:ndom)))   # same precision constraints for all domains
cv
#    DOM CV1 CV2 domainvalue
# 1 DOM1 0.1 0.1           1
# 2 DOM1 0.1 0.1           2
# 3 DOM1 0.1 0.1           3
# 4 DOM1 0.1 0.1           4
# 5 DOM1 0.1 0.1           5
# 6 DOM1 0.1 0.1           6
# 7 DOM1 0.1 0.1           7


# Find an initial solution and a suitable number of final strata in each domain -----------
solutionKmean <- KmeansSolution2(frame = frame,      # sampling frame
                                 errors = cv,        # precision constraints
                                 maxclusters = 10)   # max number of strata to be evaluated 
# number of strata to be obtained in each domain in final solution:                             
nstrat <- tapply(solutionKmean$suggestions, solutionKmean$domainvalue,
                 FUN=function(x) length(unique(x)))
nstrat
# 1  2  3  4  5  6  7 
# 10 10  7 10 10  8 10 

# Optimization step ------------------------------------------------------------------------
solution <- optimStrata(method = "continuous",        # method
                        framesamp = frame,            # sampling frame
                        errors = cv,                  # precision constraints
                        nStrata = nstrat,             # strata to be obtained in the final stratification
                        iter = 50,                    # number of iterations
                        pops = 10)                    # number of stratifications evaluated at each iteration
# Input data have been checked and are compliant with requirements
# 
# *** Starting parallel optimization for  7  domains using  5  cores
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=19s  
# 
# *** Sample size :  209
# *** Number of strata :  55
# ---------------------------

strataStructure <- summaryStrata(solution$framenew, 
                                 solution$aggr_strata)
head(strataStructure)
# Domain Stratum Population Allocation SamplingRate Lower_X1 Upper_X1 Lower_X2 Upper_X2
# 1      1       1        318          9     0.027928       27     1048       32     1166
# 2      1       2         91          7     0.077012       99     7516       48     1202
# 3      1       3         82          8     0.091802       95     8114      436     2635
# 4      1       4         26          3     0.125069      286    29559      219     2792
# 5      1       5         49          9     0.176848      130    22454     2864     7093
# 6      1       6          8          2     0.250000       78     6261     8460     9082

# Sample selection --------------------------------------------------------------------------
s <- selectSample(frame = solution$framenew,        # frame with the indication of optimized strata
                  outstrata = solution$aggr_strata) # optimized strata with sampling units allocation 
# *** Sample has been drawn successfully ***
#   209  units have been selected from  55  strata
# 
# ==> There have been  8  take-all strata 
# from which have been selected  11 units
head(s)
# DOMAINVALUE STRATO   ID  X1  X2  Y1 Y2 LABEL  WEIGHTS        FPC
# 1           1      1 5432 439 497 230 11     1 35.33333 0.02830189
# 2           1      1 5928  87 235  47  5     1 35.33333 0.02830189
# 3           1      1 5785 412 794 335 20     1 35.33333 0.02830189
# 4           1      1 5608 300 137  13  9     1 35.33333 0.02830189
# 5           1      1 5652 560 306  54  8     1 35.33333 0.02830189
# 6           1      1 5462 339 237  55 19     1 35.33333 0.02830189
```
