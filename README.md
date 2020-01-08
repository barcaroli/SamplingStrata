
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
install.packages("devtools")
devtools::install_github("barcaroli/SamplingStrata")
```

## Three different methods for the optimization step

The optimization can be run by indicating three different methods, on
the basis of the following: A. if stratification variables are
categorical (or reduced to) then the method is the “atomic”; B. if
stratification variables are continuous, then the method is the
“continuous”; C. if stratification variables are continuous, and there
is spatial correlation among units in the sampling frame, then the
required method is the “spatial”.

## Example with the “atomic” method

``` r
library(SamplingStrata)

# Load data ---------------------------------------------------------------------------------
data("swissmunicipalities")
head(swissmunicipalities[,c(2:7,10,22)])
#   REG  COM        Nom HApoly Surfacesbois Surfacescult Airind POPTOT
# 1   4  261     Zurich   8781         2326          967    260 363273
# 2   1 6621     Geneve   1593           67           31     60 177964
# 3   3 2701      Basel   2391           97           93    213 166558
# 4   2  351       Bern   5162         1726         1041    212 128634
# 5   1 5586   Lausanne   4136         1635          714     64 124914
# 6   4  230 Winterthur   6787         2807         1827    238  90483

# Define the sampling frame -----------------------------------------------------------------
frame <-buildFrameDF(df= swissmunicipalities,
                     id = "COM",                    # unique identifier of sampling units
                     domainvalue= "REG",            # domain variable (region)
                     X = c("POPTOT","HApoly"),      # stratification variables
                     Y =c("Surfacesbois","Airind")) # target variables
head(frame)
#     id     X1   X2   Y1  Y2 domainvalue
# 1  261 363273 8781 2326 260           4
# 2 6621 177964 1593   67  60           1
# 3 2701 166558 2391   97 213           3
# 4  351 128634 5162 1726 212           2
# 5 5586 124914 4136 1635  64           1
# 6  230  90483 6787 2807 238           4

# Define precision constraints ------------------------------------------------------------
ndom <- length(unique(frame$domainvalue))
cv <- as.data.frame(list(DOM = rep("DOM1",ndom),
                         CV1 = rep(0.10,ndom),      # precision (cv=10\%) for 'Surfacesbois'
                         CV2 = rep(0.10,ndom),      # precision (cv=10\%) for 'Airind'
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

# Build atomic strata ---------------------------------------------------------------------
strata <- buildStrataDF(frame)
# Number of strata:  2895
# ... of which with only one unit:  2894> head(strata)
head(strata)
#              STRATO N   M1 M2 S1 S2 COST CENS DOM1    X1   X2
# 100*305     100*305 1   59  0  0  0    1    0    1   100  305
# 1010*1661 1010*1661 1  983  0  0  0    1    0    1  1010 1661
# 102*306     102*306 1   65  0  0  0    1    0    1   102  306
# 1020*5351 1020*5351 1 1375  2  0  0    1    0    1  1020 5351
# 10227*571 10227*571 1   73 48  0  0    1    0    1 10227  571
# 10230*330 10230*330 1   15  2  0  0    1    0    1 10230  330

# Find an initial solution and a suitable number of final strata in each domain -----------
solutionKmean <- KmeansSolution(strata = strata,    # atomic strata
                                errors = cv,        # precision constraints
                                maxclusters = 10)   # max number of strata to be evaluated 
# number of strata to be obtained in each domain in final solution:                             
nstrat <- tapply(solutionKmean$suggestions, solutionKmean$domainvalue,
                 FUN=function(x) length(unique(x)))
nstrat
# 1  2  3  4  5  6  7 
# 9 10  9  8 10  8  7 

# Optimization step ------------------------------------------------------------------------
solution <- optimStrata(method = "atomic",            # method
                        framesamp = frame,            # sampling frame
                        errors = cv,                  # precision constraints
                        nStrata = nstrat,             # strata to be obtained in the final stratification
                        suggestions = solutionKmean,  # initial solution
                        iter = 50,                    # number of iterations
                        pops = 10)                    # number of stratifications evaluated at each iteration
# Number of strata:  2895
# ... of which with only one unit:  2894
#  *** Starting parallel optimization for  7  domains using  5  cores
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=20s  
# 
# *** Sample size :  990
# *** Number of strata :  53

head(solution$aggr_strata)
# STRATO       M1        M2       S1        S2   N DOM1 COST CENS      SOLUZ
# 1      1 143.7857 29.857143 152.8819 45.077484  14    1    1    0  14.000000
# 2      2 379.6769  5.600000 734.8116 15.324490  65    1    1    0  34.089876
# 3      3 413.6063  6.196850 546.5235 15.803687 254    1    1    0 137.360990
# 4      4 574.6667  7.190476 454.4968 11.639911  21    1    1    0   8.364819
# 5      5 403.2778  1.944444 843.3638  3.135322  18    1    1    0   2.000000
# 6      6 358.1860  3.412791 716.8800  7.659047 172    1    1    0  45.117978

# Sample selection --------------------------------------------------------------------------
s <- selectSample(frame = solution$framenew,        # frame with the indication of optimized strata
                  outstrata = solution$aggr_strata) # optimized strata with sampling units allocation 
# *** Sample has been drawn successfully ***
#   990  units have been selected from  53  strata
# 
# ==> There have been  6  take-all strata 
# from which have been selected  60 units
head(s)
# DOMAINVALUE STRATO    STRATUM   ID    X1   X2  Y1  Y2 LABEL WEIGHTS FPC
# 1           1      1    277*495 5669   277  495  74   0     1       1   1
# 2           1      1    322*269 5858   322  269  34   2     1       1   1
# 3           1      1 27171*2562 6266 27171 2562 277 126     1       1   1
# 4           1      1    151*345 5674   151  345  58   0     1       1   1
# 5           1      1    288*312 5513   288  312  56   4     1       1   1
# 6           1      1    172*193 5801   172  193  14   1     1       1   1
```
