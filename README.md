
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
precision constraints on the target estimates. Functions in the package
allows to: (a) analyse the obtained results of the optimisation step;
(b) assign the new strata labels to the sampling frame; (c) select a
sample from the new frame accordingly to the best allocation. Functions
for the execution of the genetic algorithm are a modified version of the
functions in the ‘genalg’ package.

## Installation

You can install SamplingStrata from github with:

``` r
# install.packages("devtools")
devtools::install_github("barcaroli/SamplingStrata")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SamplingStrata)
data(swisserrors)
data(swissstrata)
solution <- optimizeStrata (
    errors = swisserrors,
    strata = swissstrata,
  showPlot = FALSE)
# update sampling strata with new strata labels
newstrata <- updateStrata(swissstrata, 
                          solution, 
                          writeFiles = FALSE)
# update sampling frame with new strata labels
data(swissframe)
framenew <- updateFrame(frame=swissframe,
                        newstrata=newstrata,
                        writeFile=FALSE)
samp <- selectSample(framenew,solution$aggr_strata,writeFiles=TRUE)
# evaluate the current solution
eval <- evalSolution(frame = framenew, 
                     outstrata =solution$aggr_strata, 
                     nsampl = 100, 
                     cens = NULL, 
                     writeFiles = FALSE)
eval$coeff_var
swisserrors
```
