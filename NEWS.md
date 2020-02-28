
<!-- NEWS.md is generated from NEWS.Rmd. Please edit NEWS.Rmd file -->

# SamplingStrata 1.5-1

## Major changes

  - As starting from R 4.0.0 ‘stringsAsFactors = FALSE’ becomes a
    default, in all calls to data.frame the parameter ‘stringAsFactors =
    TRUE’ has been indicated, in order to ensure the same results

  - Fixed some bugs related to handling take-all strata for ‘atomic’ and
    ‘spatial’ methods

# SamplingStrata 1.5

## Major changes

  - A new *optimStrata* function is available: this function is a
    wrapper that allows to execute the three different optimization
    functions: (i) optimizeStrata (method = *atomic*); (ii)
    optimizeStrata (method = *continuous*); (iii) optimizeStrata (method
    = *spatial*).

  - A new *optimizeStrataSpatial* function is available: this function
    optimizes the frame stratification taking into account also spatial
    correlation of frame units in a territorial context. As for
    *optimizeStrata2*, this function can be used only on continuous
    stratification variables.

  - A new *kMeansSolutionSpatial* function has been added: this function
    is the same than *KmeansSolution*, but operates in the case of
    optimization with only continuous stratification variables and has
    to be used only in conjunction with *optimizeStrataSpatial*.

  - A new *KmeansSolution2* function has been added: this function is
    the same than *KmeansSolution*, but operates in the case of
    optimization with only continuous stratification variables and has
    to be used only in conjunction with *optimizeStrata2*.

  - A new *prepareSuggestion* function has been added: this function
    operates on the result of *kMeanSolution2* or *kMeanSolutionSpatial*
    order to prepare the suggestions for the optimization with only
    continuous stratification variables.

  - A new *computeGamma* function has been added: this function allows
    to calculate a heteroscedasticity index, to be passed to
    optimization step as a *model* parameter in order to correctly
    evaluate the variance in the strata.

# SamplingStrata 1.4-1

## Major changes

  - A new *summaryStrata* function has been added, enabling to get
    structured information regarding the strata produced by the
    *optimizeStrata2* function (operating on only continuous
    stratification variables).

  - A new *assignStrataLabel* function has been added, enabling to
    assign the optimized strata label to new units to be added in the
    sampling frame.

  - Fixed a bug in the optimization step for continuous variables

# SamplingStrata 1.4

## Major changes

  - A new function *optimizeStrata2* is available. This function
    performs the same task than *optimizeStrata*, but with a different
    Genetic Algorithm, operating on real values genome, instead of an
    integer one. This pemits to operate directly on the boundaries of
    the strata, instead of aggregating the initial atomic strata. In
    some situations (limited size of sampling frame) this new function
    is much more efficient. The limitation is in the nature of the
    stratification variables, that are required to be all continuous
    (though categorical ordinal could be handled).

  - A new function *expected\_CV* has been added to calculate CVs on
    target variables in different domains that may be expected from a
    given solution, output of the *optimizeStrata* execution.

  - Fixed a bug in the optimization step when considering also take-all
    strata

# SamplingStrata 1.3

## New functions

  - The optimization of population frame is run in parallel if different
    domains are considered. To this end, the parameter *parallel* can be
    set to TRUE in the *optimizeStrata* function. If not specified, n-1
    of total available cores are used OR if number of domains \< (n-1)
    cores, then number of cores equal to number of domains are used.

  - A new function *KmeansSolution* produces an initial solution using
    the kmeans algorithm by clustering atomic strata considering the
    values of the means of target variables in them. Also, if the
    parameter *nstrata* is not indicated, the optimal number of clusters
    is determined inside each domain, and the overall solution is
    obtained by concatenating optimal clusters obtained in domains. By
    indicating this solution as a suggestion to the optimization step,
    this may greatly speed the convergence to the optimal solution.

  - A new function *selectSampleSystematic* has been added. It allows to
    select a stratified sample with the systematic method, that is a
    selection that begins selecting the first unit by an initial
    randomly chosen starting point, and proceeding in selecting other
    units by adding an interval that is the inverse of the sampling rate
    in the stratum. This selection method can be useful if associated to
    a particular ordering of the selection frame, where the ordering
    variable(s) can be considered as additional stratum variable(s).

  - It is now possible to handle *anticipated variance* by introducing a
    model linking a proxy variable whose values are available for all
    units in the sampling frame, with the target variable whose values
    are not available. In this implementation only linear model can be
    chosen. When calling the *buildStrataDF* function, a dataframe is
    given, containing three parameters (beta, sigma2 and gamma) for each
    couple target / proxy. On the basis of these parameters, means and
    standard deviations in sampling strata are calculated accordingly to
    given formulas.

# SamplingStrata 1.2

## Major changes

  - The crossover function in the genetic algorithm has been modified by
    considering the *grouping* version of this algorithm: instead of
    mixing chromosomes in an indifferentiate way, groups of them in one
    parent (representing already aggregated strata) are attributed to
    the other parent when generating a child, preserving their
    composition. Moreover, parents are selected with a probability
    proportional to their fitness (in the previous version selection was
    completely at random). In many cases this can greatly speed the
    convergence to an optimal solution.

  - A new function *adjustSize* has been added. It allows to adjust the
    sample size and related allocation in strata on the basis of an
    externally indicated overall sample size. The adjustment of the
    sample size is perfomed by increasing or decreasing it
    proportionally in each optimized stratum.

  - A new function *buildFrameDF* has been added. It allows to create a
    *sampling frame* dataframe by indicating the dataset in which the
    information on all the units are contained, the identifier, the X
    variables, the Y variable and the variable that indicates the
    domains of interest.

  - In function *optimizeStrata* now the *initialStrata* parameter is a
    vector, whose length is equal to the number of strata in the
    different domains.

  - All outputs are written to an .subdirectory.

# SamplingStrata 1.1

## Major changes

  - Function *memoise* from the same package (now required) is applied
    before each evaluation in order to save processing time. This may
    largely increase the efficiency of the algorithm.

  - A *recode* function is applied on every generated solution in order
    to recode a genotype of n genes with k\<=n distinct alleles 1, 2, …,
    k in such a way that the distinct alleles of the recoded genotype
    appear in the natural order 1, 2, …, k. This avoids to consider as
    distinct two solutions that are equivalent but make use of a
    different coding.

# SamplingStrata 1.0-4

## Bug fix for old releases

# SamplingStrata 1.0-3

## New functions and changes

  - Modified the output to the console or to the file of results: of all
    the solutions, only the optimal value for each generation is
    visualised.

  - Now the visualisation of the trend in the optimal and mean values is
    optional: the plot can be avoided by setting showPlot = FALSE when
    calling optimizeStrata. It can be advisable when the number of
    iterations is very high.

# SamplingStrata 1.0-2

## Bug fix for old releases.

# SamplingStrata 1.0-1

## Major changes

  - The object returned by function *optimizeStrata* is no more a
    dataframe but a list: (i) the first element of the list is the
    solution vector
    (solution\(indices); (ii) the second element of the list is the dataframe containing aggregated strata (solution\)aggr\_strata).

  - In all the functions that previously produced .csv files and .pdf
    plots in the working directory, as a default this is no more the
    current behaviour. To write these files, it is now necessary to set
    the *writeFiles* flag to TRUE.

# SamplingStrata 1.0

## Major changes

  - Two new functions: (i) *evalSolution*, to evaluate the found
    solution in terms of expected target variables precision and bias
    obtainable by samples drawn from the otpimized frame; (ii)
    *tuneParameters*, to determine the best combination of values to
    assign to the parameters necessary for the execution of the genetic
    algorithm used for the optimization of the frame stratification.
