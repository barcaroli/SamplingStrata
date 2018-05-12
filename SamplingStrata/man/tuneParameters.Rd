\name{tuneParameters}
\Rdversion{1.1}
\alias{tuneParameters}
\title{
Execution and compared evaluation of optimization runs
}
\description{
This function allows to execute a number of optimization runs,
varying in a controlled way the values of the parameters, in order to find their most suitable values.
by comparing the resulting solutions.
It can be applied only to a given domain per time.
Most parameters of this function are the same than those of the function 'optimizeStrata', but they are given 
in a vectorial format. The length of each vector is given by the number of optimizations to be run: 
it is therefore possible to define different combination of values of the parameters for each execution of 'optimizeStrata'.
After each optimization run, from the corrisponding optimized frame a given number of samples are drawn. For each of them, the 
estimates of the target variables Y's are computed ("precision"), together with the associated coefficients of variations, 
and the absolute differences between the values of the estimates and the true values in the population ("bias").
Information on the distribution of bias (differences) and precision (CV's) are outputted, and also boxplots 
for each of them are produced, in order to permit a compared evaluation of the different solutions found in the different runs.
As the optimal solution is stored for each run, after the evaluation it is possible to use it directly, or as a
"suggestion" for a new optimization with more iterations (in order to improve it).  
}
\usage{
tuneParameters (
		noptim,
		nsampl,
		frame,
		errors = errors, 
		strata = strata, 
		cens = NULL, 
		strcens = FALSE,
		alldomains = FALSE,
		dom = 1,			
		initialStrata, 
		addStrataFactor, 
		minnumstr, 
		iter, 
		pops, 
		mut_chance, 
		elitism_rate
		)
}
\arguments{
  \item{noptim}{
  Number of optimization runs to be performed
}
  \item{nsampl}{
  Number of samples to be drawn from the optimized population frame after each optimization
}
  \item{frame}{
  The (mandatory) dataframe containing the sampling frame
}
  \item{errors}{
  This is the (mandatory) dataframe containing the precision levels expressed in terms of Coefficients of Variation
  that estimates on target variables Y's of the survey must comply 
}
  \item{strata}{
  This is the (mandatory) dataframe containing the information related to "atomic" strata, i.e. the strata obtained by 
  the Cartesian product of all auxiliary variables X's. Information concerns the identifiability of strata
  (values of X's) and variability of Y's (for each Y, mean and standard error in strata) 
}
  \item{cens}{
  This the (optional) dataframe containing the takeall strata, those strata whose units must be selected in 
  whatever sample. It has same structure than "strata" dataframe
}
  \item{strcens}{
  Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE
}
  \item{alldomains}{
  Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains. It must be left to its default (FALSE) 
}
  \item{dom}{
  Indicates the domain on which the optimization runs must be performed. It is an integer value that has to be internal to the interval 
  (1 <--> number of domains). It is mandatory, if not indicated, the default (1) is taken.
}
  \item{initialStrata}{
  This is the initial limit on the number of strata for each solution. Default is 3000.
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}
  \item{addStrataFactor}{
  This parameter indicates the probability that at each mutation the number of strata may increase with 
  respect to the current value. Default is 0.01 (1%).
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}
  \item{minnumstr}{
  Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}
  \item{iter}{
  Indicated the maximum number of iterations (= generations) of the genetic algorithm. Default is 20.
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}
  \item{pops}{
  The dimension of each generations in terms of individuals. Default is 50.
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}
  \item{mut_chance}{
  Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
  of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
  but a slower convergence, while low values permit a faster convergence, but the final solution can be
  distant from the optimal one. Default is 0.05.
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}
  \item{elitism_rate}{
  This  parameter indicates the rate of better solutions that must be preserved from one generation
  to another. Default is 0.2 (20%).
  This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
  ( = value of parameter 'noptim')
}


}


\author{
Giulio Barcaroli
}

\examples{
#
\dontrun{
#------------------------------------------------------------
# data setting
library(SamplingStrata)
data(swissstrata)
data(swisserrors)
data(swissframe)
# As this function can be applied only to a given domain per time,
# we select the first domain
frame <- swissframe[swissframe$domainvalue == 1,]
strata <- swissstrata[swissstrata$DOM1 == 1,]
errors <- swisserrors[swisserrors$domainvalue == 1,]
#------------------------------------------------------------
# parameters setting
noptim <- 8 # Number of runs
nsampl <- 100 # Number of samples to be drawn after each optimization
initialStrata <- ceiling(c(1:noptim)*0.1*(nrow(strata))) # Number of initial strata
addStrataFactor <- rep(0.01,noptim) # Rate for increasing initial strata
minnumstr <- rep(2,noptim) # Minimum number of units per stratum
iter <- rep(200,noptim) # Number of iterations for each optimization
pops <- rep(20,noptim) # Number of solutions for each iteration
mut_chance <- rep(0.004,noptim) # Mutation chance
elitism_rate <- rep(0.2,noptim) # Elitism rate
#------------------------------------------------------------
tuneParameters (
  noptim,
  nsampl,
  frame,
  errors = errors, 
  strata = strata,
  cens = NULL, 
  strcens = FALSE,
  alldomains = FALSE,
  dom = 1,			
  initialStrata, 
  addStrataFactor, 
  minnumstr, 
  iter, 
  pops, 
  mut_chance, 
  elitism_rate
  )
}
}
\keyword{ survey }

