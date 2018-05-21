\name{KmeansSolution}
\Rdversion{1.1}
\alias{KmeansSolution}
\title{
Initial solution obtained by applying kmeans clustering of atomic strata
}
\description{
In order to speed the convergence towards the optimal solution, an initial one can be 
given as "suggestion" to "optimizeStrata" function.
The function "KmeansSolution" produces this initial solution using the kemans 
algorithm by clustering atomic 
strata considering the values of the means of target variables in them.
Also, if the parameter "nstrata" is not indicated, the optimal number of clusters 
is determined inside each domain, and the overall solution is obtained by concatenating optimal 
cluster obtained in domains.
The result is a dataframe with two columns: the first indicates the clusters, 
the second the domains.
}
\usage{
KmeansSolution(strata, errors, nstrata=NA, showPlot=FALSE)
}
\arguments{
  \item{strata}{
  The (mandatory) dataframe containing the information related to atomic strata.
}
  \item{errors}{
  The (mandatory) dataframe containing the precision constraints on target variables.
}
  \item{nstrata}{
  Number of aggregate strata (if NULL, it is optimized by varying the number of cluster
  from 2 to half number of atomic strata). Default is NA.
}
  \item{showPlot}{
Allows to visualise on a plot the different sample sizes for each number of aggregate strata. Default is FALSE.
}
}

\value{
A dataframe containing the solution
}

\author{
Giulio Barcaroli 
}

 
\examples{
#
# The following example is realistic, but is time consuming
#
\dontrun{
library(SamplingStrata)
data(swisserrors)
data(swissstrata)

# suggestion
solutionKmean <- KmeansSolution(strata=swissstrata,
								  errors=swisserrors,
								  nstrata=NULL)

# optimisation of sampling strata
solution <- optimizeStrata ( 
  errors = swisserrors, 
  strata = swissstrata, 
  cens = NULL,
  strcens = FALSE,
  initialStrata = as.numeric(table(swissstrata$DOM1)),
  addStrataFactor = 0.0,    
  minnumstr = 2,
  iter = 50,
  pops = 20,
  mut_chance = NA,
  elitism_rate = 0.2,
  highvalue = 100000000,
  suggestions = solutionKmean,
  realAllocation = TRUE,
  writeFiles = FALSE,
  showPlot = TRUE)
}
}
\keyword{ survey }

