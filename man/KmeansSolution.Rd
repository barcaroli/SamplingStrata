\name{KmeansSolution}
\Rdversion{1.1}
\alias{KmeansSolution}
\title{
Initial solution obtained by applying kmeans clustering of atomic strata
}
\description{
In order to speed the convergence towards the optimal solution, an initial one can be 
given as "suggestion" to "optimizeStrata" function.
The function "KmeansSolution" produces this initial solution using the k-means 
algorithm by clustering atomic 
strata on the basis of the values of the means of target variables in them.
Also, if the parameter "nstrata" is not indicated, the optimal number of clusters 
is determined inside each domain, and the overall solution is obtained by concatenating optimal 
clusters obtained in domains.
The result is a dataframe with two columns: the first indicates the clusters, 
the second the domains.
}
\usage{
KmeansSolution(strata, 
               errors, 
               nstrata=NA, 
			         minnumstrat=2,
               maxclusters = NA,
               showPlot=FALSE)
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
  \item{minnumstrat}{
  Minimum number of units to be selected in each stratum. Default is 2.
}
  \item{maxclusters}{
Maximum number of clusters to be considered in the 
execution of kmeans algorithm. If not indicated it will be set equal to the
number of atomic strata divided by 2.
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
  suggestions = solutionKmean
)
}
}
\keyword{ survey }

