\name{KmeansSolution2}
\Rdversion{1.3}
\alias{KmeansSolution2}
\title{
Initial solution obtained by applying kmeans clustering of frame units
}
\description{
This function has to be used only in conjunction with "optimizeStrata2", i.e. in the
case of optimizing with only continuous stratification variables.
The function "KmeansSolution2" has a twofold objective:
- to give indication about a possible best number of final strata (by fixing a 
convenient value for "maxclusters", and leaving NA to "nstrata";
- to give an initial solution fo the optimization step.
If the parameter "nstrata" is not indicated, the optimal number of clusters 
is determined inside each domain, and the overall solution is obtained by concatenating optimal 
clusters obtained in domains.
The result is a dataframe with two columns: the first indicates the clusters, 
the second the domains.
}
\usage{
KmeansSolution2(frame,
                model=NULL,
                errors=cv,
                nstrata = NA,
                minnumstrat =2,
                maxclusters = NA,
                showPlot = TRUE)
}
\arguments{
  \item{frame}{
  The (mandatory) dataframe containing the information related to each unit in the sampling frame.
}
  \item{model}{
  The (optional) dataframe containing the information related to models used to predict
  values of the target variables.
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
Allows to visualise on a plot the different sample sizes for each number of aggregate strata. Default is TRUE.
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

# kmeans clustering
solutionKmean2 <- KmeansSolution2(strata=swissframe,
								model=NULL,
								errors=swisserrors,
								nstrata=NA,
								maxclusters=20,
								showPlot=TRUE)
# preparation of initial suggestion for optimizeStrata2



# optimisation of sampling strata
solution <- optimizeStrata2 ( 
  framesamp = swissframe,
  errors = swisserrors, 
  strata = swissstrata, 
  suggestions = solutionKmean
)
}
}
\keyword{ survey }

