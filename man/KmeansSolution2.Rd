\name{KmeansSolution2}
\Rdversion{1.3}
\alias{KmeansSolution2}
\title{
Initial solution obtained by applying kmeans clustering of frame units
}
\description{
This function has to be used only in conjunction with "optimizeStrata2", 
or with "optimStrata" when method = "continuous", i.e. in the
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
                errors,
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
data("swissmunicipalities")
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
swissmunicipalities$dom <- 1
frame <- buildFrameDF(swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("Pop020", "Pop2040"),
                      Y = c("Pop020", "Pop2040")
)
cv <- NULL
cv$DOM <- "DOM1"
cv$CV1 <- 0.1
cv$CV2 <- 0.1
cv <- as.data.frame(cv)
cv <- cv[rep(row.names(cv),7),]
cv$domainvalue <- c(1:7)
cv

# Solution with kmean clustering 
kmean <- KmeansSolution2(frame,model=NULL,errors=cv,nstrata=NA,maxclusters=4)
# number of strata to be obtained in each domain in final solution	
nstrat <- tapply(solutionKmean$suggestions,
                 solutionKmean$domainvalue,
                 FUN=function(x) length(unique(x)))
# Prepare suggestion for optimization
sugg <- prepareSuggestion(kmean,frame,nstrat)
# Optimization
# Attention: number of strata must be the same than in kmeans
solution <- optimStrata (
  method = "continuous",
  cv, 
  framesamp=frame,
  iter = 50,
  pops = 20,
  nStrata = nstrat,
  suggestions = sugg,
  writeFiles = FALSE,
  showPlot = FALSE,
  parallel = FALSE
)
}
}
\keyword{ survey }

