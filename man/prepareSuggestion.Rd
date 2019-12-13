\name{prepareSuggestion}
\Rdversion{1.3}
\alias{prepareSuggestion}
\title{
Prepare suggestions for optimization with method = "continuous" or "spatial"
}
\description{
This function has to be used only in conjunction with "KmeansSolution2" or 
with "KmeansSolutionSpatial", 
i.e. in the case of optimizing with only continuous stratification variables.
This function prepares the suggestion for the optimization function in case
of continuous variables (i.e. with with "optimStrata" 
when method = "continuous" or method = "spatial").
}
\usage{
prepareSuggestion(kmean,frame)
}
\arguments{
  \item{kmean}{
  The result of the execution of function 'KmeansSolution2'.
  }
  \item{frame}{
  The (mandatory) dataframe containing the information related to each unit in the sampling frame.
  }
}

\value{
A dataframe containing the suggestions
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
# Prepare suggestion for optimization
sugg <- prepareSuggestion(kmean,frame)
# Optimization
# Attention: number of strata must be the same than in kmeans
solution <- optimizeStrata2 (
  cv, 
  framesamp=frame,
  iter = 50,
  pops = 20,
  nStrata = 4,
  suggestions = sugg,
  writeFiles = FALSE,
  showPlot = FALSE,
  parallel = FALSE
)
}
}
\keyword{ survey }

