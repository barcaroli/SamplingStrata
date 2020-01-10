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
prepareSuggestion(kmean = kmean, frame = frame, nstrat = nstrat)
}
\arguments{
  \item{kmean}{
  The result of the execution of function 'KmeansSolution2'.
  }
  \item{frame}{
  The dataframe containing the information related to each unit in the sampling frame.
  }
    \item{nstrat}{
  The vector of number of strata identified as the best for each domain.
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
                      X = c("POPTOT", "HApoly"),
                      Y = c("Surfacesbois", "Airind")
)
ndom <- length(unique(frame$domainvalue))
cv <- as.data.frame(list(
              DOM = rep("DOM1",ndom),
              CV1 = rep(0.1,ndom),
              CV2 = rep(0.1,ndom),
              domainvalue = c(1:ndom)))

# Solution with kmean clustering 
kmean <- KmeansSolution2(frame,model=NULL,errors=cv,nstrata=NA,maxclusters=4)
# Number of strata for each domain
nstrat <- tapply(kmean$suggestions,
                 kmean$domainvalue,
                 FUN=function(x) length(unique(x)))
# Prepare suggestion for optimization step                 
sugg <- prepareSuggestion(kmean = kmean,
                          frame = frame,
                          nstrat = nstrat)
# Optimization
solution <- optimStrata (
  method="continuous",
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = nstrat,
  suggestions = sugg)
}
}
\keyword{ survey }

