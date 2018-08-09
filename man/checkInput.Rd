\name{checkInput}
\Rdversion{1.3}
\alias{checkInput}
\title{
Checks the inputs to the package: dataframes "errors", "strata" and "sampling frame"
}
\description{
This functions checks the internal structure of the different input dataframes ("errors", "strata" and "sampling frame"), 
and also the correctness of the relationships among them.
 
}
\usage{
checkInput(errors=NULL, strata=NULL, sampframe=NULL)
}
\arguments{
  \item{errors}{
  Dataframe containing the precision levels expressed in terms of 
  maximum acceptable coefficients of variation that 
  estimates of target variables Y's of the survey must comply. 
}
  \item{strata}{
  Dataframe containing the information related to strata.
}
  \item{sampframe}{
  Dataframe containing the information related to all the units belonging to the population of interest.
}
}

\author{
Giulio Barcaroli
}


\examples{
library(SamplingStrata)
data(swisserrors)
data(swissstrata)
data(swissframe)
checkInput(swisserrors,swissstrata,swissframe)
checkInput(strata=swissstrata,sampframe=swissframe)
checkInput(strata=swissstrata)
}
\keyword{ survey }

