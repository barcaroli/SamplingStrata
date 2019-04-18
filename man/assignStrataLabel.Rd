\name{assignStrataLabel}
\alias{assignStrataLabel}
\title{Function to assign the optimized strata labels}
\description{
  Function to assign the optimized strata labels to new sampling units in the frame
  on the basis of the strata structure obtained by executing the function
  'summaryStrata' after optimizing with 'optimizeStrata2'
}
\usage{
assignStrataLabel(dataset, s)
}
\arguments{
  \item{dataset}{dataset with new sampling units in the frame}
  \item{s}{structure of the strata}
 }
\value{
The same dataset in input with the label of the optimized stratum
}
\examples{
\dontrun{
library(SamplingStrata)
data("swissmunicipalities")
data("errors")
errors$CV1 <- 0.1
errors$CV2 <- 0.1
errors <- errors[rep(row.names(errors),7),]
errors$domainvalue <- c(1:7)
errors
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
swissmunicipalities$domain = 1
frame <- buildFrameDF(swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("Surfacesbois","Surfacescult"),
                      Y = c("Pop020", "Pop2040")
)
solution <- optimizeStrata2 (
  errors, 
  frame,
  nStrata = 5,
  iter = 10, 
  pops = 10, 
  writeFiles = FALSE, 
  showPlot = TRUE, 
  parallel = FALSE) 
strataStructure <- summaryStrata(solution$framenew, solution$aggr_strata)
strataStructure


newset <- assignStrataLabel(solution$framenew,strataStructure)
}
}
\keyword{survey}
