\name{summaryStrata}
\alias{summaryStrata}
\title{Information on strata structure}
\description{
  Information on strata (population, allocation, sampling rate 
  and X variables ranges)
}
\usage{
summaryStrata(x,outstrata,progress,writeFiles)
}
\arguments{
  \item{x}{the sampling frame}
  \item{outstrata}{the optimized strata}
  \item{progress}{progress bar}
  \item{writeFiles}{csv output of the strata structure}
 }
\value{
A formatted output containing information on the strata in the given domain
}
\examples{
\dontrun{
library(SamplingStrata)
data("swissmunicipalities")
swissmunicipalities = swissmunicipalities[swissmunicipalities$REG==1,]
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
            framecens = NULL, 
            strcens = FALSE, 
            alldomains = TRUE, 
            dom = NULL, 
            nStrata = 5, 
            minnumstr = 2, 
            iter = 50, 
            pops = 20, 
            mut_chance = NA, 
            elitism_rate = 0.2, 
            highvalue = 1e+08, 
            suggestions = NULL, 
            realAllocation = TRUE, 
            writeFiles = FALSE, 
            showPlot = TRUE, 
            parallel = FALSE) 
infoStrata <- summaryStrata(solution$framenew, 
<<<<<<< HEAD
                            solution$aggr_strata)
infoStrata
=======
                            solution$aggr_strata,
                            progress = TRUE,
                            writeFiles = FALSE)
infostrata
>>>>>>> 28a8b56e6245e722be8611a495db9b410336ecef
}
}
\keyword{survey}
