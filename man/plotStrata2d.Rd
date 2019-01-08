\name{plotStrata2d}
\alias{plotStrata2d}
\title{Plot bivariate distibutions in strata}
\description{
  Plots a 2d graph showing obtained strata 
}
\usage{
plotStrata2d (x,outstrata,domain,vars,labels)
}
\arguments{
  \item{x}{the sampling frame}
  \item{outstrata}{the optimized strata}
  \item{domain}{a domain in the frame}
  \item{vars}{vars to appear in x and y axis}
  \item{labels}{labels to appear in x and y axis}
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
p <- plotStrata2d(solution$framenew, 
                  solution$aggr_strata,
                  domain = 2, 
                  vars = c("X1","X2"),
                  labels = c("Surfacesbois","Surfacescult"))
p
}
}
\keyword{survey}
